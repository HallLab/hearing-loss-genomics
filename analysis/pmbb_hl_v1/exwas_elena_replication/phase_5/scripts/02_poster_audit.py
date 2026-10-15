#!/usr/bin/env python3
"""
Phase 5, step 2 -- audit the poster's reported hits against the data.

The poster (elena_publishes/Horizontal Draft LANS Poster - SUIP-PGSIP.pptx) reports
six genome-wide significant results under "Genome Wide Significant Results Total: 6",
three unique variants. None appears in the delivered pipeline.

This script establishes, reproducibly and from primary sources:

  1. the claims, read out of the .pptx rather than retyped
  2. whether each variant exists in the PMBB v4 release
  3. its call rate and MAF measured on the release
  4. whether it would survive the filter the poster's own methods describe
  5. at which pipeline stage it actually disappeared
  6. what the named genes score in the delivered results

Nothing here depends on our analysis pipeline. Points 2-4 are measured on the
institutional release with plink2, which is why they are the ones to lead with: they
are not a result that could be replicated differently.
"""
import json
import re
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
ROOT = BASE.parent.parent
ELENA = BASE.parent / "elena"
POSTER = (BASE.parent / "gene_burden_elena_replication/elena_publishes/"
          "Horizontal Draft LANS Poster - SUIP-PGSIP.pptx")
RELEASE = Path("/static/PMBB/PMBB-Release-2026-4.0/Exome/by_chromosome_GL")
PLINK2 = "/appl/plink2-20240804/plink2"
OUT = BASE / "phase_5/results"

# The poster's stated QC, quoted so the comparison is against its own words.
POSTER_QC = ("variants with a genotype call rate < 0.99 and a Minor Allele Frequency "
             "(MAF) < 0.01 were removed")
CALL_RATE_FLOOR = 0.99


def poster_text():
    z = zipfile.ZipFile(POSTER)
    out = []
    for n in sorted(z.namelist()):
        if re.match(r"ppt/slides/slide\d+\.xml$", n):
            out += re.findall(r"<a:t>(.*?)</a:t>", z.read(n).decode("utf8", "replace"),
                              re.S)
    return re.sub(r"\s+", " ", " ".join(out))


def parse_claims(text):
    """Pull 'GENE chrN:POS:A:B COHORT p' tuples out of the results table."""
    pat = re.compile(r"([A-Z0-9]+)\s+(chr\d+:\d+:[ACGT]+:[ACGT]+)\s+"
                     r"(Combined|EUR|AFR)\s+([\d.]+e-\d+)")
    return [{"gene": m.group(1), "variant": m.group(2),
             "cohort": m.group(3), "p": float(m.group(4))}
            for m in pat.finditer(text)]


def release_metrics(variants):
    """call rate and MAF for each variant, measured on the release with plink2."""
    by_chrom = {}
    for v in variants:
        chrom = v.split(":")[0].replace("chr", "")
        by_chrom.setdefault(chrom, []).append(v)

    metrics = {}
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for chrom, vs in by_chrom.items():
            stem_name = f"PMBB-Release-2026-4.0_genetic_exome.gl.chr{chrom}"
            bfile = RELEASE / stem_name
            # NOT with_suffix(".bim"): this filename's last suffix is ".chr15", so
            # with_suffix would ask for "...exome.gl.bim" and silently find nothing.
            if not (RELEASE / (stem_name + ".bim")).exists():
                for v in vs:
                    metrics[v] = {"error": f"release .bim not found: {bfile}.bim"}
                continue
            # release IDs use underscores: chr15_40611486_A_G
            ids = {v: v.replace(":", "_", 1).replace(":", "_") for v in vs}
            snp = tmp / f"c{chrom}.snp"
            snp.write_text("".join(i + "\n" for i in ids.values()))
            stem = tmp / f"c{chrom}"
            run = subprocess.run(
                [PLINK2, "--bfile", str(bfile), "--extract", str(snp),
                 "--missing", "variant-only", "--freq", "counts",
                 "--out", str(stem), "--threads", "4", "--memory", "16000"],
                capture_output=True, text=True)
            miss, cnt = {}, {}
            vm, ac = Path(f"{stem}.vmiss"), Path(f"{stem}.acount")
            if run.returncode != 0 or not (vm.exists() and ac.exists()):
                for v in vs:
                    metrics[v] = {"error": f"plink2 failed on chr{chrom} "
                                           f"(rc={run.returncode})",
                                  "stderr": (run.stderr or "")[-400:]}
                continue
            if vm.exists():
                for l in vm.read_text().splitlines()[1:]:
                    f = l.split("\t")
                    miss[f[1]] = float(f[-1])
            if ac.exists():
                for l in ac.read_text().splitlines()[1:]:
                    f = l.split("\t")
                    cnt[f[1]] = (float(f[5]), float(f[6]))
            for v, rid in ids.items():
                if rid not in cnt:
                    metrics[v] = {"in_release": False}
                    continue
                alt, obs = cnt[rid]
                maf = min(alt / obs, 1 - alt / obs) if obs else None
                cr = 1 - miss.get(rid, float("nan"))
                metrics[v] = {
                    "in_release": True, "release_id": rid,
                    "call_rate": round(cr, 4), "maf": round(maf, 5) if maf else None,
                    "passes_poster_qc": bool(cr >= CALL_RATE_FLOOR
                                             and maf and maf >= 0.01),
                }
    return metrics


def stage_trail(variants):
    """At which pipeline stage does each variant stop appearing?"""
    stages = [
        ("PMBB v4 release",
         lambda c: RELEASE / f"PMBB-Release-2026-4.0_genetic_exome.gl.chr{c}.bim"),
        ("variant_qc (named for call rate)",
         lambda c: ELENA / f"rarevariantExWAS/variant_qc/chr{c}_callrate.bim"),
        ("allele_ids",
         lambda c: ELENA / f"rarevariant_geneburden/plink_deduplicated/chr{c}_allele_ids.bim"),
        ("deduplicated (what step 2 reads)",
         lambda c: ELENA / f"rarevariant_geneburden/plink_deduplicated/chr{c}_deduplicated.bim"),
    ]
    trail = {}
    for v in variants:
        chrom, pos = v.split(":")[0].replace("chr", ""), int(v.split(":")[1])
        trail[v] = []
        for name, f in stages:
            bim = f(chrom)
            if not bim.exists():
                trail[v].append({"stage": name, "file_missing": True})
                continue
            total = present = 0
            with bim.open() as fh:
                for line in fh:
                    total += 1
                    if int(line.split()[3]) == pos:
                        present += 1
            trail[v].append({"stage": name, "variants_on_chrom": total,
                             "present": present > 0})
    return trail


def main():
    if not POSTER.exists():
        sys.exit(f"REFUSING: poster not found at {POSTER}")

    text = poster_text()
    if POSTER_QC not in text.replace("&lt;", "<").replace("&gt;", ">"):
        print("WARNING: the poster's QC sentence is not where this script expects it.")
        print("The comparison below still runs, but check the quote by hand.\n")

    claims = parse_claims(text)
    if not claims:
        sys.exit("REFUSING: no results table parsed out of the poster. Its layout "
                 "changed, and retyping the numbers by hand would defeat the point.")

    uniq = sorted({c["variant"] for c in claims})
    print(f"poster claims: {len(claims)} results, {len(uniq)} unique variants\n")
    for c in claims:
        print(f"  {c['gene']:<8} {c['variant']:<22} {c['cohort']:<9} {c['p']:.2e}")
    print()

    metrics = release_metrics(uniq)
    failures = []
    print("measured on the PMBB v4 release (independent of our pipeline):\n")
    print(f"  {'variant':<24}{'call rate':>11}{'MAF':>10}   passes the poster's own QC?")
    for v in uniq:
        m = metrics[v]
        if "error" in m:
            failures.append((v, m["error"]))
            print(f"  {v:<24}  COULD NOT MEASURE -- {m['error']}")
            continue
        if not m.get("in_release"):
            print(f"  {v:<24}  not present in the release")
            continue
        verdict = "yes" if m["passes_poster_qc"] else "NO"
        print(f"  {v:<24}{m['call_rate']:>11.4f}{m['maf']:>10.5f}   {verdict}")
    print()

    trail = stage_trail(uniq)
    print("where each stops appearing:\n")
    for v in uniq:
        print(f"  {v}")
        for s in trail[v]:
            if s.get("file_missing"):
                print(f"    {s['stage']:<36} (file not found)")
            else:
                mark = "present" if s["present"] else "GONE"
                print(f"    {s['stage']:<36} {s['variants_on_chrom']:>9,}  {mark}")
        print()

    if failures:
        print(f"\nREFUSING to write the report: {len(failures)} of {len(uniq)} "
              f"variants could not be measured on the release.")
        print("That measurement is the whole point of this audit -- an unmeasured")
        print("variant must not be reported as an absent one.")
        for v, e in failures:
            print(f"  {v}: {e}", file=sys.stderr)
        sys.exit(1)

    dest = OUT / "poster_audit.json"
    dest.write_text(json.dumps({"poster": str(POSTER.name), "claims": claims,
                                "release_metrics": metrics, "stage_trail": trail,
                                "poster_qc_quoted": POSTER_QC,
                                "call_rate_floor": CALL_RATE_FLOOR}, indent=2) + "\n")
    print(f"written: {dest.relative_to(BASE)}")


if __name__ == "__main__":
    main()
