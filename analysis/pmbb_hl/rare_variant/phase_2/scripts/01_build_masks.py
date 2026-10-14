#!/usr/bin/env python3
"""
Phase 2 -- build the burden masks for one chromosome, from the release annotation.

Run as a 22-task array by 01_build_masks.bsub; 02_merge_masks.py assembles the output.
The first pass did this serially over 106 GB of VEP; splitting it by chromosome turns
hours into minutes and lets a single bad chromosome fail loudly instead of poisoning a
whole run.

ONE BUILD SERVES BOTH PHENOTYPES. A mask says which variants belong to which gene. That
depends on the annotation and on the genotype set, and both arms read the same
`chr*_deduplicated` files, so the masks are identical for the broad and the bilateral
cohorts. Building them twice would only create a chance for them to differ.

THE THREE MASKS

  pLOF      exact term match on {frameshift, stop_gained, start_lost, stop_lost,
            splice_acceptor, splice_donor, transcript_ablation}
            OR any other splice annotation WITH SpliceAI >= 0.2
  pDM       am_class in {pathogenic, likely_pathogenic} OR max REVEL >= 0.5
  pLOF_pDM  the union, pLOF taking precedence in the annotation column

all three AND BIOTYPE == protein_coding.

`ALL` is not built. It was dropped from the analysis in the first pass -- an unfiltered
burden test answers no question this project asks -- so building it would leave a mask
on disk that nothing tests. Per-gene variant counts are emitted separately for QC, which
is the only thing `ALL` was being used for.

THREE THINGS THIS FIXES, each a defect the first pass established in the original run:

  terms matched EXACTLY, not as substrings. The original matched "splice" anywhere in
  the consequence string, which admitted three IMPACT=LOW terms unconditionally and put
  roughly 64% of the pLOF mask outside its own documented definition. The SpliceAI >= 0.2
  gate its analysis plan specifies was never applied.

  REVEL read as the per-transcript LIST VEP writes, maximum taken. The original coerced
  multi-valued entries to NaN, which silently dropped roughly three quarters of pDM.
  Found by Nikki Palmiero on 2026-10-05, after this replication had missed it.

  BIOTYPE == protein_coding. The original's gene set carries 651 antisense, 88 LINC and
  14 microRNA genes. `pLOF` means predicted loss of function, and a gene with no protein
  has no function to lose. Its top burden gene in two of three cohorts was TMC3-AS1, an
  antisense lncRNA.

Terms are matched against the consequence list and NOT against VEP's IMPACT column,
which is per row: a row reading `stop_gained,splice_polypyrimidine_tract_variant` is
HIGH as a whole, so filtering on IMPACT would readmit exactly what is being excluded.
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
REL = Path("/static/PMBB/PMBB-Release-2026-4.0/Exome/vep_annotations")
PLINK = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena"
             "/rarevariant_geneburden/plink_deduplicated")
OUT = HERE / "data/by_chrom"

REAL_LOF = {"stop_gained", "frameshift_variant", "splice_acceptor_variant",
            "splice_donor_variant", "start_lost", "stop_lost",
            "transcript_ablation"}
AM_DM = {"pathogenic", "likely_pathogenic"}
SPLICEAI_GATE = 0.2
REVEL_GATE = 0.5


def max_revel(s):
    best = None
    for x in s.split(","):
        if x in (".", "", "-"):
            continue
        try:
            v = float(x)
        except ValueError:
            continue
        best = v if best is None else max(best, v)
    return best


def max_spliceai(s):
    best = 0.0
    for pred in s.split(","):
        q = pred.split("|")
        if len(q) >= 5:
            for x in q[1:5]:
                try:
                    best = max(best, float(x))
                except ValueError:
                    pass
    return best


def main(chrom):
    bim = PLINK / f"chr{chrom}_deduplicated.bim"
    vep = REL / f"PMBB-Release-2026-4.0_genetic_exome.vep_annotations.chr{chrom}.tsv"
    for f in (bim, vep):
        if not f.exists():
            sys.exit(f"REFUSING: {f} not found")

    # variant ids present in the cohort genotypes, and the id form the mask writes
    bim_ids = {}
    with bim.open() as fh:
        for line in fh:
            p = line.split()[1].split(":")
            if len(p) == 4:
                bim_ids[f"chr{p[0]}_{p[1]}_{p[2]}_{p[3]}"] = ":".join(p)
    if not bim_ids:
        sys.exit(f"REFUSING: no parseable variant ids in {bim}")

    stats = defaultdict(int)
    per_var = defaultdict(lambda: {"plof": False, "pdm": False})

    with vep.open() as fh:
        h = fh.readline().rstrip("\n").split("\t")
        need = ["#Uploaded_variation", "Consequence", "am_class", "REVEL_score",
                "SpliceAI_pred", "SYMBOL", "BIOTYPE"]
        missing = [c for c in need if c not in h]
        if missing:
            sys.exit(f"REFUSING: chr{chrom} VEP lacks columns {missing}")
        iu, ic, iam, irv, isp, isy, ibt = (h.index(c) for c in need)

        for line in fh:
            f = line.rstrip("\n").split("\t")
            vid = f[iu]
            if vid not in bim_ids:
                continue
            stats["rows_in_cohort"] += 1
            if f[ibt] != "protein_coding":
                stats["dropped_not_protein_coding"] += 1
                continue
            sym = f[isy]
            if not sym or sym == "-":
                stats["dropped_no_symbol"] += 1
                continue
            r = per_var[(vid, sym)]
            terms = set(f[ic].split(","))
            if terms & REAL_LOF:
                r["plof"] = True
                stats["plof_by_exact_term"] += 1
            elif any("splice" in t for t in terms):
                if max_spliceai(f[isp]) >= SPLICEAI_GATE:
                    r["plof"] = True
                    stats["plof_by_spliceai_gate"] += 1
                else:
                    stats["splice_rejected_by_gate"] += 1
            rev = max_revel(f[irv])
            if f[iam] in AM_DM:
                r["pdm"] = True
                stats["pdm_by_alphamissense"] += 1
            elif rev is not None and rev >= REVEL_GATE:
                r["pdm"] = True
                stats["pdm_by_revel"] += 1

    masks = {m: defaultdict(list) for m in ("pLOF", "pDM", "pLOF_pDM")}
    counts = defaultdict(int)
    for (vid, sym), r in per_var.items():
        mid = bim_ids[vid]
        counts[sym] += 1
        if r["plof"]:
            masks["pLOF"][sym].append((mid, "pLOF"))
        if r["pdm"]:
            masks["pDM"][sym].append((mid, "pDM"))
        if r["plof"] or r["pdm"]:
            masks["pLOF_pDM"][sym].append((mid, "pLOF" if r["plof"] else "pDM"))

    OUT.mkdir(parents=True, exist_ok=True)
    out = {"chrom": chrom, "stats": dict(stats), "masks": {}}
    for m, genes in masks.items():
        path = OUT / f"{m}.chr{chrom}.txt"
        n = 0
        with path.open("w") as fh:
            for g in sorted(genes):
                items = sorted(set(genes[g]))
                fh.write("\t".join([g, "var"] + [v for v, _ in items]) + "\n")
                fh.write("\t".join([g, "anno"] + [a for _, a in items]) + "\n")
                n += len(items)
        out["masks"][m] = {"genes": len(genes), "entries": n}
    # per-gene variant totals, the QC the ALL mask used to provide
    with (OUT / f"gene_variant_counts.chr{chrom}.tsv").open("w") as fh:
        for g in sorted(counts):
            fh.write(f"{g}\t{counts[g]}\n")
    (OUT / f"manifest.chr{chrom}.json").write_text(json.dumps(out, indent=2) + "\n")

    print(f"chr{chrom}: {stats['rows_in_cohort']:,} annotation rows in cohort, "
          f"{len(counts):,} protein-coding genes")
    for m, v in out["masks"].items():
        print(f"    {m:<10} {v['genes']:>6,} genes  {v['entries']:>8,} entries")


if __name__ == "__main__":
    main(int(sys.argv[1]))
