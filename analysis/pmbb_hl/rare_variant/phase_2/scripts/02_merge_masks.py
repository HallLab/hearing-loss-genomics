#!/usr/bin/env python3
"""
Phase 2 -- merge the per-chromosome masks, with the checks worth having.

Refuses rather than producing a partial mask. A mask that is quietly short does not
fail downstream: SAIGE happily tests whatever genes it is given, and the result looks
like a finished analysis with fewer genes in it.

WHAT IS CHECKED

  all 22 manifests present, and LSF asked rather than mtimes trusted
  every mask file parses as alternating var/anno lines with equal lengths
  entry counts match what each chromosome's manifest recorded
  no gene appears on two chromosomes -- a symbol on two chromosomes means either an
    annotation problem or a name collision, and silently merging the two would pool
    unrelated variants into one burden test
  pLOF_pDM is exactly the union of pLOF and pDM, per gene

The last one is cheap and catches a whole class of mistake: if the union were built
wrong, every pLOF_pDM test would be subtly incorrect and nothing else would notice.
"""
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
BY_CHROM = HERE / "data/by_chrom"
OUT = HERE / "results/masks"
MASKS = ["pLOF", "pDM", "pLOF_pDM"]
CHROMS = range(1, 23)


def in_flight():
    try:
        out = subprocess.run(["bjobs", "-J", "masks"], capture_output=True,
                             text=True, timeout=60).stdout
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    return [l for l in out.splitlines()[1:] if l.strip()]


def read_mask(path):
    """{gene: [(variant, annotation), ...]}, refusing on any structural surprise."""
    genes = {}
    lines = path.read_text().splitlines()
    if len(lines) % 2:
        raise ValueError(f"{path.name}: odd number of lines, so a var/anno pair is "
                         f"incomplete")
    for i in range(0, len(lines), 2):
        v = lines[i].split("\t")
        a = lines[i + 1].split("\t")
        if v[1] != "var" or a[1] != "anno":
            raise ValueError(f"{path.name} line {i+1}: expected var/anno, got "
                             f"{v[1]!r}/{a[1]!r}")
        if v[0] != a[0]:
            raise ValueError(f"{path.name} line {i+1}: gene {v[0]!r} != {a[0]!r}")
        if len(v) != len(a):
            raise ValueError(f"{path.name} gene {v[0]}: {len(v)-2} variants but "
                             f"{len(a)-2} annotations")
        genes[v[0]] = list(zip(v[2:], a[2:]))
    return genes


def main():
    live = in_flight()
    if live:
        sys.exit(f"REFUSING: {len(live)} mask tasks still in flight.")
    if live is None:
        print("WARNING: could not query LSF; relying on the per-file checks\n")

    missing = [c for c in CHROMS
               if not (BY_CHROM / f"manifest.chr{c}.json").exists()]
    if missing:
        sys.exit(f"REFUSING: manifests absent for chromosomes {missing}")

    merged = {m: {} for m in MASKS}
    gene_chrom = {}
    problems = []
    stats = defaultdict(int)

    for c in CHROMS:
        man = json.loads((BY_CHROM / f"manifest.chr{c}.json").read_text())
        for k, v in man["stats"].items():
            stats[k] += v
        for m in MASKS:
            path = BY_CHROM / f"{m}.chr{c}.txt"
            if not path.exists():
                problems.append(f"{m} chr{c}: file absent")
                continue
            try:
                genes = read_mask(path)
            except ValueError as e:
                problems.append(str(e))
                continue
            entries = sum(len(v) for v in genes.values())
            if entries != man["masks"][m]["entries"]:
                problems.append(f"{m} chr{c}: {entries} entries on disk, manifest "
                                f"says {man['masks'][m]['entries']}")
                continue
            for g, items in genes.items():
                prev = gene_chrom.setdefault(g, c)
                if prev != c:
                    problems.append(f"gene {g} appears on chr{prev} and chr{c}")
                    continue
                merged[m].setdefault(g, []).extend(items)

    if problems:
        print(f"REFUSING -- {len(problems)} structural problems:", file=sys.stderr)
        for p in problems[:20]:
            print(f"  {p}", file=sys.stderr)
        sys.exit(1)

    # pLOF_pDM must be exactly the union
    bad_union = []
    for g, items in merged["pLOF_pDM"].items():
        want = {v for v, _ in merged["pLOF"].get(g, [])} | \
               {v for v, _ in merged["pDM"].get(g, [])}
        got = {v for v, _ in items}
        if want != got:
            bad_union.append(f"{g}: union has {len(want)}, pLOF_pDM has {len(got)}")
    if bad_union:
        print(f"REFUSING: pLOF_pDM is not the union of pLOF and pDM for "
              f"{len(bad_union)} genes:", file=sys.stderr)
        for b in bad_union[:10]:
            print(f"  {b}", file=sys.stderr)
        sys.exit(1)

    OUT.mkdir(parents=True, exist_ok=True)
    manifest = {"source": "release VEP, rebuilt per chromosome",
                "annotation_rows_in_cohort": stats["rows_in_cohort"],
                "fixes": {
                    "dropped_not_protein_coding": stats["dropped_not_protein_coding"],
                    "splice_rejected_by_spliceai_gate":
                        stats["splice_rejected_by_gate"],
                    "pdm_from_revel": stats["pdm_by_revel"],
                    "pdm_from_alphamissense": stats["pdm_by_alphamissense"],
                    "plof_by_exact_term": stats["plof_by_exact_term"],
                    "plof_by_spliceai_gate": stats["plof_by_spliceai_gate"],
                },
                "masks": {}}

    print(f"{stats['rows_in_cohort']:,} annotation rows inside the cohort genotypes\n")
    for m in MASKS:
        path = OUT / f"{m}.txt"
        n = 0
        with path.open("w") as fh:
            for g in sorted(merged[m]):
                items = sorted(set(merged[m][g]))
                fh.write("\t".join([g, "var"] + [v for v, _ in items]) + "\n")
                fh.write("\t".join([g, "anno"] + [a for _, a in items]) + "\n")
                n += len(items)
        manifest["masks"][m] = {"genes": len(merged[m]), "entries": n}
        print(f"  {m:<10} {len(merged[m]):>6,} genes   {n:>9,} entries")

    dest = HERE / "results/02_masks_manifest.json"
    dest.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nwritten: results/masks/, {dest.name}")
    print(f"\nfixes, genome-wide:")
    for k, v in manifest["fixes"].items():
        print(f"  {k:<36} {v:>10,}")


if __name__ == "__main__":
    main()
