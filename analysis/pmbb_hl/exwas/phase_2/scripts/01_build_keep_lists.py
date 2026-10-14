#!/usr/bin/env python3
"""
Phase 2, step 1 -- plink2 --keep lists, one per cohort.

The MAF census has to be computed on exactly the samples SAIGE will test, because
MAF is a property of the sample set, not of the variant. A variant at 0.9% in the
combined cohort can sit above 1% in AFR alone and be tested there and not here.

Reads the IID column of the inherited covariate files (premise P1) and writes a
two-column FID/IID file for plink2. FID == IID in these genotypes, verified below.

Refuses rather than proceeding if any cohort sample is absent from the .fam --
a silent --keep miss would shrink the cohort and nothing downstream would notice.
"""
import sys
from pathlib import Path

BURDEN = Path("/project/hall/analysis/hearing-loss-genomics/analysis/pmbb_hl/rare_variant")
COV = BURDEN / "phase_3/results/covariates"
FAM = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena/"
           "rarevariant_geneburden/plink_deduplicated/chr21_deduplicated.fam")
OUT = Path(__file__).resolve().parents[1] / "data"

# cohort -> (covariate filename, expected sample count from PREMISES P1)
# Two phenotypes now, not one. MAF is a property of the sample set, so the census must
# be computed on each cohort SAIGE will actually test -- and the bilateral arm is a
# different sample set, not merely a relabelling of the broad one.
COHORTS = {
    ("broad", "combined"):     ("broad_combined.txt",     57_498),
    ("broad", "EUR"):          ("broad_EUR.txt",          42_779),
    ("broad", "AFR"):          ("broad_AFR.txt",          11_334),
    ("bilateral", "combined"): ("bilateral_combined.txt", 53_910),
    ("bilateral", "EUR"):      ("bilateral_EUR.txt",      40_143),
    ("bilateral", "AFR"):      ("bilateral_AFR.txt",      10_539),
}

def read_fam(path):
    fid_of = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        f = line.split()
        fid_of[f[1]] = f[0]
    return fid_of


def main():
    if not FAM.exists():
        sys.exit(f"REFUSING: .fam not found: {FAM}")
    fid_of = read_fam(FAM)
    print(f".fam: {len(fid_of):,} samples  ({FAM.name})")

    mismatched_fid = [i for i, f in fid_of.items() if f != i]
    if mismatched_fid:
        sys.exit(f"REFUSING: {len(mismatched_fid)} samples have FID != IID; "
                 f"the two-column assumption below is wrong. First: {mismatched_fid[:3]}")
    print("FID == IID for every sample -- two-column --keep is safe")

    OUT.mkdir(parents=True, exist_ok=True)
    problems = []
    for (pheno, cohort), (fname, expected) in COHORTS.items():
        tag = f"{pheno}_{cohort}"
        src = COV / fname
        if not src.exists():
            problems.append(f"{tag}: covariate file missing: {src}")
            continue

        lines = src.read_text().splitlines()
        header = lines[0].split("\t")
        if header[0] != "IID":
            problems.append(f"{tag}: first column is {header[0]!r}, expected 'IID'")
            continue

        iids = [ln.split("\t")[0] for ln in lines[1:] if ln.strip()]
        if len(iids) != len(set(iids)):
            problems.append(f"{tag}: {len(iids) - len(set(iids))} duplicate IIDs")
            continue
        if len(iids) != expected:
            problems.append(f"{tag}: {len(iids):,} samples, PREMISES P1 says {expected:,}")
            continue

        absent = [i for i in iids if i not in fid_of]
        if absent:
            problems.append(f"{tag}: {len(absent)} samples absent from .fam, "
                            f"first: {absent[:3]}")
            continue

        dest = OUT / f"{tag}.keep"
        dest.write_text("".join(f"{i}\t{i}\n" for i in iids))
        print(f"  {tag:<22} {len(iids):>7,} samples -> {dest.name}")

    if problems:
        print("\nREFUSING -- nothing written for:", file=sys.stderr)
        for p in problems:
            print(f"  {p}", file=sys.stderr)
        sys.exit(1)
    print("\nall six cohorts written, every sample present in the genotypes")


if __name__ == "__main__":
    main()
