#!/usr/bin/env python3
"""
Power control, step 1 -- five random size-matched case draws per cohort.

Takes the BROAD arm's covariate files and writes, for each draw, a copy in which a
random subset of cases keeps PHENO=1 and the rest are DROPPED -- not relabelled as
controls. Dropping matches what the bilateral restriction did: the 3,588 excluded
participants have hearing loss and are not controls in the real arm either.

Seed is fixed and recorded. Re-running reproduces the same five draws.
"""
import json
import random
from pathlib import Path

BASE = Path(__file__).resolve().parents[3]
BROAD_COV = (BASE.parent / "gene_burden_elena_replication/phase_3/results/covariates")
OUT = Path(__file__).resolve().parents[1] / "data"
SEED = 20261009
N_DRAWS = 5

# cohort -> (broad covariate file, cases to keep = the bilateral arm's count)
PLAN = {
    "combined": ("corrected_combined_covariates.txt", 3_164),
    "EUR": ("corrected_EUR_covariates.txt", 2_547),
    "AFR": ("corrected_AFR_covariates.txt", 490),
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    manifest = {"seed": SEED, "n_draws": N_DRAWS, "draws": {}}

    for cohort, (fname, keep_n) in PLAN.items():
        src = BROAD_COV / fname
        lines = src.read_text().splitlines()
        header, body = lines[0], [l for l in lines[1:] if l.strip()]
        cases = [l for l in body if l.split("\t")[1] == "1"]
        ctrls = [l for l in body if l.split("\t")[1] != "1"]

        if len(cases) <= keep_n:
            raise SystemExit(f"REFUSING: {cohort} has {len(cases)} cases, "
                             f"cannot draw {keep_n} as a strict subset")

        manifest["draws"][cohort] = {
            "source": str(src.relative_to(BASE.parent.parent)),
            "broad_cases": len(cases), "controls": len(ctrls),
            "cases_kept": keep_n, "cases_dropped": len(cases) - keep_n,
            "files": [],
        }
        for d in range(1, N_DRAWS + 1):
            drawn = rng.sample(cases, keep_n)
            dest = OUT / f"{cohort}_draw{d}.txt"
            dest.write_text("\n".join([header] + drawn + ctrls) + "\n")
            n = len(drawn) + len(ctrls)
            manifest["draws"][cohort]["files"].append(
                {"draw": d, "file": dest.name, "samples": n, "cases": keep_n})
            print(f"  {cohort:<9} draw {d}  {keep_n:>6,} cases + "
                  f"{len(ctrls):,} controls = {n:,}  -> {dest.name}")

    (Path(__file__).resolve().parents[1] / "results/01_draws.json").write_text(
        json.dumps(manifest, indent=2) + "\n")
    print(f"\nseed {SEED}, {N_DRAWS} draws per cohort, manifest written")


if __name__ == "__main__":
    main()
