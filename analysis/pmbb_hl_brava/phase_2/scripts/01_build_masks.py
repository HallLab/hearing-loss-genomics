#!/usr/bin/env python3
"""
Phase 2 (BRaVa) -- build masks for one chromosome under the BRaVa annotation standard.

WHAT THIS IS NOT. It is not a correction of `pmbb_hl/rare_variant/phase_2`. That arm's
masks are defensible and reproduce the first pass byte for byte. This is a SECOND
annotation of the same genotypes, built so that (a) our results can be placed beside a
BRaVa-annotated biobank without the difference being unattributable, and (b) the ClinGen
enrichment can be tested against a mask set that shares only a fifth of its variants.

Rules, from `docs/papers/BRaVa annotation recommendations.pdf`. Categories are in
PRECEDENCE order and each variant receives exactly one:

  1  pLOF        LOFTEE HC
  2  pDM         not (1), and any of
                   a. missense / start_lost / stop_lost / inframe indel
                      AND (REVEL >= 0.773 OR CADD_PHRED >= 28.1)
                   b. SpliceAI DS >= 0.2,  DS = max(DS_AG, DS_AL, DS_DG, DS_DL)
                   c. LOFTEE LC
  4  synonymous  synonymous_variant AND SpliceAI DS < 0.2        -- the CONTROL set

Category 3 (other missense) is not built; nothing here tests it.

FOUR DIFFERENCES FROM `pmbb_hl`, each deliberate and each measured in phase 5:

  pLOF is LOFTEE, not consequence terms. Ours admits any row carrying a LoF term; BRaVa
  admits only what LOFTEE calls high-confidence, which excludes LoF in the last exon,
  rescued splice sites, single-exon genes and NAGNAG acceptors. On chr21 ours calls
  12,181 rows pLOF, of which LOFTEE says HC for 10,555, LC for 1,111, and declines to
  call 515. The LC ones are not discarded here -- BRaVa moves them to pDM.

  pDM uses CADD and a REVEL cut of 0.773. Ours uses AlphaMissense and REVEL >= 0.5.
  AlphaMissense postdates the BRaVa document, so this is not BRaVa rejecting it.

  ONE TRANSCRIPT, not all of them. BRaVa specifies MANE Select, falling back to
  CANONICAL. Ours reads every transcript row for a variant and takes the maximum REVEL
  across them, which is more permissive by construction: a variant benign on its MANE
  transcript can enter our pDM on a minor one.

  A gnomAD FREQUENCY FILTER. BRaVa drops any variant above 1% in any broad continental
  gnomAD population. Ours filters on INTERNAL MAF at test time, which is a different
  cut: a variant rare in PMBB and common in gnomAD passes ours and fails theirs.

THE SYNONYMOUS MASK IS A CONTROL AND IS EXPECTED TO SHOW NOTHING. Synonymous variants
without splice signal should not associate with hearing loss. If they do, the problem is
structural -- something about how cases and controls differ -- and no other check in
this project would catch it. See PREMISES.md P4.
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

# BRaVa 2a -- the protein-altering consequences the REVEL/CADD gate applies to. Narrower
# than "anything missense-ish": the document lists these four kinds and no others.
PROTEIN_ALTERING = {"missense_variant", "start_lost", "stop_lost",
                    "inframe_insertion", "inframe_deletion"}
REVEL_GATE = 0.773          # BRaVa, from biorxiv 2022.03.17.484479
CADD_GATE = 28.1            # same source
SPLICEAI_GATE = 0.2
GNOMAD_MAX_AF = 0.01
# "major continental" populations. gnomAD's exome set is used; ASJ, FIN, MID and
# REMAINING are deliberately excluded -- they are not broad continental groups and
# including them would let a drift-elevated frequency in a small population drop a
# variant that is rare everywhere else.
GNOMAD_POPS = ["gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_EAS_AF",
               "gnomADe_NFE_AF", "gnomADe_SAS_AF"]

NEED = ["#Uploaded_variation", "Consequence", "SYMBOL", "BIOTYPE", "Feature",
        "LoF", "REVEL_score", "CADD_PHRED", "SpliceAI_pred",
        "MANE_SELECT", "CANONICAL"] + GNOMAD_POPS


def num(s):
    """One float, or None. VEP writes '-' and '' for absent."""
    if s in (".", "", "-", None):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def max_field(s):
    """Max over VEP's comma-separated per-transcript list. Used for REVEL and CADD,
    which can carry one value per transcript even on a single row."""
    best = None
    for x in str(s).split(","):
        v = num(x)
        if v is not None:
            best = v if best is None else max(best, v)
    return best


def max_spliceai(s):
    """BRaVa's DS = max(DS_AG, DS_AL, DS_DG, DS_DL).

    THE FIELD OFFSET IS THE TRAP. This release writes

        SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL

    so the four delta scores are fields 1..4. Reading 2..5 picks up DP_AG, a POSITION
    in base pairs -- values like -7, 14, 28 -- which clears any threshold on sight. The
    first run of this script did exactly that and called 60,516 variants spliceogenic on
    chr21 alone, against 2,394 LOFTEE HC. A delta score is a probability, so the range
    check below makes that failure loud instead of plausible."""
    best = 0.0
    for pred in str(s).split(","):
        q = pred.split("|")
        if len(q) >= 5:
            for x in q[1:5]:
                v = num(x)
                if v is None:
                    continue
                if not 0.0 <= v <= 1.0:
                    raise ValueError(
                        f"SpliceAI delta score out of range: {v} in {pred!r}. "
                        f"The field offset is wrong -- a position was read as a score.")
                best = max(best, v)
    return best


def gnomad_max(fields, idx):
    best = 0.0
    for i in idx:
        v = num(fields[i])
        if v is not None:
            best = max(best, v)
    return best


def main(chrom):
    bim = PLINK / f"chr{chrom}_deduplicated.bim"
    vep = REL / f"PMBB-Release-2026-4.0_genetic_exome.vep_annotations.chr{chrom}.tsv"
    for f in (bim, vep):
        if not f.exists():
            sys.exit(f"REFUSING: {f} not found")

    bim_ids = {}
    with bim.open() as fh:
        for line in fh:
            p = line.split()[1].split(":")
            if len(p) == 4:
                bim_ids[f"chr{p[0]}_{p[1]}_{p[2]}_{p[3]}"] = ":".join(p)
    if not bim_ids:
        sys.exit(f"REFUSING: no parseable variant ids in {bim}")

    stats = defaultdict(int)
    # BRaVa wants ONE transcript per variant. Rows are ranked: MANE Select beats
    # CANONICAL beats anything else, and the best row seen wins. Keeping the row rather
    # than merging across rows is the whole difference from pmbb_hl, so it is done
    # explicitly rather than by a groupby that would hide the choice.
    chosen = {}

    with vep.open() as fh:
        h = fh.readline().rstrip("\n").split("\t")
        missing = [c for c in NEED if c not in h]
        if missing:
            sys.exit(f"REFUSING: chr{chrom} VEP lacks columns {missing}")
        I = {c: h.index(c) for c in NEED}
        gidx = [I[c] for c in GNOMAD_POPS]

        for line in fh:
            f = line.rstrip("\n").split("\t")
            vid = f[I["#Uploaded_variation"]]
            if vid not in bim_ids:
                continue
            stats["rows_in_cohort"] += 1
            if f[I["BIOTYPE"]] != "protein_coding":
                stats["dropped_not_protein_coding"] += 1
                continue
            sym = f[I["SYMBOL"]]
            if not sym or sym == "-":
                stats["dropped_no_symbol"] += 1
                continue

            mane = f[I["MANE_SELECT"]] not in ("", "-")
            canon = f[I["CANONICAL"]] == "YES"
            rank = 2 if mane else (1 if canon else 0)
            if rank == 0:
                stats["dropped_not_mane_or_canonical"] += 1
                continue

            prev = chosen.get(vid)
            if prev is not None and prev["rank"] >= rank:
                stats["rows_superseded"] += 1
                continue
            chosen[vid] = {"rank": rank, "sym": sym,
                           "cq": set(f[I["Consequence"]].split(",")),
                           "lof": f[I["LoF"]],
                           "revel": max_field(f[I["REVEL_score"]]),
                           "cadd": max_field(f[I["CADD_PHRED"]]),
                           "ds": max_spliceai(f[I["SpliceAI_pred"]]),
                           "gaf": gnomad_max(f, gidx)}

    masks = {m: defaultdict(list) for m in ("pLOF", "pDM", "pLOF_pDM", "synonymous")}
    for vid, r in chosen.items():
        if r["gaf"] > GNOMAD_MAX_AF:
            stats["dropped_gnomad_common"] += 1
            continue
        mid, sym = bim_ids[vid], r["sym"]

        # precedence: 1 beats 2 beats 4, and each variant gets exactly one
        if r["lof"] == "HC":
            cat = "pLOF"
            stats["cat1_loftee_hc"] += 1
        elif (r["cq"] & PROTEIN_ALTERING
              and ((r["revel"] is not None and r["revel"] >= REVEL_GATE)
                   or (r["cadd"] is not None and r["cadd"] >= CADD_GATE))):
            cat = "pDM"
            stats["cat2a_revel_or_cadd"] += 1
        elif r["ds"] >= SPLICEAI_GATE:
            cat = "pDM"
            stats["cat2b_spliceai"] += 1
        elif r["lof"] == "LC":
            cat = "pDM"
            stats["cat2c_loftee_lc"] += 1
        elif "synonymous_variant" in r["cq"] and r["ds"] < SPLICEAI_GATE:
            cat = "synonymous"
            stats["cat4_synonymous"] += 1
        else:
            stats["uncategorised"] += 1
            continue

        masks[cat][sym].append((mid, cat if cat != "synonymous" else "synonymous"))
        if cat in ("pLOF", "pDM"):
            masks["pLOF_pDM"][sym].append((mid, cat))

    OUT.mkdir(parents=True, exist_ok=True)
    out = {"chrom": chrom, "stats": dict(stats), "masks": {}}
    for m, genes in masks.items():
        path = OUT / f"{m}.chr{chrom}.txt"
        n = 0
        with path.open("w") as fh:
            for g in sorted(genes):
                items = sorted(set(genes[g]))
                fh.write(f"{g}\tvar\t" + "\t".join(i[0] for i in items) + "\n")
                fh.write(f"{g}\tanno\t" + "\t".join(i[1] for i in items) + "\n")
                n += len(items)
        out["masks"][m] = {"genes": len(genes), "entries": n}

    (OUT / f"manifest.chr{chrom}.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"chr{chrom}")
    for m, e in out["masks"].items():
        print(f"  {m:<12} {e['genes']:>6,} genes  {e['entries']:>9,} entries")
    for k in sorted(stats):
        print(f"    {k:<32} {stats[k]:>12,}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: 01_build_masks.py <chrom 1-22>")
    main(sys.argv[1])
