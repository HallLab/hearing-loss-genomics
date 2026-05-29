#!/usr/bin/env python3
"""
znf175_signal_loss_diagnostic.py — Chapter 2 · Phase 2.1 + 2.2

Diagnose the ZNF175 signal-loss between PMBB v1 (~11k, Joe Park's analysis)
and PMBB v2 (~43k, Daniel's failed reproduction).

Inputs:
- data/PMBB_Exome/ZNF175/Joe_analyses/znf175_variants.txt.gz
    Joe's preserved variant list (hg38 coordinates, with PMBB v1 Hom/Het counts)
- data/PMBB_Exome/ZNF175/ZNF175_annot_genes_full.txt.gz
    All ZNF175 variants in PMBB v2 annotation (with v2 Hom/Het counts, REVEL,
    gnomAD, ClinVar, etc.)

Outputs (in analysis/daniel/outputs/phase8_signal_diagnostic/):
- znf175_inventory.tsv      — full variant inventory with all filters applied
- joe_variants_status.tsv   — Joe's original variants with their status in v2
- variant_115_candidates.tsv — variants with ~115 carriers in v2 (Doug's flag)
- summary.md                — written summary of findings
"""

import argparse
import gzip
import sys
from pathlib import Path


PMBB_V2_N = 43_731

# Joe's 10 hg38 positions + their hg19 lifts (from Daniel's runbook lines 213-232)
# Format: hg38_position → (hg19_position, runbook_note)
JOE_HG38_TO_HG19 = {
    "52084689-52084690":  ("51581436-51581437", "yes"),
    "52084742-52084743":  ("51581489-51581490", "yes"),
    "52090597-52090598":  ("51587344-51587345", "no, but 51587346 is SNP stopgain"),
    "52090980-52090981":  ("51587727-51587728", "not in -- multiallelic"),
    "52091066-52091067":  ("51587813-51587814", "yes"),
    "52091163-52091164":  ("51587910-51587911", "not in -- multiallelic"),
    "52091467-52091468":  ("51588214-51588215", "in, is indel"),
    "52091567-52091568":  ("51588314-51588315", "yes"),
    "52091634-52091635":  ("51588381-51588382", "not in -- not in annotation"),
    "52091681-52091682":  ("51588428-51588429", "yes"),
}


def read_joe_variants(path: Path):
    """Read Joe's original PMBB v1 variant list (hg38)."""
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            row = dict(zip(header, parts))
            rows.append(row)
    return header, rows


def read_v2_annotation(path: Path):
    """Read PMBB v2 ZNF175 annotation."""
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            row = dict(zip(header, parts))
            rows.append(row)
    return header, rows


def safe_int(x, default=0):
    try:
        return int(x)
    except (ValueError, TypeError):
        return default


def safe_float(x, default=None):
    try:
        if x in ("", ".", "NA", "nan"):
            return default
        return float(x)
    except (ValueError, TypeError):
        return default


def classify_function(func_refgene, exonic_func):
    """Classify variant function class."""
    plof_classes = {
        "frameshift substitution", "frameshift_substitution",
        "frameshift insertion", "frameshift_insertion",
        "frameshift deletion", "frameshift_deletion",
        "stopgain", "stoploss", "splicing",
    }
    if func_refgene == "splicing":
        return "pLOF"
    if exonic_func in plof_classes:
        return "pLOF"
    if exonic_func == "nonsynonymous SNV":
        return "missense"
    if exonic_func in ("synonymous SNV", "synonymous_SNV"):
        return "synonymous"
    if func_refgene in ("UTR3", "UTR5", "intronic", "intergenic", "downstream", "upstream"):
        return "non-coding"
    return f"other:{exonic_func or func_refgene}"


def v2_maf(het, hom, missing, n_total=PMBB_V2_N):
    """Compute MAF from PMBB v2 counts."""
    called = n_total - missing
    if called <= 0:
        return None
    alt_alleles = het + 2 * hom
    return alt_alleles / (2 * called)


def fmt(x, prec=4):
    if x is None:
        return "—"
    if isinstance(x, float):
        if abs(x) < 1e-3 and x != 0:
            return f"{x:.2e}"
        return f"{x:.{prec}f}"
    return str(x)


def analyze_v2(rows):
    """Build full v2 ZNF175 inventory."""
    out = []
    for r in rows:
        het = safe_int(r.get("HET_REF_ALT_CTS", 0))
        hom = safe_int(r.get("TWO_ALT_GENO_CTS", 0))
        missing = safe_int(r.get("MISSING_CT", 0))
        carriers = het + hom
        maf = v2_maf(het, hom, missing)
        revel = safe_float(r.get("REVEL_score"))
        gnomad = safe_float(r.get("gnomAD_exome_ALL"))
        clnsig = r.get("CLNSIG", "") or ""
        func_class = classify_function(r.get("Func.refGene", ""),
                                       r.get("ExonicFunc.refGene", ""))
        passes_maf_001 = maf is not None and maf < 0.001
        passes_maf_01 = maf is not None and maf < 0.01
        passes_revel_06 = revel is not None and revel > 0.6
        is_burden_eligible = func_class == "pLOF" or (func_class == "missense" and passes_revel_06)
        out.append({
            "id": r.get("ID", ""),
            "pos": int(r.get("Start", 0) or 0),
            "ref": r.get("Ref", ""),
            "alt": r.get("Alt", ""),
            "func_class": func_class,
            "exonic_func": r.get("ExonicFunc.refGene", ""),
            "het": het,
            "hom": hom,
            "missing": missing,
            "carriers": carriers,
            "v2_maf": maf,
            "gnomad_maf": gnomad,
            "revel": revel,
            "clnsig": clnsig,
            "passes_maf_001": passes_maf_001,
            "passes_maf_01": passes_maf_01,
            "passes_revel_06": passes_revel_06,
            "burden_eligible": is_burden_eligible,
        })
    return out


def analyze_joe(rows):
    """Build Joe's PMBB v1 variant list."""
    out = []
    for r in rows:
        het = safe_int(r.get("Het", 0))
        hom = safe_int(r.get("Hom", 0))
        missing = safe_int(r.get("Missing", 0))
        carriers = het + hom
        revel = safe_float(r.get("REVEL"))
        gnomad = safe_float(r.get("gnomAD_exome_ALL"))
        clnsig = r.get("CLNSIG", "") or ""
        out.append({
            "constant_id": r.get("Constant_ID", ""),
            "pos_hg38": int(r.get("Start", 0) or 0),
            "ref": r.get("Ref", ""),
            "alt": r.get("Alt", ""),
            "func": r.get("Func.refGene", ""),
            "exonic_func": r.get("ExonicFunc.refGene", ""),
            "v1_het": het,
            "v1_hom": hom,
            "v1_missing": missing,
            "v1_carriers": carriers,
            "revel": revel,
            "gnomad_maf": gnomad,
            "clnsig": clnsig,
        })
    return out


def write_tsv(rows, path, cols):
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--joe-list", required=True, type=Path,
                    help="data/PMBB_Exome/ZNF175/Joe_analyses/znf175_variants.txt.gz")
    ap.add_argument("--v2-annot", required=True, type=Path,
                    help="data/PMBB_Exome/ZNF175/ZNF175_annot_genes_full.txt.gz")
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== Loading Joe's PMBB v1 variant list ===")
    _, joe_rows = read_joe_variants(args.joe_list)
    print(f"  {len(joe_rows)} variants in Joe's list")

    print(f"\n=== Loading PMBB v2 ZNF175 annotation ===")
    _, v2_rows = read_v2_annotation(args.v2_annot)
    print(f"  {len(v2_rows)} variants in v2 annotation")

    print(f"\n=== Analyzing v2 variants ===")
    v2 = analyze_v2(v2_rows)

    # Filter summaries
    plof = [v for v in v2 if v["func_class"] == "pLOF"]
    missense = [v for v in v2 if v["func_class"] == "missense"]
    plof_maf001 = [v for v in plof if v["passes_maf_001"]]
    missense_burden = [v for v in missense if v["burden_eligible"] and v["passes_maf_001"]]
    print(f"  Total v2 variants: {len(v2)}")
    print(f"  pLOF: {len(plof)}")
    print(f"  Missense: {len(missense)}")
    print(f"  pLOF passing MAF<0.001: {len(plof_maf001)}")
    print(f"  Missense passing MAF<0.001 AND REVEL>0.6: {len(missense_burden)}")

    print(f"\n=== Analyzing Joe's v1 variants ===")
    joe = analyze_joe(joe_rows)
    joe_with_carriers = [j for j in joe if j["v1_carriers"] > 0]
    print(f"  Total Joe variants: {len(joe)}")
    print(f"  Joe variants with v1 carriers > 0: {len(joe_with_carriers)}")

    # ----- Output 1: Full v2 inventory -----
    print(f"\n=== Writing v2 inventory ===")
    cols = ["id", "pos", "ref", "alt", "func_class", "exonic_func",
            "het", "hom", "carriers", "v2_maf", "gnomad_maf", "revel", "clnsig",
            "passes_maf_001", "passes_revel_06", "burden_eligible"]
    write_tsv(sorted(v2, key=lambda r: r["pos"]),
              args.out_dir / "znf175_v2_inventory.tsv", cols)

    # ----- Output 2: Joe's variants with current v2 status -----
    print(f"=== Writing Joe variants ===")
    joe_cols = ["constant_id", "pos_hg38", "ref", "alt", "func", "exonic_func",
                "v1_het", "v1_hom", "v1_carriers", "revel", "gnomad_maf", "clnsig"]
    write_tsv(sorted(joe, key=lambda r: r["pos_hg38"]),
              args.out_dir / "joe_v1_variants.tsv", joe_cols)

    # ----- Output 3: The "115 carrier" variant search -----
    print(f"\n=== Searching for the '115 occurrence' variant (Doug's flag) ===")
    candidates = sorted([v for v in v2 if v["carriers"] >= 50],
                        key=lambda v: -v["carriers"])
    print(f"  Variants with ≥50 carriers (top 10):")
    print(f"  {'ID':<25} {'Carriers':>9} {'v2 MAF':>10} {'gnomAD':>10} {'Func':<15} {'ClinVar':<30}")
    for c in candidates[:10]:
        print(f"  {c['id']:<25} {c['carriers']:>9} {fmt(c['v2_maf']):>10} "
              f"{fmt(c['gnomad_maf']):>10} {c['func_class']:<15} {c['clnsig'][:28]:<30}")
    write_tsv(candidates,
              args.out_dir / "variant_115_candidates.tsv", cols)

    # ----- Output 4: pLOF + missense burden-eligible -----
    print(f"\n=== Burden-eligible variants (PMBB v2) ===")
    print(f"\n  ALL pLOF (regardless of MAF):")
    for v in sorted(plof, key=lambda r: r["pos"]):
        print(f"    {v['id']:<25} carriers={v['carriers']:>3} MAF={fmt(v['v2_maf']):>10} "
              f"passes_MAF<0.001={v['passes_maf_001']!s:<5}  {v['exonic_func']}")

    print(f"\n  Missense with REVEL>0.6 (regardless of MAF):")
    miss_high_revel = [v for v in missense if v["passes_revel_06"]]
    for v in sorted(miss_high_revel, key=lambda r: r["pos"]):
        print(f"    {v['id']:<25} carriers={v['carriers']:>3} MAF={fmt(v['v2_maf']):>10} "
              f"REVEL={fmt(v['revel'])} passes_MAF<0.001={v['passes_maf_001']!s:<5}")

    # ----- Output 5: Variants Joe had with v1 carriers > 0 -----
    print(f"\n=== Joe's PMBB v1 variants with carriers > 0 (top by v1_carriers) ===")
    joe_sorted = sorted(joe_with_carriers, key=lambda r: -r["v1_carriers"])
    print(f"  {'hg38_pos':<25} {'v1_carriers':>11} {'gnomAD':>10} {'Function':<15} {'ExonicFunc':<25}")
    for j in joe_sorted[:20]:
        print(f"  {j['constant_id'][:24]:<25} {j['v1_carriers']:>11} "
              f"{fmt(j['gnomad_maf']):>10} {j['func'][:14]:<15} {j['exonic_func'][:24]:<25}")

    # ----- Output 6: Cross-reference Joe v1 ↔ PMBB v2 -----
    # ZNF175 hg38 → hg19 has consistent shift of -503,253 in this region
    # (verified against Daniel's 10 manual lifts in runbook lines 213-232)
    HG38_TO_HG19_SHIFT = -503_253
    v2_by_pos = {}
    for v in v2:
        v2_by_pos.setdefault(v["pos"], []).append(v)

    print(f"\n=== Cross-referencing Joe v1 list with PMBB v2 (hg38→hg19 shift = {HG38_TO_HG19_SHIFT}) ===")
    cross_ref = []
    for j in joe:
        hg19_pos = j["pos_hg38"] + HG38_TO_HG19_SHIFT
        v2_matches = v2_by_pos.get(hg19_pos, [])
        # Try to find exact ref/alt match
        exact_match = None
        for v2v in v2_matches:
            if (v2v["ref"] == j["ref"] and v2v["alt"] == j["alt"]):
                exact_match = v2v
                break
        if exact_match is None and v2_matches:
            exact_match = v2_matches[0]

        is_plof_joe = classify_function(j["func"], j["exonic_func"]) == "pLOF"
        is_missense_joe = j["exonic_func"] == "nonsynonymous SNV"

        cross_ref.append({
            "joe_id": j["constant_id"],
            "pos_hg38": j["pos_hg38"],
            "pos_hg19": hg19_pos,
            "joe_func": j["func"],
            "joe_exonic_func": j["exonic_func"],
            "joe_is_plof": is_plof_joe,
            "joe_is_missense_high_revel": is_missense_joe and (j["revel"] or 0) > 0.6,
            "v1_carriers": j["v1_carriers"],
            "joe_revel": j["revel"],
            "joe_gnomad": j["gnomad_maf"],
            "v2_match_id": exact_match["id"] if exact_match else "—",
            "v2_carriers": exact_match["carriers"] if exact_match else None,
            "v2_maf": exact_match["v2_maf"] if exact_match else None,
            "v2_passes_maf_001": exact_match["passes_maf_001"] if exact_match else False,
            "v2_burden_eligible": exact_match["burden_eligible"] if exact_match else False,
            "in_v2_annotation": exact_match is not None,
        })

    # Filter to variants that mattered for burden test (pLOF or missense + REVEL>0.6) AND had v1 carriers
    burden_relevant = [c for c in cross_ref
                       if (c["joe_is_plof"] or c["joe_is_missense_high_revel"]) and c["v1_carriers"] > 0]

    print(f"\n  Joe's burden-relevant variants (pLOF or missense REVEL>0.6, v1 carriers > 0): {len(burden_relevant)}")
    print(f"  {'hg38_pos':<22} {'hg19_pos':>10} {'type':<8} {'v1_car':>6} {'v2_car':>6} {'v2_maf':>10} {'in_v2':<6} {'burden_OK':<9}")
    for c in sorted(burden_relevant, key=lambda x: -x["v1_carriers"]):
        joe_type = "pLOF" if c["joe_is_plof"] else "missense"
        print(f"  {c['joe_id'][:21]:<22} {c['pos_hg19']:>10} {joe_type:<8} "
              f"{c['v1_carriers']:>6} {str(c['v2_carriers'] or '—'):>6} "
              f"{fmt(c['v2_maf']):>10} "
              f"{'YES' if c['in_v2_annotation'] else 'NO':<6} "
              f"{'YES' if c['v2_burden_eligible'] else 'NO':<9}")

    # Critical finding: variants Joe used that are NOT burden-eligible in v2
    dropped_in_v2 = [c for c in burden_relevant if not c["v2_burden_eligible"]]
    print(f"\n  *** {len(dropped_in_v2)} burden-relevant variants from Joe's list are NOT burden-eligible in v2 ***")
    for c in sorted(dropped_in_v2, key=lambda x: -x["v1_carriers"]):
        reason = []
        if not c["in_v2_annotation"]:
            reason.append("NOT_IN_V2_ANNOT")
        elif not c["v2_passes_maf_001"]:
            reason.append(f"MAF_EXCEEDS_0.001 (v2_MAF={fmt(c['v2_maf'])})")
        else:
            reason.append("filter_unclear")
        print(f"    {c['joe_id'][:30]:<32} v1_carriers={c['v1_carriers']:>5}  "
              f"v2_carriers={c['v2_carriers']}  reason={','.join(reason)}")

    cols_xref = ["joe_id", "pos_hg38", "pos_hg19", "joe_func", "joe_exonic_func",
                 "joe_is_plof", "joe_is_missense_high_revel",
                 "v1_carriers", "joe_revel", "joe_gnomad",
                 "v2_match_id", "v2_carriers", "v2_maf",
                 "v2_passes_maf_001", "v2_burden_eligible", "in_v2_annotation"]
    write_tsv(sorted(cross_ref, key=lambda r: r["pos_hg38"]),
              args.out_dir / "joe_v1_v2_crossref.tsv", cols_xref)

    print(f"\n✓ Outputs written to {args.out_dir}")


if __name__ == "__main__":
    sys.exit(main())
