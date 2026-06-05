# NB 03 — Biofilter 4 annotation of the ZNF175 union table

Report: annotation_master_variant (most_severe_only). Input: 384 variants (chr:pos:ref:alt; '*' excluded).
Output: cmp_union_annotated_bf4.csv (union table + BF4 columns).

## What BF4 delivered (strengths)
- found: 97 / 386 variants  (rich precomputed scores: REVEL, CADD, SpliceAI, gnomAD AF, VEP consequence)
- pLOF (LOFTEE): 6 variants flagged (1 HC + 5 LC) -- frameshift / stop_gained

## What BF4 did NOT answer (weaknesses — to revisit with VEP)
- not_found: 287 variants (75% of submitted). ROOT CAUSE: the loaded
  BF4 DB was filtered to **gnomAD AC >= 5** (a build-time choice) -> variants with AC < 5 are excluded. They DO
  exist in gnomAD (10 not_found spot-checked on the gnomAD site -> all AC < 5); BF4-found variants have min ac = 5.
  So this is a CONFIGURABLE cutoff, not an inherent gnomAD limit -> reversible by reloading BF4 without the AC filter.
  (low cohort-MAF tracks low gnomAD AC, hence found-rate ~11% in the <0.01% bin.)
- AlphaMissense: NOT populated for ZNF175 in this DB build -> 0 scores despite found missense variants.
- '*' spanning-deletion alleles (2) not submitted (not real alleles).

## Implication
BF4 is strong for variants ABOVE its AC cutoff (pLOF + precomputed scores). To include the rare/private burden
candidates, either reload BF4 without the AC filter, or annotate by coordinate. NEXT: run VEP + dbNSFP4.5a on the
region VCF (annotates ALL variants) and compare against BF4 to quantify exactly what BF4's AC cutoff excluded.

## Outputs
- cmp_union_annotated_bf4.csv, bf4_coverage_by_rarity.csv, znf175_bf4_annot.csv
