# NB 04 — VEP + dbNSFP4.5a annotation of the ZNF175 union (no Biofilter)

Pipeline: VEP (cache 109, --pick) for consequence+IMPACT; ANNOVAR dbNSFP4.5a(canonical) for AlphaMissense.
Annotates by COORDINATE -> covers ALL variants (incl. private). Output: cmp_union_annotated_vep.csv.

## Coverage
- VEP consequence: 384/386 variants annotated (vs BF4 found 97/386)
- pLOF (IMPACT=HIGH): 24 (vs BF4 6)
- AlphaMissense scores: 177 (vs BF4 0 for ZNF175)

## Why VEP beats BF4 here
- BF4 = gnomAD-known-variant DB -> 75% not_found (private/ultra-rare escape; ~11% found <0.01% MAF).
- VEP/dbNSFP compute consequence/AlphaMissense from the GENOME by coordinate -> ~100% coverage at every rarity.

## Next
- NB 05 (planned): full VEP-vs-BF4 reconciliation -- agreement on shared variants, what each adds,
  and whether LOFTEE HC/LC (needs human_ancestor/gerp data) changes the pLOF set.

## Outputs
- cmp_union_annotated_vep.csv, znf175_vep.tab, znf175_annovar.hg38_multianno.txt, znf175_union.vcf
