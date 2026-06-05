# NB 05 — VEP/dbNSFP vs Biofilter 4 reconciliation

## Coverage
- VEP annotated: 384/386 | BF4 found: 97/386
- VEP-only (BF4 missed): 287 variants
  -> of which pLOF: 19,
     with AlphaMissense: 131

## pLOF
- VEP (IMPACT=HIGH): 24 | BF4 (LOFTEE): 6
- 5 of 6 BF4-LoF are also VEP-HIGH. The 1 difference (51573407:CAG:C) is a TRANSCRIPT-CHOICE
  artifact: frameshift on one transcript (BF4/LOFTEE) vs splice_region/intron on VEP --pick's transcript.
- extra pLOF only VEP found: 19 (mostly BF4 not_found, ultra-rare private burden candidates)
- NOTE: VEP --pick can UNDER-call pLOF (one transcript only); --most_severe would catch 51573407.
  For a pLOF-focused pass, re-run VEP with --most_severe.

## AlphaMissense
- BF4: 0 | VEP/dbNSFP: 177

## Verdict
VEP/dbNSFP is strictly superior for THIS task: same calls where both annotate, but ~100% coverage
(by coordinate) vs BF4's ~25%, 4x more pLOF (the private burden candidates), and AlphaMissense where
BF4 had none. IMPORTANT: BF4's low coverage here is NOT inherent — the loaded BF4 DB was filtered to
gnomAD AC >= 5, so AC < 5 variants were excluded (they exist in gnomAD; BF4-found min ac = 5). Reloading
BF4 without that cutoff would recover most not_found. BF4's standing edge = convenience + extra
precomputed scores (SpliceAI/Pangolin/CADD) for variants in its DB.

## Caveats
- VEP --pick chooses ONE transcript -> consequence disagreements with dbNSFP (e.g. intron on the picked
  transcript but missense -> AlphaMissense on another). Also under-calls pLOF (see 51573407). --most_severe fixes this.
- LOFTEE HC/LC not run (human_ancestor/gerp data absent) -> VEP pLOF here = IMPACT=HIGH, not LOFTEE-filtered.

## Output
- vep_vs_bf4_reconciliation.csv
