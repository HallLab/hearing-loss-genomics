# NB 02 — ZNF175 carrier-cases via canonical phecode (v1, 11K)

Carriers: 34 (matches Park's ~35). In phecode table (linked): 27 | unlinked: 7.

## Carrier-cases by definition
- our earlier RAW ICD (388.3x/H93.1x, rule-of-2): 4
- **canonical PHECODE 389.4 (tinnitus, + control exclusions): 4**
- phecode 389 (hearing loss, broader): 4
- phecode 389.1 (sensorineural HL): 2

## Sensitivity grid (does any reasonable definition reach 8?)
- tinnitus 389.4 is **4** under BOTH distinct-date≥2 and event-count≥2 — event-vs-date counting does NOT move it.
- it reaches 5 only with rule-of-1 (presence); broadening to hearing-loss 389 gives 5 (event≥2) or 6 (rule-of-1);
  tinnitus-OR-HL rule-of-1 = 6.
- **No definition on the linked carriers reaches 8** (max = 6).

## Annotation robustness (§2c) — VEP vs ANNOVAR/REVEL/AlphaMissense
- ceiling = 11 tinnitus cases carry a rare ZNF175 variant; only 4 are qualifying pLOF.
- the other 7 carry synonymous / intron / UTR / **benign** missense (REVEL≪0.5, AlphaMissense B) —
  no annotator (ANNOVAR LOF, REVEL≥0.5, AlphaMissense) upgrades them. Carrier-cases are **robust to annotation tool**, not just counting rule.

## v1 ↔ v2 same individuals (§2d)
- the 4 carrier-cases are the SAME 4 people in v1 and v2 (variant fingerprint {'19:51588428:CA:C': 2, '19:51587727:CAAAG:C': 1, '19:51588214:CAG:C': 1} identical across freezes; IDs only recoded).
- "no new cases in v2" is NOT strange: carriers stay ~3x enriched (real, modest); observed 4 matches OR~3.5,
  sits above null (~1.4) and far below Park's OR~14.6 (which predicts ~16). Winner's curse made concrete —
  the discovery effect was anchored on these 4, and quadrupling the cohort added ~0 cases.

## Reading
- Canonical phecode tinnitus gives 4 carrier-cases (vs our raw-ICD 4) — robust; event-vs-date and the rollup do
  not change it.
- The gap to the PI's 8 is NOT closed by the phenotype definition. The remaining lever is the 7 carriers
  with no phenotype linkage (their tinnitus status is unknown) — need the full PMBB ID map to assess them.

## Output
- znf175_carriers_phecode.csv (per-carrier phecode status)
- carrier_sensitivity_grid.csv (definition × carrier-case count)
