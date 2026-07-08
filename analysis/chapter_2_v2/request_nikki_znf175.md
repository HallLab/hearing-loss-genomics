# Draft message to Nikki — ZNF175/tinnitus: reconciling 4 vs 8 carrier-cases

**Subject:** ZNF175 → tinnitus replication (Freeze One): can we find the code behind the "8 carrier-cases"?

Hi Nikki,

Quick update and one main request on the ZNF175 / tinnitus replication (Park 2021, Freeze One WES).

**What we did.** We applied the Park pipeline to Freeze One and get **34 qualifying rare-pLOF carriers** (matching Park's ~35) and **4 tinnitus carrier-cases**. We stress-tested that 4 and it's robust to every choice we can vary — variant set, ICD code list, rule-of-2, event-vs-date counting, the canonical phecode definition (PheWAS `createPhenotypes` + control exclusions), and the annotation tool (VEP vs ANNOVAR/REVEL/AlphaMissense). We also confirmed the **same 4 individuals** carry the signal in both v1 (11K) and v2 (44K), by variant fingerprint.

**Where we think the 8 comes from.** Our best reconstruction is **8 = 6 + 2**:
- **6** = a looser phenotype — broad ear/hearing phecodes at **rule-of-1** (≥1 date) instead of strict rule-of-2 tinnitus — on the linked carriers. We reproduce this exactly.
- **+2** = ~2 of **7 carriers we cannot phenotype**: their WES GENO_IDs aren't in the Freeze One `Demographics` crosswalk (and we can't find any UPENN↔PMBB_ID bridge on disk). At the linked-carrier rate (~22%), we'd expect ~1.5–2 of those 7 to be cases — which lands the total at 8.

**Main ask.** Could you point us to the **original code/scripts that produced the "8 carrier-cases"** (Daniel's ZNF175 pipeline, or wherever that number came from)? Seeing the exact qualifying-variant filter, phenotype definition, and how the linkage was handled would let us reconcile in one step — and confirm (or correct) the 6 + 2 reconstruction.

**Secondary ask (if the code isn't handy).** The GENO_ID ↔ PT_ID crosswalk for these **7 GENO_IDs** (attached, *List A*), so we can check their tinnitus status ourselves. Prediction: ~2 are ear/tinnitus cases. We also have *List B* — 77 tinnitus cases with no genotype — if the reverse mapping is easier on your end.

Happy to share the full write-up (`preconclusion_znf175_4_vs_8.md`) and the two lists. Thanks!

Andre

---
*Attachments: `nb03_list_A_wes_carriers_no_phenotype.csv` (7 GENO_IDs + variants), `nb03_list_B_tinnitus_no_genotype.csv` (77 PT_IDs).*
