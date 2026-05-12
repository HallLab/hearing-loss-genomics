# Kickoff Meeting Summary — Hearing Loss Genomics

**Date:** Early 2026 (project kickoff)
**Participants:** Molly Hall (PI), Douglas Epstein (collaborator), Nikki Palmiero, Andre Rico
**Reference paper:** Hui et al., PLOS Genetics 2023 (doi: 10.1371/journal.pgen.1010584)

---

## 1. Scientific background (Doug Epstein)

### Hearing loss landscape
- ~50% of Americans aged 70+ have some hearing loss (HL); ~50% disabling.
- ~50% of children failing newborn screening have a genetic basis for HL.
- ~100 known non-syndromic congenital HL genes (more if syndromic forms counted).
- Age-related HL heritability: 30–70%; 89 independent GWAS loci known, all with small effects.

### Adult-onset HL is understudied
Adult patients are rarely assessed genetically because results don't currently change clinical management. But:
- Translational therapy may benefit from linking adult HL to well-characterized congenital HL genes.
- Biobank cohorts allow large-scale study without targeted recruitment.

### The ZNF175 story (priority of this project)
- **Origin:** Hui et al. ran a gene-burden test of putative loss-of-function (pLoF) variants across PMBB phenotypes. Standard positive controls validated the approach (CFTR/CF, BRCA/cancer, TTN/cardiomyopathy).
- **Novel hit:** ZNF175 — a KRAB-zinc finger transcription factor (KRAB-ZFP) on chr19. KRAB-ZFPs silence transposable elements; they evolve rapidly because transposons do.
- **Mouse work (Epstein lab):** The mouse syntenic region has a family of 4–5 paralogs. Only one — *Zfp719* — produces an HL phenotype when knocked out. But there's a catch:
  - pLoF on a **true wild-type background → no phenotype**.
  - The common lab mouse strain carries a known hypomorphic mutation in a separate deafness gene; the ZNF175 ortholog phenotype only emerges in that sensitized background.
  - When corrected, *Zfp719* loss causes a progressive HL phenotype with a worked-out mechanism.
  - **The phenotype has now been rescued** (newer than the last conversation with Molly).
- **Targets:** Human ZNF175 binds older transposon classes; the mouse ortholog binds mouse-specific (newer) transposons. Species-specific transposon regulation.

### The PMBB puzzle to solve
- The original GWAS/burden signal came from **8 cases**.
- A curated case-control cohort (~140 individuals with ZNF175 mutations + matched controls) showed roughly equal HL rates (~18%) between cases and controls.
- The genetic *cause* of HL in those carriers was never elaborated.
- **Hypothesis:** ZNF175 acts as a modifier; HL phenotype requires a second hit in a known HL gene. The 8 driving cases likely carry such second hits.

---

## 2. PMBB resources available

- **40,627 individuals** total in Hui et al.'s analysis (75% European, 25% African American ancestry).
- **16,000 individuals with audiograms** in PMBB total, but only **~4,000 with both audiograms AND exomes** — this 4K subset is the high-value cohort for quantitative-trait analysis.
- Audiogram data lives **outside the EHR** but the ENT team has matched audiograms to PMBB IDs.
- Building a "database of hearing loss variants" is an ongoing parallel effort that should incorporate these 4K cases.

---

## 3. Immediate priorities

### Priority 1 — Deep-dive the 8 signal-driving exomes
- Scan their exomes for variants in *other* HL genes that could explain their phenotype.
- This is a **screen / hypothesis-generating step**, not a formal statistical test on its own.
- If hits emerge, use them to design a larger-scale analysis (potentially in other biobanks).

### Priority 2 — Statistical follow-up on the curated 140-case cohort
- Use Doug's curated cohort (~140 ZNF175-mutation carriers + matched controls).
- Split into: (a) ZNF175-mutation + HL, (b) ZNF175-mutation, no HL.
- Look for differential burden of variants in other HL genes between (a) and (b).
- Caveat: small N — power will be limited.
- Could be run as: exome-wide rare variant burden / GWAS scan on this subcohort, after removing the ZNF175 signal itself.

### Priority 3 — Document & port the original pipeline
- Daniel Hui ran the original pipeline. Already shared the scripts.
- Joe Park provided methodological recommendations.
- Newer rare-variant approaches (developed within the last ~6 months) should also be applied as a methodological comparator.
- Nikki/Andre should review, document, and re-run with current methods.

---

## 4. Future / later phases (mentioned but deferred)

- **UK Biobank replication of ZNF175:** signal was maintained in UKBB, but UKBB HL phenotype is self-reported (less precise). Still potentially worth doing because participant numbers will exceed 8.
- **All of Us extension:** ZNF175 not yet tested in AoU. Nikki has AoU access/experience.
- **Quantitative-trait analysis** on the ~4K audiogram+exome subset of PMBB.
- **Targeted PRS / gene-based risk scores** as predictive tools — Nikki has expertise here.
- **Gene × environment:** Molly mentioned a separate UKBB analysis where early-life trauma showed strong association with HL phecodes. Calcium dysregulation also surfaced. Doug noted the inner ear / brain split. Marked as exciting but not now.

---

## 5. Andre's role (per the meeting)

- Primary expertise: programming, software infrastructure, database design, PRS.
- Will pair with Nikki — they decide between themselves who runs what.
- Junior on the biological side but supporting all technical aspects.

## 6. Ops cadence

- Hall Lab has **weekly meetings** with Nikki/Andre; Doug joins when needed for analysis-plan reviews.
- Next concrete deliverable: an **analysis plan** to review with Doug before executing.
