# Meeting Summary — Research Vision & Elena Onboarding

**Date:** 2026-06-03
**Participants:** Molly Hall, Douglas Epstein, Nikki Palmiero, Andre Rico, Elena (incoming summer undergrad)
**Format:** Two halves — (1) Doug's full research vision for the new lab member, (2) Molly framing how our biobank approach fits in, plus Elena's onboarding logistics.
**Reference papers:** Park et al., *Nature Medicine* 2021 (original ExWAS, Dan Rader & Marilyn Ritchie senior authors); Hui et al., *PLOS Genetics* 2023 (Daniel's replication).

---

## 1. Scientific background (Doug Epstein)

### Origin of the signal
- The project traces back to a PMBB putative-loss-of-function (pLoF) screen — ~1,000 genes carrying pLoF in ≥30 individuals, associated against broad phenotypes (Park et al., *Nature Medicine* 2021). Standard positive controls validated the approach (BRCA/breast cancer, CFTR/cystic fibrosis).
- The novel hit: **ZNF175** (a KRAB zinc-finger transcription factor) associated with **tinnitus**. Some carriers had audiometrically-defined hearing loss (HL), some self-reported HL, some neither — i.e., **incomplete penetrance**, and the original association was driven by a small number of cases.

### KRAB-ZFP biology
- KRAB-ZFPs are the largest TF family (~400 members), repressors whose primary job is silencing **transposable elements (TEs)** — mostly in the germline, but also in somatic tissue.
- These genes evolve rapidly because the transposons they police are species-specific.

### Mouse model
- Mouse has **no 1:1 ortholog** of ZNF175 — a 4-to-1 expansion on the syntenic region. Only **one paralog, *Zfp719***, is meaningfully expressed in the ear and the only one with a phenotype → treated as the functional ortholog.
- Mice are deaf; progressive. **Homozygotes:** profound HL at all frequencies by ~3 months. **Heterozygotes:** progressive but high-frequency only.
- **Mechanism (mouse):** a transposon-derived (LTR) promoter was co-opted into a **lncRNA** expressed *antisense* to a **calcium efflux pump** critical for clearing Ca²⁺/K⁺ from cochlear hair cells. Normally the lncRNA falls and the efflux pump rises after birth; in mutants the lncRNA stays high → pump stuck at ~50% → calcium dyshomeostasis → loss of hair-cell synapses over time. **Deleting the lncRNA fully restores hearing** in both het and homozygous mice.

### Three interaction layers (all observed in mouse)
- **Gene–gene:** mechanism was only uncovered on a calcium-sensitized background. A second mutation is what drove the mechanistic insight.
- **Gene–environment:** noise exposure. Corrected mutants are ~fine until exposed to loud noise and then fail to recover (parallels a 21-y-o clinic patient who lost hearing after a concert).
- **Gene–immune:** mutation activates innate immunity; tissue-resident macrophages invade the organ of Corti and are **protective** (removing them makes the phenotype worse, faster). A specific common variant in **CX3CR1** (chemokine receptor, ~20% of population) abolishes that protective effect — crossing it onto the mutant background lost protection.

### ⚠ Critical translation caveat
- Mouse and human KRAB-ZFPs bind **different (species-specific) TEs**. The mouse lncRNA mechanism is **mouse-specific and does not exist in humans**. So the human *downstream* mechanism is very likely different.
- What is **shared** is the *upstream* regulatory network keeping the gene expressed in the ear / hair cells (similar promoter architecture). Selective pressure is on keeping *some* KRAB-ZFP expressed in the ear to silence TEs; which gene and which downstream targets differ by species.
- **Implication for us:** the mouse does not hand us a concrete human molecular target to test directly. It hands us the **"modifier + second hit" framework**. The one concrete, named, testable human handoff is **CX3CR1**.

---

## 2. The statistical puzzle & our approach (Molly Hall)

### The 11K → 45K signal loss (the crux)
- ZNF175–tinnitus was significant at ~11K exomes (Joe Park, the original ExWAS) but the signal was **lost at ~45K exomes** (Daniel's later analysis — negative/non-significant).
- The statistical-genetics view (Ian Matheson): signal lost = not important in humans.
- **Molly's counter:** the two analyses were run by **different people with possibly different pipelines**. We must rule out a *technical* explanation before invoking biology. (Noted detail: original analysis was pLoF; replications in BioMe / DiscovEHR / UK Biobank shifted to *missense* and would not have passed initial discovery thresholds.)

### Where we are
- Andre has **reproduced Daniel's (45K) results** — confirmed negative/non-significant.
- **Next step:** track down and reproduce **Joe Park's original 11K pipeline** to get a true apples-to-apples comparison. Joe is hard to reach (left Penn, fellowship in NY). Look for a GitHub in the *Nature Medicine* 2021 paper, or the Ritchie Lab files on LPC. Note: Joe did reanalyze on 45K and reproduced Daniel's (negative) result.

### New direction — PMBB v4
- PMBB **version 4** is officially out (more patients than the v3 work). Plan a fresh **exome-wide rare-variant study of hearing loss** (all genes, not just ZNF175), plus **tinnitus**; **Meniere's** tabled (likely too few cases in PMBB for rare-variant; GWAS already done via meta-analysis of other biobanks).
- **GWAS vs XWAS:** leaning XWAS (exome focus matches the gene-centric interest and pairs naturally with the rare-variant work on the same genes). GWAS of HL already published and "not that interesting." Methodology/QC is shared, so the choice can be settled later.
- **Rare-variant method:** Molly is moving away from the word "burden" — likely **SKAT**-type methods, since the toolkit has evolved since Daniel's work. Nikki to lead method choice (she's been doing cutting-edge rare-variant work for the Cardiovascular Institute).

### Long-game framework — "variant status" (from the HCM/sarcomeric work)
- Modeled on Nikki/Charlene Day's hypertrophic-cardiomyopathy approach: take an **adjudicated ClinVar list** of HL genes/variants → assign each individual a **carrier (variant) status** across the curated set → ask what fraction carry an **additional hit** in a modifier gene (e.g., a macrophage-associated gene) → and whether that second hit **changes phenotype** vs carriers without it.
- This is the tractable, **descriptive** route to the gene–gene interaction question — it sidesteps the (admittedly under-powered) formal interaction test. Flagged as several steps down the line but the most promising long-term direction.

### Age-cohort idea — considered and dropped
- Mouse shows an age component (single-mutant mice start losing hearing between 12–15 months); idea was to stratify the biobank by age cohort, and to check whether an **age shift between the 11K and 45K sets** could explain the signal loss.
- Dropped for onset analysis: ICD codes give **diagnosis date ≠ age of onset** — "too dirty." (Checking the age *distribution* difference between 11K and 45K is still a reasonable thing to look at as a possible driver.)

---

## 3. Decisions & action items

| # | Action | Owner |
|---|--------|-------|
| 1 | Reproduce **Joe Park's original 11K pipeline** (find GitHub in *Nat Med* 2021, or Ritchie Lab files on LPC); compare apples-to-apples with Daniel's 45K | Andre (+ Nikki) |
| 2 | Add the **Park et al. *Nature Medicine* 2021** paper to project documentation for Elena | Andre |
| 3 | Check the **age distribution** difference between the 11K and 45K sets (possible confounder of signal loss) | Andre |
| 4 | **Onboarding tutorial:** Elena reproduces Andre's reproduction of Daniel's burden test (173-gene set is the easy entry point) — purely to get fluent with PMBB / LPC / scripts | Elena (Andre supports) |
| 5 | After tutorial: build a **new analysis plan** with state-of-the-art methods — start with **GWAS/XWAS** (simpler QC), then **rare-variant (SKAT)** of HL on PMBB v4 | Elena (Nikki mentors) |
| 6 | Provide Elena with example **GWAS and rare-variant analysis plans**; send foundational GWAS + rare-variant papers (incl. Nikki's HCM rare-variant references) | Nikki / Andre |
| 7 | **Long-term:** match audiograms to PMBB IDs (currently unmatched) — gates any quantitative-trait / age-onset work | Andre (with ENT coordination) |

## 4. Open questions / deferred

- GWAS vs XWAS final call (leaning XWAS) — table until analysis-plan stage.
- Whether/how to fold **CX3CR1** into the rare-variant analysis (the cleanest gene–immune second-hit test). Worth prioritizing as the one concrete mouse→human handoff.
- Meniere's rare-variant — likely under-powered in PMBB; coordinate with Bogdan/Ian (drivers of the Meniere's GWAS) before stepping in.
- UK Biobank is **off the table** (data breach / access frozen since early May; plan expected ~June).

## 5. Logistics

- **Elena:** LPC access being approved (pending Molly's OK on the systems ticket); no data access until training (GCP/HIPAA, ~3-hr) is complete — can explore Project Hall on LPC meanwhile. In the Hall lab Tue/Wed; will shadow wet-lab (Sophie starts ~24th) and possibly the HL clinic (Tiff). May use Andre's desk while he's home.
- **Andre:** working mostly from home in June (wife recovering from a back injury; Brazil trip canceled; green card received ~2 months ago).
- No Hall lab meetings for most of June; Molly remote in July.

---

## My read (Andre) — analytical notes

- The whole near-term plan hinges on **technical artifact vs real null** for the 11K→45K loss. Ruling out the pipeline difference is correct and necessary — but a small-case signal that vanishes at 4× the sample is the textbook *winner's-curse* pattern, so even a perfect pipeline match may simply not replicate. That's a valid answer, not a failure.
- The mouse story is compelling but **does not give us a human target** except CX3CR1 (different TEs, mouse-specific lncRNA). Treat the mouse as motivation for the *second-hit framework*, and CX3CR1 as the one concrete testable variant.
- **Power is the elephant** — Doug acknowledges the interaction analysis won't be statistically powered. The honest, tractable path is the **descriptive "variant-status" framework**, not formal interaction testing.
- **Scope risk:** this was a brainstorm, not a plan. The three real near-term deliverables are (a) reproduce Joe's 11K, (b) Elena's tutorial, (c) v4 HL ExWAS. Everything else is later-phase.
