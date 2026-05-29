# Phase 7 — The Carrier Numbers, Explained Simply

**Audience:** anyone who wants to understand what's in our analysis without opening the spreadsheets.
**Companion to:** [`phase7_replication_report.md`](phase7_replication_report.md) (the technical version).
**Date:** 2026-05-15

---

## What is a "carrier"?

Every person has 2 copies of the ZNF175 gene (one from each parent). In some people, **one of those copies has a mutation that breaks the gene** — a "loss-of-function" (pLoF). We call those people **carriers**.

ZNF175 can have **many different mutations**, each at a specific position in the DNA. Each mutation has an "address":

```
chr19 : 51587727 : C : T
  │       │        │   │
  │       │        │   └── new letter (alt allele)
  │       │        └────── original letter (ref allele)
  │       └─────────────── position on the chromosome
  └─────────────────────── chromosome (19)
```

**Different people can carry different mutations in the same gene.** This is where the number confusion starts — counts depend on which mutation I'm looking at.

---

## The numbers that appear across our docs

Different reports cite different numbers (8, 90, 142, 3, 6, 1). **They're all correct — they just refer to different subsets.**

```
                    PMBB (43,731 sequenced individuals)
                              │
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
   RARE variant          COMMON variant       Other variants
   chr19:51587727        chr19:51581437       (Joe's curated list)
        │                     │                     │
     8 carriers          90 carriers          various
        │                     │                     │
     1 HL case           2 HL cases           ?
     (PMBB2106…)         (PMBB7501…,
                          PMBB1245…)
        │                     │                     │
        └─────────────────────┴─────────────────────┘
                              │
                              ▼
              Daniel aggregated everything into a single list:
                       
                       142 carriers
                  (= the "140-case cohort"
                     mentioned in kickoff,
                     approximately)
                              │
                              ▼
              Of the 142, with audiogram phenotype:
                              
                  3 are HL cases (SNHL=1)
                  129 are HL controls (SNHL=0)
                  10 unknown status (no audiogram)
```

### Each number, in one line

| Number | What it represents |
|---:|---|
| **43,731** | Individuals sequenced in PMBB total (the paper's cohort) |
| **8** | People carrying the **rare** variant `chr19:51587727` (= the "8 carriers" of the runbook) |
| **1** | Of those 8, how many have audiogram-confirmed HL → **PMBB2106731298975** |
| **90** | People carrying the **more common** variant `chr19:51581437` |
| **2** | Of those 90, how many have audiogram-confirmed HL → PMBB7501686571326, PMBB1245988577461 |
| **142** | People carrying **any** of the 10 pLoF variants Joe Park curated (= Daniel's canonical list) |
| **3** | Of those 142, how many have audiogram-confirmed HL (1 + 2 from the two variants above) |
| **6** | Specific people Joe Park **named explicitly** in an email to Daniel for follow-up |
| **6** | Of those 6, how many are in our list of 142 carriers (**ALL** of them) |
| **1** | Of those 6, how many have audiogram-confirmed HL → PMBB2106731298975 |

---

## The person who appears in EVERY subset: PMBB2106731298975 ⭐

This individual is the protagonist of the ZNF175 story:

| Subset | PMBB2106731298975 present? |
|---|---|
| 8 carriers of chr19:51587727 | ✓ yes (the only HL case among the 8) |
| 142 carriers (Joe's curated full list) | ✓ yes |
| 3 HL-cases among the 142 | ✓ yes |
| 6 named individuals from Joe Park | ✓ yes |
| The only Joe-named with confirmed HL | ✓ yes |
| The "1 case" noted in Daniel's runbook ("only case is PMBB2106731298975") | ✓ yes |

**In other words: all ZNF175 + HL findings in PMBB converge on this single individual.** The other 2 HL-cases (PMBB7501686571326, PMBB1245988577461) are additional findings that emerge when we aggregate more variants — including chr19:51581437, which Daniel noted may have been borderline-too-common and "should have been removed".

---

## So where does the "8 cases" of the kickoff come from?

The kickoff meeting summary states:

> *"the original GWAS/burden signal came from 8 cases"*

But what we actually see is **8 carriers**, of which **1 is a case**. Two possible explanations:

**Interpretation A (most likely):** In my notes from the kickoff meeting, Doug’s original phrasing may have been something like *“the signal came from 8 driving carriers, one of which was a case”* — but in the notes it ended up as *“8 cases”*.

**Interpretation B:** There's an additional filter producing 8 actual cases that we haven't yet recovered from Daniel's preserved files.

---

## Joe Park's 6 named individuals — clearing up the confusion

Joe Park (Daniel's collaborator on the paper) emailed Daniel **6 specific PMBB IDs** of ZNF175 carriers for detailed follow-up. These 6 are preserved in Daniel's runbook (lines 348-353).

### The 6 individuals and where they stand in our analysis

| Joe's named individual | In our 142 carriers? | HL status (audiogram) | Second-hits found |
|---|---|---|---:|
| PMBB8949342677388 | ✓ yes | NA (not audiogrammed) | 1 |
| PMBB4664277909557 | ✓ yes | NA | 1 |
| **PMBB2106731298975** ⭐ | ✓ yes | **HL case (SNHL=1)** | **2** |
| PMBB7083231520332 | ✓ yes | NA | 4 |
| PMBB9701760809542 | ✓ yes | NA | 4 |
| PMBB9508968070076 | ✓ yes | NA | 0 |

**Key points:**

1. **All 6 are in our 142-carrier list.** Nobody is missing.

2. **Only 1 of the 6 has audiogram-confirmed HL** (PMBB2106731298975). The other 5 show up as `NA` in our table because **they were not audiogrammed** — Joe likely identified them through other sources (clinical referral, EMR, otolaryngology records) that aren't part of our `cases_control.txt` phenotype file.

3. **NA does NOT mean control.** When someone is NA, it does **NOT** mean they don't have HL. It just means **we don't have audiogram data** for that person. To confirm or rule out HL in those 5, we'd need a chart review by the ENT team.

4. **Daniel preserved a re-contactable subset.** A file `with_email_living/ZNF175_pLOF_Joes.txt.gz` exists with 116 contactable carriers — this was likely the basis for a re-contact effort that wasn't completed.

---

## What we found in Phase 7

The core analysis: **for each of the 142 ZNF175 carriers, we looked at whether they ALSO carry deleterious variants in any of the 176 other known HL genes** (the gene set used in Phase 1).

### The 3 HL-case carriers — all have second-hits

| Person | Second-hit HL genes | Associated disease |
|---|---|---|
| **PMBB2106731298975** ⭐ | USH2A, MYO3A | Usher syndrome + DFNB30 |
| PMBB7501686571326 | COL11A1, COL11A2, GPSM2 | Stickler syndrome + DFNB82 |
| PMBB1245988577461 | **GJB2**, MYO7A, COL11A1, COL11A2 | **DFNB1 (the #1 HL gene worldwide)** + Usher + Stickler |

**All 3 carry hits in classic Mendelian HL genes.** 100% of cases. This is exactly the pattern Doug predicted from his mouse work: ZNF175 alone doesn't cause HL, but ZNF175 + a second hit in another HL gene = HL phenotype.

### And the 129 controls?

| Metric | Cases (N=3) | Controls (N=129) |
|---|---:|---:|
| Mean second-hits per person | 3.0 | 2.5 |
| % with at least 1 second-hit | 100% (3/3) | 93% (120/129) |

**The controls also have second-hits — almost all of them (93%).** This is informative:
- ZNF175 + second-hit alone is **not sufficient** to cause HL (otherwise the 93% of controls with second-hits should all have HL)
- But it might be **necessary or contributing** — all 3 cases have at least one

Hypotheses worth exploring:
- The "second-hit" definition is too permissive (we're using all 9,667 deleterious variants; restricting to only ClinVar pathogenic/likely-pathogenic might separate cases from controls more cleanly)
- Other factors enter the model (environmental, age of noise exposure, etc.)
- The audiogram-based phenotyping may miss subtle HL cases (some of the 93% of controls with second-hits may have mild HL that wasn't captured)

### The main limitation: N=3 is small

With only 3 cases, any formal statistical test ends up non-significant (our Mann-Whitney gave p=0.24). The qualitative story is strong, but for confirmation we'd need more cases. Paths forward:
- **UKBB** (UK Biobank — Daniel ran a partial UKBB analysis preserved in walkthrough Phase 15)
- **All of Us** (Nikki has access, mentioned in the kickoff)
- **PMBB v3 or v4** (newer releases — more sequenced participants)

---

## What we'll ask Daniel (3 questions in the email)

1. **The direction of the ZNF175 story** — did Doug bring this to you (based on his mouse Zfp719 work), or did ZNF175 first emerge from one of your burden tests?
2. **"8 signal-driving cases"** — are these the 8 carriers of chr19:51587727 (with 1 actually being an HL case), or is it a different subset?
3. **The "140-case cohort"** — is that the 142 carriers from Joe's list? Or is it something different and where is it documented?

---

## 30-second TL;DR

1. **ZNF175 has multiple mutations.** Different mutations have different carrier counts.
2. **8 carriers** = people with one specific rare variant (1 of them has HL).
3. **142 carriers** = aggregated list across Joe's 10 curated variants (3 of them have HL).
4. **3 cases** = the HL-cases among the 142, and **all 3 carry deleterious variants in classic Mendelian HL genes** (GJB2, USH2A, MYO7A, COL11A1, COL11A2, etc.).
5. **PMBB2106731298975** = the individual who appears in every subset; our "perfect example" of the second-hit hypothesis.
6. **Joe's 6** = all 6 are among the 142 carriers, but only 1 (PMBB2106731298975) has audiogram-confirmed HL — the other 5 were never audiogrammed (status unknown, not "no HL").
7. **Biological conclusion:** consistent with Doug's hypothesis, but N=3 is too small for formal statistical confirmation. Next steps: sensitivity analysis (ClinVar P/LP only), matched-controls comparison, external replication in UKBB / AoU.
