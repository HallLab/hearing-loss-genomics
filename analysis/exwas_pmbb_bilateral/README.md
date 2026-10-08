# ExWAS PMBB — bilateral sensorineural hearing loss

**Status:** Phase 1 complete · Phases 2–3 inherited · Phases 4–5 to run
**Release:** PMBB-Release-2026-4.0 · **Test:** SAIGE-GENE+ rare-variant gene burden

---

## What this is

The analysis this group stands behind: the corrected pipeline from
`analysis/elena_replication/`, run on the phenotype the clinical lead considers the phenotype.

It is **not** a replication and carries no second arm. The replication exists to answer *what was
done, and what was wrong with it*, and it stays frozen for that purpose. This folder answers a
different question: **what does the data say, done the way we would defend.**

## Why it is separate

Every correction here is already justified somewhere else, and mixing them would destroy both.
A replication that also changes the phenotype can no longer attribute a divergence to anything —
you cannot tell a fixed defect from a redefined question. So the two live apart, and the dependency
runs one way: this folder was seeded from the replication's outputs and never writes back to it.

## Why "bilateral sensorineural"

From the meeting of 2026-10-02, Douglas Epstein:

> "we typically exclude one-sided hearing loss, and we just focus on the bilateral sensorineural
> hearing loss"
>
> "unilateral hearing loss we exclude, because it's less likely genetic, more likely environment
> related"

The replication's case definition is `SO_396`, the phecodeX parent, which bundles conductive, mixed,
unilateral and unspecified hearing loss together. That is a faithful reproduction of what the
pipeline did. It is not the phenotype anyone intended to study.

The cost is large and has to be stated up front: **cases fall from 6,752 to 3,164.** See
[`PREMISES.md`](PREMISES.md) for what each premise moves, and
[`pipeline.md`](pipeline.md) for the stages.

## Self-contained, on purpose

Everything this analysis needs is **inside this folder**, including artifacts that were copied
rather than recomputed. Masks, the ancestry PCA and the GRMs are phenotype-independent, so they are
reused as-is from the replication — but they live here, byte-identical, not referenced across
directories.

The point is that this folder can be handed to someone, or moved, or archived with a manuscript,
without a trail of references to a sibling directory that may change. What was copied and what was
built is recorded in [`PROVENANCE.md`](PROVENANCE.md), and every copied artifact was verified
against its source.

**One exception, and it is stated rather than hidden.** SAIGE step 2 reads the per-chromosome exome
genotypes at
`analysis/elena/rarevariant_geneburden/plink_deduplicated/`, which is **438 GB**. Copying that is not
sensible and would not make anything safer — it is a mechanical conversion of the release, the same
class of input as `/static/PMBB/` itself, and it is the only path in this folder that points
outside. If it ever moves, one line in `phase_4/scripts/03_saige_step2.bsub` changes.

Everything that encodes a *decision* — masks, PCA, GRMs, the phenotype — is local.

## Layout

```
phase_1/   phenotype and cohort            BUILT HERE
phase_2/   variant masks                   copied, validated upstream, not re-run
phase_3/   ancestry PCA + covariates       PCA copied; covariates built here
phase_4/   SAIGE step 1 and step 2         BUILT HERE
phase_5/   results, figures, tables        BUILT HERE
andre_notes/  reference and plain-language notes
```

## What it inherits without re-deciding

These were established in the replication, each against evidence recorded there, and are not
re-litigated:

- the cohort is framed on **exome** samples, not the imputed `.fam`
- tinnitus exclusion reads the OMOP `observation` table, where PMBB v4 moved the standard codes
- masks are rebuilt: SpliceAI gate enforced, REVEL parsed as a list, protein-coding genes only
- PC counts are **5 / 4 / 3**, from the scree
- `Batch` is declared categorical
- the `ALL` mask is not run — which the clinical lead independently objected to on 2026-10-02
