#!/usr/bin/env Rscript
# Build phecode phenotypes from ICD codes using the PheWAS package (Park's canonical method):
#   ICD-9/ICD-10 -> phecode mapping (Phecode Map 1.2) + rule-of-2 + control exclusions + sex restriction.
# Usage: Rscript run_phewas_createphenotypes.R <input.csv> <sex.csv> <pop_ids.txt> <out.csv>
#   input.csv : columns id, vocabulary_id (ICD9CM|ICD10CM), code, index
#               index = the RAW diagnosis DATE (one row per id×code×date), kept as a STRING.
#               ⚠️ Do NOT pre-aggregate to an integer count. createPhenotypes' rule-of-2 counts
#               distinct dates itself: default_code_agg does length(unique(index)) for CHARACTER
#               index, but sum(index) for NUMERIC — and mapCodesToPhecodes' distinct() collapses
#               identical (id,phecode,count) rows across ICD codes, which UNDERCOUNTS cases when a
#               pre-aggregated integer count is fed. Feeding raw dates (strings) is the correct form.
#   sex.csv   : columns id, sex (M|F)
#   pop_ids.txt: one id per line — the full population (defines controls)
#   out.csv   : phecode matrix (id + one column per phecode; TRUE=case, FALSE=control, NA=excluded)
suppressMessages(library(PheWAS))
a <- commandArgs(trailingOnly = TRUE)

d <- read.csv(a[1], colClasses = "character")   # index stays CHARACTER -> rule-of-2 counts distinct dates
sx <- read.csv(a[2], colClasses = "character")
pop <- readLines(a[3])

ph <- createPhenotypes(
  d,
  min.code.count = 2,            # rule-of-2: case = phecode-mapped codes on >= 2 distinct dates
  add.phecode.exclusions = TRUE, # PheWAS control exclusions (exclude related conditions from controls)
  translate = TRUE,              # map ICD -> phecode using the built-in Phecode Map 1.2
  id.sex = data.frame(id = sx$id, sex = sx$sex),
  full.population.ids = unique(pop)
)

write.csv(ph, a[4], row.names = FALSE)
cat("wrote", a[4], "-", nrow(ph), "participants x", ncol(ph) - 1, "phecodes\n")
