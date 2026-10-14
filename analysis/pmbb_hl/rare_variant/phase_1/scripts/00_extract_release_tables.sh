#!/usr/bin/env bash
# Phase 1 / check 00 — cache the two release tables the later checks need.
#
# Both source files are multi-GB and are scanned once each, so the extracts are
# cached under phase_1/data/ rather than re-read. They are regenerable and untracked.
#
#   _ear_family.tsv    all SO_39x (ear-family) diagnosis rows from conditions_phecode_x
#   _tinnitus_obs.tsv  tinnitus rows from the OMOP observation table
#
# The second one exists because PMBB v4 relocated the standard tinnitus ICD codes
# (388.3x, H93.1x) into `observation`. The phecode files keep only the ~3,016
# pulsatile events; the bulk (25,094) is here. Omitting this source makes the
# ear-family evidence incomplete — see results/FINDINGS.md, sub-finding A.
#
# Expected output (PMBB-Release-2026-4.0):
#   _ear_family.tsv    655,946 rows
#   _tinnitus_obs.tsv   25,094 rows

set -euo pipefail

REL=/static/PMBB/PMBB-Release-2026-4.0/Phenotype/4.0
PFX=PMBB-Release-2026-4.0_phenotype
OUT="$(cd "$(dirname "$0")/.." && pwd)/data"
mkdir -p "$OUT"

echo "==> ear-family diagnoses (SO_39x) from conditions_phecode_x"
awk -F'\t' '
  NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
  $h["condition_source_value"] ~ /^SO_39[0-9]($|\.)/ {
    print $h["person_id"]"\t"$h["condition_start_date"]"\t"$h["condition_source_value"]
  }' "$REL/${PFX}_conditions_phecode_x.txt" > "$OUT/_ear_family.tsv"
echo "    $(wc -l < "$OUT/_ear_family.tsv") rows  (expected 655,946)"

echo "==> tinnitus from the OMOP observation table (388.3x, H93.1x)"
awk -F'\t' '
  NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
  $h["observation_source_value"] ~ /^(388\.3|H93\.1)/ {
    print $h["person_id"]"\t"$h["observation_date"]"\t"$h["observation_source_value"]
  }' "$REL/${PFX}_observation.txt" > "$OUT/_tinnitus_obs.tsv"
echo "    $(wc -l < "$OUT/_tinnitus_obs.tsv") rows  (expected 25,094)"

echo "done -> $OUT"
