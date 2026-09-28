#!/usr/bin/env python3
"""
map_clingen_diseases_to_phecodes.py — list the diseases attached to the adopted
ClinGen HL gene set, and relate them to PheWAS phecodes.

Runs on the output of fetch_clingen_hl_genes.py. Two products:

  1. The disease list — straightforward. ClinGen already gives one curated
     DISEASE LABEL + MONDO id per gene-disease pair.

  2. The phecode relation — an inference chain, not a lookup:

         MONDO term --xref--> ICD-9-CM / ICD-10-CM --phecode map--> phecode

     Every link loses coverage. MONDO carries an ICD xref for only a minority
     of terms (~2.1K ICD10CM xrefs over 63K terms), because most Mendelian
     syndromes have no billing code of their own. Where the term itself has no
     xref we walk up `is_a` to the nearest ancestor that does, and record how
     far we walked. **Read the MATCH_LEVEL column before using a row**: an
     ancestor hop of 2-3 usually means "some broad category", not this disease.

     This is a triage aid for phenotype design, not an equivalence table. The
     expected finding is that most of these diseases do NOT have a phecode --
     which is itself the argument for Cycle 2's Tier B audiometric phenotype
     (`cycle_2/README.md` §3).

Sources:
  MONDO    https://purl.obolibrary.org/obo/mondo.obo         (cached in --outdir)
  phecodes analysis/chapter_2_v2/results/phecode_map12.csv   (Phecode 1.2 ICD map)
  labels   PheWAS R package `pheinfo`                        (same table createPhenotypes uses)
"""

import argparse
import csv
import json
import re
import subprocess
import sys
import urllib.request
from collections import defaultdict, deque
from pathlib import Path

MONDO_URL = "https://purl.obolibrary.org/obo/mondo.obo"
PHEINFO_RDA_URL = "https://raw.githubusercontent.com/PheWAS/PheWAS/master/data/pheinfo.rda"
ICD_VOCABS = ("ICD10CM", "ICD9CM", "ICD9")
MAX_ANCESTOR_DEPTH = 3

# MONDO's is_a hierarchy mixes disease categories with *mode of inheritance* and
# top-level grouping terms. Several of them carry an ICD xref of their own
# ("autosomal dominant disease" -> ICD9:758.5, a chromosomal-anomaly code), so an
# unguarded ancestor walk drags a third of the gene set onto phecode 758.1. These
# terms say nothing about what the disease *is* — never resolve through them.
EXCLUDED_ANCESTORS = {
    "MONDO:0000001",  # disease
    "MONDO:0000425",  # X-linked disease
    "MONDO:0000426",  # autosomal dominant disease
    "MONDO:0000429",  # autosomal genetic disease
    "MONDO:0002254",  # syndromic disease
    "MONDO:0003847",  # hereditary disease
    "MONDO:0006025",  # autosomal recessive disease
    "MONDO:0019040",  # chromosomal disorder
    "MONDO:0020604",  # X-linked dominant disease
    "MONDO:0020605",  # X-linked recessive disease
    "MONDO:0700096",  # human disease
}

# --- confidence triage -------------------------------------------------------
# Not every phecode a disease reaches is worth quoting. Three buckets, assigned
# by rule so the classification is reviewable and moves with the data:
#
#   clinical    — a real manifestation of the syndrome. Usable.
#   wastebasket — the mapping is CORRECT but the phecode is a residual
#                 "other/unspecified" bin, or the MONDO ancestor is a
#                 biochemical/system category rather than a presentation.
#                 Nobody is coded this way in an EHR for this reason.
#   suspect     — do not quote without checking. Either the ancestor walk went
#                 too far, or the underlying xref is known to be wrong.
#
# `suspect` means "needs review", not "wrong" — rule 1 has false positives
# (Brown-Vialetto-van Laere really is a motor neuron disorder at depth 3).

# MONDO xrefs that are demonstrably incorrect. Curated, with the reason, because
# no heuristic detects them: they are specific codes, at shallow depth, that
# simply point at the wrong organ.
KNOWN_BAD_XREF = {
    # MONDO xrefs "primary ovarian failure" to ICD9:253.4 = disorders of the
    # ANTERIOR PITUITARY. Primary ovarian failure is gonadal, not pituitary.
    # The correct phecode for these genes (Perrault syndrome) is 627.5, which
    # the same MONDO term also yields.
    ("primary ovarian failure", "253.4"),
}

# Residual phecode bins: correct, uninformative.
_WASTEBASKET_DESC = re.compile(r"^Other and unspecified|^Unspecified|\bNEC\b")
# MONDO ancestors that classify by biochemistry or organ system rather than by
# clinical presentation ("inherited lipid metabolism disorder"), as opposed to a
# named disease ("Bartter syndrome", "long QT syndrome").
_CATEGORY_ANCESTOR = re.compile(
    r"metabolism disorder$|^inborn disorder of|system disorder$", re.I)


def confidence(phecode: str, desc: str, match_level: str, via: str) -> str:
    if (via, phecode) in KNOWN_BAD_XREF:
        return "suspect"
    if _WASTEBASKET_DESC.search(desc) or (via and _CATEGORY_ANCESTOR.search(via)):
        return "wastebasket"
    if match_level == "ancestor_d3":
        return "suspect"
    return "clinical"


_CONF_RANK = {"clinical": 0, "wastebasket": 1, "suspect": 2}


def repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").exists():
            return parent
    return Path.cwd()


def download(url: str, dest: Path) -> None:
    print(f"[fetch] {url}", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=600) as resp:
        dest.write_bytes(resp.read())
    print(f"[fetch] wrote {dest} ({dest.stat().st_size:,} bytes)", file=sys.stderr)


def parse_obo(path: Path) -> dict[str, dict]:
    """Minimal OBO reader — id, name, xrefs, is_a. Skips obsolete terms."""
    terms: dict[str, dict] = {}
    cur: dict | None = None
    in_term = False

    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("["):
                if cur and in_term and not cur["obsolete"]:
                    terms[cur["id"]] = cur
                in_term = line == "[Term]"
                cur = {"id": "", "name": "", "xref": defaultdict(list),
                       "is_a": [], "obsolete": False} if in_term else None
                continue
            if not in_term or cur is None or ": " not in line:
                continue
            tag, _, val = line.partition(": ")
            # OBO trailing qualifiers: `xref: ICD10CM:H90.3 {source="..."}`
            val = val.split(" {", 1)[0].strip()
            if tag == "id":
                cur["id"] = val
            elif tag == "name":
                cur["name"] = val
            elif tag == "is_obsolete" and val == "true":
                cur["obsolete"] = True
            elif tag == "is_a":
                cur["is_a"].append(val.split(" !", 1)[0].strip())
            elif tag == "xref":
                vocab, _, code = val.partition(":")
                if code:
                    cur["xref"][vocab].append(code.strip())
    if cur and in_term and not cur["obsolete"]:
        terms[cur["id"]] = cur
    return terms


def load_phecode_map(path: Path) -> dict[tuple[str, str], set[str]]:
    """(vocabulary, ICD code) -> {phecode}. MONDO says ICD9, phecode map says ICD9CM."""
    out: dict[tuple[str, str], set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            vocab, code, phecode = (row["vocabulary_id"].strip(),
                                    row["code"].strip(), row["phecode"].strip())
            out[(vocab, code)].add(phecode)
            if vocab == "ICD9CM":
                out[("ICD9", code)].add(phecode)
    return out


def ensure_pheinfo(csv_path: Path, rda_path: Path) -> dict[str, tuple[str, str]]:
    """phecode -> (description, group). Converts the PheWAS .rda via R if needed."""
    if not csv_path.exists():
        if not rda_path.exists():
            download(PHEINFO_RDA_URL, rda_path)
        r = ('e <- new.env(); load("%s", envir=e); '
             'write.csv(get(ls(e)[1], envir=e), "%s", row.names=FALSE)'
             % (rda_path, csv_path))
        print("[pheinfo] converting .rda via Rscript", file=sys.stderr)
        try:
            subprocess.run(["Rscript", "-e", r], check=True,
                           capture_output=True, text=True)
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            print(f"[pheinfo] WARNING: no phecode labels ({exc}).\n"
                  f"          Run `module load R/4.5` first, or pass --pheinfo-csv.",
                  file=sys.stderr)
            return {}
    with csv_path.open(newline="", encoding="utf-8") as f:
        return {r["phecode"].strip(): (r["description"], r.get("group", ""))
                for r in csv.DictReader(f)}


def icd_to_phecodes(vocab: str, code: str,
                    pmap: dict) -> list[tuple[str, str, str]]:
    """(vocab, code) -> [(phecode, matched_code, how)]. Falls back to the parent
    ICD code (H90.3 -> H90), which is a real widening — flagged as such."""
    hits = pmap.get((vocab, code))
    if hits:
        return [(p, code, "exact") for p in sorted(hits)]
    if "." in code:
        parent = code.split(".", 1)[0]
        hits = pmap.get((vocab, parent))
        if hits:
            return [(p, parent, "icd_parent") for p in sorted(hits)]
    return []


def resolve_term(mondo_id: str, terms: dict, pmap: dict
                 ) -> tuple[list[dict], str]:
    """Find phecodes for a MONDO term, walking up is_a only if it has no ICD
    xref of its own. Returns (rows, note)."""
    if mondo_id not in terms:
        return [], "mondo_id_not_in_ontology"

    def direct(tid: str) -> list[dict]:
        rows = []
        for vocab in ICD_VOCABS:
            for code in terms[tid]["xref"].get(vocab, []):
                for phecode, matched, how in icd_to_phecodes(vocab, code, pmap):
                    rows.append({"icd_vocab": vocab, "icd_code": code,
                                 "icd_matched": matched, "phecode": phecode,
                                 "via_mondo": tid, "icd_match": how})
        return rows

    rows = direct(mondo_id)
    if rows:
        for r in rows:
            r["match_level"] = "direct"
        return rows, ""

    # Breadth-first up the is_a graph. Depth is reported, not hidden: an
    # ancestor 3 hops up is a category, not this disease.
    seen = {mondo_id}
    frontier = deque((p, 1) for p in terms[mondo_id]["is_a"])
    while frontier:
        tid, depth = frontier.popleft()
        if tid in seen or depth > MAX_ANCESTOR_DEPTH or tid not in terms:
            continue
        if tid in EXCLUDED_ANCESTORS:
            continue  # do not resolve through it, and do not climb past it
        seen.add(tid)
        rows = direct(tid)
        if rows:
            for r in rows:
                r["match_level"] = f"ancestor_d{depth}"
                r["via_mondo_name"] = terms[tid]["name"]
            return rows, ""
        frontier.extend((p, depth + 1) for p in terms[tid]["is_a"])
    return [], f"no_icd_xref_within_{MAX_ANCESTOR_DEPTH}_ancestors"


def main() -> int:
    repo = repo_root()
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, default=repo / "cycle_2/data/clingen")
    ap.add_argument("--pairs-tsv", type=Path, default=None,
                    help="default: <outdir>/clingen_hl_gcep.tsv")
    ap.add_argument("--gene-list", type=Path, default=None,
                    help="restrict to these genes (default: <outdir>/genes_definitive_strong.txt); "
                         "pass 'all' for every HL GCEP gene")
    ap.add_argument("--phecode-map", type=Path,
                    default=repo / "analysis/chapter_2_v2/results/phecode_map12.csv")
    ap.add_argument("--pheinfo-csv", type=Path, default=None)
    ap.add_argument("--mondo-obo", type=Path, default=None)
    args = ap.parse_args()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    pairs_tsv = args.pairs_tsv or outdir / "clingen_hl_gcep.tsv"
    mondo_obo = args.mondo_obo or outdir / "mondo.obo"
    pheinfo_csv = args.pheinfo_csv or outdir / "phecode_definitions_pheinfo.csv"

    if not pairs_tsv.exists():
        sys.exit(f"{pairs_tsv} not found — run fetch_clingen_hl_genes.py first")
    if not args.phecode_map.exists():
        sys.exit(f"phecode map not found: {args.phecode_map}")
    if not mondo_obo.exists():
        download(MONDO_URL, mondo_obo)

    genes_keep: set[str] | None = None
    gl = args.gene_list or outdir / "genes_definitive_strong.txt"
    if str(gl) != "all":
        if not Path(gl).exists():
            sys.exit(f"gene list not found: {gl}")
        genes_keep = {g.strip() for g in Path(gl).read_text().split() if g.strip()}

    with pairs_tsv.open(newline="", encoding="utf-8") as f:
        pairs = [r for r in csv.DictReader(f, delimiter="\t")]
    if genes_keep is not None:
        pairs = [r for r in pairs if r["GENE SYMBOL"] in genes_keep]

    print(f"[load] {len(pairs)} gene-disease pairs"
          f"{f' over {len(genes_keep)} genes' if genes_keep else ''}", file=sys.stderr)
    print("[load] parsing mondo.obo ...", file=sys.stderr)
    terms = parse_obo(mondo_obo)
    print(f"[load] {len(terms):,} non-obsolete MONDO terms", file=sys.stderr)
    pmap = load_phecode_map(args.phecode_map)
    pheinfo = ensure_pheinfo(pheinfo_csv, outdir / "pheinfo.rda")

    # 1 — disease list
    diseases_tsv = outdir / "hl_diseases.tsv"
    dcols = ["GENE SYMBOL", "DISEASE LABEL", "DISEASE ID (MONDO)", "MOI",
             "CLASSIFICATION", "MONDO_NAME", "OMIM", "ORPHANET"]
    with diseases_tsv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(dcols)
        for r in sorted(pairs, key=lambda r: (r["GENE SYMBOL"], r["DISEASE LABEL"])):
            t = terms.get(r["DISEASE ID (MONDO)"], {})
            xref = t.get("xref", {})
            w.writerow([r["GENE SYMBOL"], r["DISEASE LABEL"], r["DISEASE ID (MONDO)"],
                        r["MOI"], r["CLASSIFICATION"], t.get("name", ""),
                        ";".join(xref.get("OMIM", [])),
                        ";".join(xref.get("Orphanet", []))])

    # 2 — disease -> phecode
    resolved: dict[str, tuple[list[dict], str]] = {}
    map_rows, unmapped = [], []
    for r in pairs:
        mid = r["DISEASE ID (MONDO)"]
        if mid not in resolved:
            resolved[mid] = resolve_term(mid, terms, pmap)
        hits, note = resolved[mid]
        if not hits:
            unmapped.append([r["GENE SYMBOL"], r["DISEASE LABEL"], mid,
                             r["CLASSIFICATION"], note or "no_phecode"])
            continue
        # One MONDO term often xrefs several ICD codes that collapse onto the
        # same phecode (H90, 389, 389.8, 389.9 -> 389). Emit the phecode once
        # and keep the codes as evidence, rather than four near-identical rows.
        collapsed: dict[str, dict] = {}
        for h in hits:
            e = collapsed.setdefault(h["phecode"], {"icd": [], "level": h["match_level"],
                                                    "via": h.get("via_mondo_name", ""),
                                                    "match": set()})
            e["icd"].append(f'{h["icd_vocab"]}:{h["icd_code"]}')
            e["match"].add(h["icd_match"])
        for phecode, e in collapsed.items():
            desc, group = pheinfo.get(phecode, ("", ""))
            map_rows.append([r["GENE SYMBOL"], r["DISEASE LABEL"], mid,
                             r["MOI"], r["CLASSIFICATION"],
                             phecode, desc, group,
                             confidence(phecode, desc, e["level"], e["via"]),
                             ";".join(sorted(set(e["icd"]))),
                             ";".join(sorted(e["match"])),
                             e["level"], e["via"]])

    map_tsv = outdir / "hl_disease_phecode_map.tsv"
    with map_tsv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["GENE", "DISEASE_LABEL", "MONDO", "MOI", "CLASSIFICATION",
                    "PHECODE", "PHECODE_DESC", "PHECODE_GROUP", "CONFIDENCE",
                    "ICD_EVIDENCE", "ICD_MATCH", "MATCH_LEVEL", "VIA_ANCESTOR"])
        w.writerows(sorted(map_rows, key=lambda r: (r[0], r[5])))

    unmapped_tsv = outdir / "hl_diseases_unmapped.tsv"
    with unmapped_tsv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["GENE", "DISEASE_LABEL", "MONDO", "CLASSIFICATION", "REASON"])
        w.writerows(sorted(unmapped))

    # 3 — phecode-centric view: which phecodes does this gene set implicate?
    by_phecode: dict[str, dict] = defaultdict(
        lambda: {"genes": set(), "diseases": set(), "levels": set(), "conf": set()})
    for row in map_rows:
        e = by_phecode[row[5]]
        e["genes"].add(row[0]); e["diseases"].add(row[1]); e["levels"].add(row[11]); e["conf"].add(row[8])
    summary_tsv = outdir / "hl_phecode_summary.tsv"
    with summary_tsv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["PHECODE", "PHECODE_DESC", "PHECODE_GROUP", "CONFIDENCE",
                    "N_GENES", "N_DISEASES", "MATCH_LEVELS", "GENES"])
        for p, e in sorted(by_phecode.items(),
                           key=lambda kv: (_CONF_RANK[min(kv[1]["conf"], key=_CONF_RANK.get)],
                                           -len(kv[1]["genes"]), kv[0])):
            desc, group = pheinfo.get(p, ("", ""))
            best = min(e["conf"], key=_CONF_RANK.get)
            w.writerow([p, desc, group, best, len(e["genes"]), len(e["diseases"]),
                        ";".join(sorted(e["levels"])), ";".join(sorted(e["genes"]))])

    n_pairs_mapped = len({(r[0], r[2]) for r in map_rows})
    n_pairs_total = len({(r["GENE SYMBOL"], r["DISEASE ID (MONDO)"]) for r in pairs})
    direct = len({(r[0], r[2]) for r in map_rows if r[11] == "direct"})
    stats = {
        "gene_list": str(gl), "n_pairs": len(pairs),
        "n_unique_gene_disease": n_pairs_total,
        "n_mapped_to_any_phecode": n_pairs_mapped,
        "n_mapped_direct_xref": direct,
        "n_mapped_via_ancestor": n_pairs_mapped - direct,
        "n_unmapped": n_pairs_total - n_pairs_mapped,
        "n_distinct_phecodes": len(by_phecode),
        "phecodes_by_confidence": {
            c: sum(1 for e in by_phecode.values()
                   if min(e["conf"], key=_CONF_RANK.get) == c)
            for c in ("clinical", "wastebasket", "suspect")},
        "mondo_terms_parsed": len(terms),
        "phecode_map": str(args.phecode_map),
        "phecode_labels": "PheWAS::pheinfo" if pheinfo else "UNAVAILABLE",
    }
    (outdir / "phecode_mapping_stats.json").write_text(
        json.dumps(stats, indent=2) + "\n", encoding="utf-8")

    print(f"\ngene-disease pairs        : {n_pairs_total}")
    print(f"  mapped to >=1 phecode   : {n_pairs_mapped}"
          f"  ({direct} direct xref, {n_pairs_mapped - direct} via ancestor)")
    print(f"  unmapped                : {n_pairs_total - n_pairs_mapped}")
    conf_counts = collections_counter = {
        c: sum(1 for e in by_phecode.values()
               if min(e["conf"], key=_CONF_RANK.get) == c)
        for c in ("clinical", "wastebasket", "suspect")}
    print(f"distinct phecodes reached : {len(by_phecode)}"
          f"  ({conf_counts['clinical']} clinical,"
          f" {conf_counts['wastebasket']} wastebasket,"
          f" {conf_counts['suspect']} suspect)")
    print(f"\nwrote:\n  {diseases_tsv.name}\n  {map_tsv.name}"
          f"\n  {summary_tsv.name}\n  {unmapped_tsv.name}\n  phecode_mapping_stats.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
