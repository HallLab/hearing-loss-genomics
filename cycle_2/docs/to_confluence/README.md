# `to_confluence/` — pages staged for Confluence

Confluence is where this project is presented: plans, progress, results, and the decisions taken
at each meeting. This folder holds those pages **as source**, in the repo, so that what was shown
to the room and what the code actually did stay traceable to each other.

## Convention

**One page per step, numbered.** `NN_slug.md`, where `NN` is the order of the step in the cycle,
not the date. A step is a unit someone can approve or reject — "select the gene set", "build the
combined phenotype", "run Tier A". Meeting notes do **not** belong here; those stay in
`docs/meetings/`.

**Every page opens with `Decision requested`.** These pages exist to get an answer out of a
meeting. If a page has nothing to decide, it is a result write-up and belongs in
`cycle_2/results/` instead.

**English.** Same as `cycle_2/README.md` — the audience is Molly, Doug and Nikki.

**Numbers are generated, never transcribed.** Each page ends with a `Reproducibility` table naming
the scripts that produced its figures. When the underlying data is refreshed, regenerate the page
rather than editing numbers in place — that is the whole reason these live next to the code.

**State what is provisional.** A decision taken by Andre alone to keep work moving is marked
provisional with its date, so the room can tell it apart from something the group already ratified.

## Pasting into Confluence

Confluence Cloud converts pasted Markdown: headings, tables, code blocks and lists all survive.
Two things to watch:

- Paste into an **empty** page body. Pasting into an existing paragraph can drop the conversion.
- Blockquotes (`>`) and nested tables convert inconsistently — check them after pasting.

Keep the Confluence page title identical to the `# ` heading of the file so the two can be matched
later.

## Pages

| # | Page | Step | Status |
|---|---|---|---|
| 01 | [`01_gene_set_and_phecode_map.md`](01_gene_set_and_phecode_map.md) | Gene set selection + disease→phecode mapping | for review |

Page 01 is **generated**, not hand-written: regenerate with
`./venv/bin/python3 cycle_2/analysis/clingen_fetcher/render_confluence_page.py` after refreshing
the ClinGen snapshot. Edit the prose in that script, never in the `.md`.
