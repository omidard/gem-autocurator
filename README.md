# GEM Autocurator

**Drop a genome-scale metabolic model, get it curated — in your browser.**
Live: https://omidard.github.io/gem-autocurator/ · everything runs client-side; your model never leaves the page.

Upload an SBML or COBRA-JSON model and the autocurator:

1. **Identifier mapping** — resolves every metabolite & reaction id against **BiGG**, and mints a deterministic **BiGG-like** id for compounds/reactions that only KEGG / ModelSEED / MetaNetX / ChEBI / RHEA carry, so the model stays internally consistent and portable.
2. **Discrepancies** — the key step: detects when two different ids in your model refer to the *same* entity (duplicate metabolites/reactions) and merges them to one canonical id.
3. **Structural QC** — dead-end metabolites, orphan reactions.
4. **Mass & charge @ pH** — element and charge conservation per reaction, at your chosen simulation pH.
5. **Thermodynamics & cycles** *(v2)* — blocked reactions (FVA), energy-generating-cycle detection, ΔG-consensus directionality.
6. **Lab validation** *(v2)* — fuzzy strain match against **GrowthDB**, media from **MediaDB**, GAM/NGAM fitting, prediction-vs-experiment.

You **supervise** every fix (approve / reject) and export a curated **SBML / JSON** (MAT roadmap) plus a full analytical report.

## Reference backbone
`tools/build_id_maps.py` builds `docs/data/{metabolite_map,reaction_map,bigg_met_props}.json` from BiGG + MetaNetX MNXref + KEGG/ModelSEED/ChEBI cross-references (9,403 BiGG + 35k BiGG-like metabolites; 28k BiGG + 124k BiGG-like reactions). Regenerate with `python3 tools/build_id_maps.py`.

## Stack
Pure client-side (GitHub Pages): SBML/JSON parsing, id canonicalisation, QC and export in vanilla JS; GLPK-WASM for the flux-based checks. Design shared with [Flux Studio](https://omidard.github.io/FluxStudio/).
