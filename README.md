# GEM Autocurator

**Drop a genome-scale metabolic model, get it curated — in your browser.**
Live: https://omidard.github.io/gem-autocurator/ · everything runs client-side; your model never leaves the page.

Upload an SBML or COBRA-JSON model and the autocurator:

1. **Identifier mapping** — resolves every metabolite & reaction id against **BiGG**, mints a deterministic **BiGG-like** id for compounds/reactions that only KEGG / ModelSEED / MetaNetX / ChEBI / RHEA carry, and *synthesises* a readable, BiGG-styled id from the name for entities absent from every database — so the model stays internally consistent and portable.
2. **Discrepancies** — the key step: detects when two different ids in your model refer to the *same* entity (duplicate metabolites/reactions) and merges them to one canonical id.
3. **Structural QC** — dead-end metabolites, orphan reactions.
4. **Mass & charge @ pH** — element conservation, and charge conservation with metabolite charges **recomputed at your simulation pH** (Henderson-Hasselbalch over ModelSEED pKa sites, anchored to the model's pH-7.2 charge); missing-proton/water fixes proposed and applied on export.
5. **Thermodynamics & cycles** — blocked reactions (FVA), energy-generating-cycle detection, and **ΔG-consensus directionality** (ModelSEED reversibility + Δ<sub>r</sub>G′ᵐ), all solved live with the bundled GLPK-WASM LP solver.
6. **Lab validation** — the organism is **auto-detected** from the model id/name (strain name, GCF/GCA assembly accession with or without prefix, or BV-BRC genome id; taxa resolved offline, accessions via the NCBI Datasets API — only the id string is sent, never the model). Its measured growth from **GrowthDB** (5,163 literature-cited records) and medium from the **Media DB** load into an **editable condition** — every exchange flux is tunable, rows from measured rates are tagged *experimental*, the rest *pre-set*. If there's no GrowthDB data, you enter your own μ and formulate the medium by hand. Then FBA — predicted μ<sub>max</sub> vs measured growth with an honest verdict and a secretion-pattern check.

You **supervise** every fix (approve / reject) and export a curated **SBML / JSON / MAT** (COBRA-Toolbox struct) plus a full analytical report with a Plotly dashboard.

## The identifier backbone — one stable id per compound/reaction

`tools/build_backbone.py` ingests the full **MetaNetX (MNXref)** cross-reference universe — which already clusters every compound and reaction from BiGG · KEGG · ModelSEED · ChEBI · MetaCyc · Rhea · HMDB by **InChIKey / structure** — and keeps the GEM-relevant clusters. Each cluster is assigned **one canonical id**: its BiGG id, or (when BiGG has none) **one stable BiGG-like id shared by every database id in the cluster**, frozen in `tools/backbone_freeze.json` so it never changes between builds. Every DB id, InChIKey and synonym then maps to that one id, so `glc__D`, `C00031`, `cpd00027` and `CHEBI:4167` all resolve to `glc__D` — and nobody can introduce a discrepancy with their own id logic. Coverage: **47,930 metabolite clusters** (8,885 in BiGG + 39,045 BiGG-like) and **68,294 reaction clusters** (20,731 in BiGG + 47,563 BiGG-like). Outputs are gzipped (`docs/data/*.json.gz`, ~4.6 MB total) and decompressed in-browser. Get the MetaNetX inputs from `https://www.metanetx.org/ftp/4.4/` into `/data/mnx_tmp`. (The older `tools/build_id_maps.py` is superseded.)
- `tools/build_thermo_pka.py` → `docs/data/{pka_table,thermo_rxn}.json` from the ModelSEED biochemistry DB (compound pKa sites; reaction reversibility + Δ<sub>r</sub>G′ᵐ).
- `tools/build_validation.py` → `docs/data/{growthdb,media_ex}.json` from GrowthDB + the Media DB.

Raw upstream inputs (ModelSEED TSVs, GrowthDB dump) are build-time only and kept outside the repo — see each script's header for the download commands.

## Stack
Pure client-side (GitHub Pages): SBML/JSON parsing, id canonicalisation, QC, pH charge titration, id synthesis, MAT/SBML/JSON export and the Plotly dashboard in vanilla JS; **GLPK-WASM** for every flux-based check (blocked reactions, EGCs, directionality, validation FBA). Design shared with [Flux Studio](https://omidard.github.io/FluxStudio/).
