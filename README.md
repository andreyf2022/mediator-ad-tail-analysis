# Mediator activation-domain–Tail analysis: processed data and code

Processed source data and analysis code supporting quantitative analyses of activation-domain engagement with the human Mediator Tail.

Version 1.0.0 provides frozen processed Source Data plus deterministic validation, summary and plotting code supporting quantitative analyses of activation-domain engagement with the human Mediator Tail in the accompanying Mediator manuscript. It is not an end-to-end analysis repository: upstream AF3 processing, contact calling, deletion mapping/classification, bootstrap resampling, STRIDE/pLDDT mapping and structural-alignment calculations are not rerun here. It uses processed inputs and does not contain AF3 coordinates, PAE/full-data files, MSAs, templates, structural sessions, Illustrator files, upstream structural-analysis code or manual structural-rendering machinery.

## Supported results

- Multisubunit model and recurrent-AD statistics: validates and redraws the exact released counts table.
- Deletion-centered contact profiles and WFYL-zero proportions: redraws the released analyzed source tables, applying only the stated five-residue display smoothing to the contact profiles.
- Interpeak interval distributions: recomputes descriptive n/median/IQR summaries from the frozen interval table and redraws the panel.
- pLDDT/helicity summaries: redraws frozen estimates/intervals and recomputes only the displayed Benjamini-Hochberg adjustment from the supplied contrast P values.
- Receptor-site co-engagement network: redraws frozen node values, all 27 plotted nonzero edges and the exact eight projected node positions.
- The reported MED23/MED24 numerical geometry values as a frozen coordinate-free summary. This package does not recalculate them.
- The 293-job AF3-versus-8GXS per-job RMSDs and reported median/IQR as frozen coordinate-free Source Data. This package does not recalculate them.
- Extended Data Table 1 as a machine-readable extraction of the frozen final table.
- Exact construct sequences, seeds, model IDs, ranking scores, model-retention state, source intervals, noncontiguous-junction annotations, local-mask assignments, and receptor-site envelopes defining the analyzed ensemble.

## Environment and rendering

Create the environment:

```bash
conda env create -f config/environment.yml
conda activate mediator_coordinate_free_release
```

The environment uses `conda-forge` with `nodefaults`, so user-configured Conda channels are not inherited.

From the package root, render the coordinate-free quantitative panels into new output paths:

```bash
python code/render_figure1_multisubunit.py --output outputs/figure1_multisubunit.png
python code/render_figure2_panels_CD.py --output-dir outputs/figure2_panels_CD
python code/render_extended_data_figure1.py --output-dir outputs/extended_data_figure1
python code/render_extended_data_figure2.py --output-dir outputs/extended_data_figure2
python code/render_figure3b_network.py --out-dir outputs/figure3b_network
```

All scripts expose `--help` and refuse to overwrite existing output files or directories. Their default input paths are package-relative and can be replaced with documented arguments. The retained renderers reproduce their canonical repository render outputs byte-for-byte in the tested environment.

## What the code does

- No retained script performs upstream scientific analysis from coordinates or raw measurements.
- `render_extended_data_figure1.py` performs a small descriptive aggregation (n, median and IQR) before plotting.
- `render_extended_data_figure2.py` recomputes only the displayed Benjamini-Hochberg adjustment across the 36 supplied contrast P values before plotting.
- `render_figure1_multisubunit.py` validates frozen summary/per-AD records and plots them.
- `render_figure2_panels_CD.py` and `render_figure3b_network.py` are plotting scripts; Figure 2C additionally applies its documented five-residue display smoothing.

## Dataset definition

`data/dataset/construct_model_manifest.csv` replaces 293 job-request JSONs and 1,465 summary-confidence JSONs. It contains the release-relevant information from those files: exact modeled AD sequence, AlphaFold Server job name and seed, model IDs 0-4, all five `ranking_score` values, the model chosen for the 8GXS comparison, and the final contact-ensemble status of every model. Its validated totals are 293 constructs, 1,465 models, 1,139 retained models and 234 represented AD constructs; 22 constructs contain noncontiguous concatenated source intervals.

The initial and postmask whole-model exclusions are encoded per model in that manifest. Exact partial-inpainting and local-threading contact masks are supplied in `inpainting_operational_ledger.csv` and `threading_local_mask_residues.csv`. The authoritative reviewed receptor envelopes are in `receptor_site_definitions.csv`.

The three receptor-site co-engagement tables were exported by one authorized execution of the unchanged scientific calculation, with only a staging-local CSV export hook. The run reproduced the canonical PNG exactly (SHA-256 `2f806ecccf42d52a29783ee6679840701040dbcd047a86e9fb33b89d1c828e14`). The released renderer reads those tables without coordinates and also reproduces that PNG byte-for-byte. The frozen graph preserves the released eight-site set, including quantitative treatment of MED24 D653.

## Coordinate free limitations and blockers

- No AF3-versus-8GXS executable is included. The audited scientific procedure is known, but the exact semiglobal sequence-matching implementation was not persisted; reconstructing it would require guessing. The per-job and summary values are included as Source Data.
- The reported MED23/MED24 geometry values are supplied as a frozen coordinate-free summary. The exact matched-atom and transformation outputs were not persisted, so atom-level reproduction requires the private coordinate/session archive.
- Manual structural panels and the conceptual model in the accompanying Mediator manuscript are outside this coordinate-free quantitative package.

## Frozen scientific scope

The package preserves the frozen datasets, definitions and scientific choices encoded in the released tables, including the distinct deletion populations used for the contact-profile and WFYL summaries and the receptor-network inclusion of MED24 D653.

`release_manifest.csv` records every payload file's stable release provenance, SHA-256 and purpose; the manifest omits only its own self-hash. Derived compact files name all frozen inputs used to create them.

## Authors and affiliation

- Andrey Feklistov, Department of Structural Biology, Stanford University School of Medicine, Stanford, California 94305, USA; ORCID 0000-0001-6207-1160; correspondence: andreyf@stanford.edu.
- Roger D. Kornberg, Department of Structural Biology, Stanford University School of Medicine, Stanford, California 94305, USA; ORCID 0000-0002-2425-7519.

## Citation

Release citation metadata are provided in `CITATION.cff`. A DOI for the accompanying Mediator manuscript can be added after publication without changing the scientific payload.

## Licenses

- Code under `code/` and the environment definition under `config/`: MIT License (`LICENSE-CODE`).
- Processed data under `data/`, README and metadata: Creative Commons Attribution 4.0 International (`LICENSE-DATA`).
