# Processed data and plotting scripts accompanying *A Structural Model of Transcriptional Activation by Mediator*

This repository contains processed source data, consistency checks, and five Python renderers supporting quantitative analyses of activation-domain (AD) engagement with the human Mediator Tail. It is a coordinate-free companion to the manuscript, not an end-to-end analysis pipeline.

The deposited ensemble comprises 293 AlphaFold Server jobs (five models each), with 1,139 models from 234 AD constructs retained for contact analyses after whole-model artifact exclusions and residue-level masks. The deposit excludes AF3 coordinates, PAE/full-data files, MSAs, templates, PyMOL sessions, upstream contact/filtering/statistical code, and manual structural-rendering files.

## Contents and supported outputs

| Result | Processed input | Script | What is recomputed |
|---|---|---|---|
| Figure 1B/C quantitative multisubunit summaries | `data/figure1/final_counts_used.csv`; construct manifest | `code/render_figure1_multisubunit.py` | Internal count/fraction checks and a summary plot; not the assembled donut layout |
| Figure 2C/D | `data/figure2/*.csv` | `code/render_figure2_panels_CD.py` | Five-position display smoothing for panel C; charts from stored estimates and intervals |
| Extended Data Figure 1 | `data/extended_data_figure1/*` | `code/render_extended_data_figure1.py` | n, median, and IQR from identified intervals; distribution plot |
| Extended Data Figure 2 | `data/extended_data_figure2/*` | `code/render_extended_data_figure2.py` | Combined Benjamini–Hochberg adjustment across the 36 supplied P values; charts |
| Figure 3B | `data/figure3/figure3b_{nodes,edges,positions}.csv` | `code/render_figure3b_network.py` | Schematic network from stored node values, 27 nonzero edges, and eight fixed positions |

`data/figure3/figure3c_reported_geometry.csv`, `data/rmsd/`, and `data/extended_data_table1/` are coordinate-free source-data summaries. They cannot regenerate the underlying structural alignments. `data/dataset/construct_model_manifest.csv` retains selected job-request, confidence, and filtering fields; it is not a lossless substitute for the original AF3 inputs or outputs. See [`data/README.md`](data/README.md) for definitions and joins.

## Install, validate, and render

Tested versions are pinned in `config/environment.yml`:

```bash
conda env create -f config/environment.yml
conda activate mediator-ad-tail-analysis
python tests/validate_release.py
```

From the repository root, write outputs to new paths:

```bash
python code/render_figure1_multisubunit.py --output outputs/figure1_multisubunit.png
python code/render_figure2_panels_CD.py --output-dir outputs/figure2_panels_CD
python code/render_extended_data_figure1.py --output-dir outputs/extended_data_figure1
python code/render_extended_data_figure2.py --output-dir outputs/extended_data_figure2
python code/render_figure3b_network.py --out-dir outputs/figure3b_network
```

Every renderer exposes `--help`, uses package-relative default inputs, and refuses to overwrite an existing output. DejaVu Sans is used directly or as the portable fallback when Arial/Helvetica is unavailable. Numerical checks are deterministic in the pinned environment; exact image bytes can vary with fonts or rendering libraries.

## Scope limits

The scripts validate, aggregate, or redraw deposited processed tables. They do not rerun AlphaFold, call contacts from coordinates, repeat deletion mapping or bootstrap resampling, calculate STRIDE/pLDDT features, reproduce structural illustrations, or recalculate the AF3-versus-8GXS RMSDs and MED23/MED24 geometry. The exact semiglobal sequence-matching executable used for the RMSDs was not persisted.

## Citation and licenses

Citation metadata are in `CITATION.cff`; the original v1.0.0 archive is [Zenodo DOI 10.5281/zenodo.22819859](https://doi.org/10.5281/zenodo.22819859). Version 1.0.1 is an unreleased maintenance update until a new archive identifier is assigned.

Code and the environment definition are licensed under MIT (`LICENSE-CODE`). Processed data and documentation are licensed under CC BY 4.0 (`LICENSE-DATA`).
