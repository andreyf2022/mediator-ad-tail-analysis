# Data dictionary

All residue positions are one-based and intervals are inclusive unless stated otherwise. Blank fields mean not applicable or not recorded; they do not imply a zero measurement or a reviewed negative result.

## Dataset and filtering tables

### `dataset/construct_model_manifest.csv`

One row per modeled AD construct (293 rows; `construct_index` 2–294). `gene` is the lowercase join key; `display_gene` is the label; `af3_job_name` preserves the submitted name (including the E2F1 and ZBTB44 `_RERUN` suffixes); `fold_dir` is encoded within the per-model status IDs and must not be derived by lowercasing the job name. A full model ID is `fold_dir__model_N`, where N is 0–4. `ranking_score_model_N` is the AlphaFold Server ranking score, not a contact score.

`modeled_ad_sequence` is the exact chain-H sequence and `ad_length` its residue count. `source_components` summarizes source-coordinate annotations. Those annotations are not guaranteed to have the same length as the experimentally supplied tile sequence; `dataset/construct_source_components.csv` records the actual assembly. Stable indices, sequences, seeds, scores, and selection states are unchanged from v1.0.0.

Production files identify the platform build as `_software.version` and `_ma_model_list.model_group_name` = `AlphaFold-beta-20231127 (7a14c934-fe78-43c1-a552-7fb305b0ebd1 @ 2024-11-04 18:04:41)`. The shorter build label is repeated in the manifest.

Each `model_N_contact_status` is `RETAINED`, an initial heavy-atom-overlap exclusion, or a subsequent whole-model exclusion. For the initial screen, `n_AD_heavy_atoms_within_1p2A` counts AD heavy atoms whose nearest receptor heavy atom was at most 1.2 Å; it is not an AF3 confidence value. The manifest contains 86 initial exclusions, 240 later whole-model exclusions, and 1,139 retained models from 234 constructs.

`model_N_inpainting_class=NOT_LISTED_IN_REVIEW_LEDGER` means only that the model is absent from the special inpainting operational ledger; it does not assert either human approval or lack of other review. `n_retained_contact_models` and `retained_contact_model_ids` are construct-level summaries.

### `dataset/construct_source_components.csv`

One row per original activating tile (374 rows). Rows are keyed by `construct_index`, `gene`, and `component_order`. Tiles were ordered by source start/end/domain; the first tile was copied completely, and for each subsequent overlapping tile the already represented N-terminal overlap was removed before appending. Separated tiles were appended directly, omitting native intervening sequence. `synthetic_junction_after_previous_component=yes` marks those direct noncontiguous joins (22 constructs). Concatenating `appended_sequence` by construct exactly reconstructs every modeled sequence.

### `dataset/inpainting_operational_ledger.csv`

One row per model listed in the special inpainting ledger. `PARTIAL_CONTACT_MASK` supplies a chain-H range for models that remain in contact analyses; ranges in `masked_region_note` are inclusive and total 913 masked residues across the 21 retained models. Eight additional partial-mask ledger rows belong to models already removed by the initial whole-model screen, so their ranges do not contribute to retained contact calculations. `FULLY_EXCLUDE_PRIMARY_DATASET` and `INTERFACE_EXCLUDE_KEEP_SS` are whole-model exclusions for contact analyses. `NOT_INPAINTED` is an explicit ledger class.

### `dataset/threading_local_mask_residues.csv`

One row per masked chain-H residue (61 rows), keyed by model ID and `ad_resnum`. These are residue-level contact masks, not whole-model exclusions.

### `dataset/receptor_site_definitions.csv`

One row per final analysis site/candidate. `receptor_site_id` is the stable join identifier; `figure3b_hotspot_key` is the stable graph key. Receptor coordinates are one-based construct-local AF3 positions. `display_marker` uses native protein numbering where available; notably MED15 A:76 maps to native L595 in the MED15 520–653 construct. Core plus shoulder is the contact-counting envelope. General shoulders used within-site support ≥0.50 or global support ≥0.70; the MED16/MED25 composite used the recorded stricter ≥0.65/≥0.80 exception. Eight rows are included in Figure 3B. MED23 G986 is retained as a provenance candidate but not displayed. MED24 D653 is included without downweighting in the network; its historical haze/downweighted designation applied to a separate linker/warhead-remapping workflow.

## Figure 1 and Figure 2

`figure1/final_counts_used.csv` contains model-valency summaries and one row per retained construct per definition. Proximity means any AD heavy atom within 5 Å of a Tail residue. Strict hydrophobic contact requires a side-chain atom distance ≤4.5 Å and hydrophobic identities W/F/Y/L/I/V/M on both sides. Multisubunit means at least two physical Tail subunits; recurrent means this occurs in at least two retained models of one AD.

`figure2/contact_profile_source_data.csv` has one row per series and offset (−30…30 residues). `mean_contact_percentile` is the upstream per-residue Tail-contact score converted to a within-AD percentile before deletion-centered aggregation. Stored CIs are upstream gene-cluster bootstrap intervals; `n_valid_at_offset` varies because boundary positions are unavailable for some windows. The renderer applies only a five-position centered display mean. Panel C activation-impairing is LoF OR Necessary (n=1,643). Comparison to the stored deletion-coordinate table found no deletion interval crossing a synthetic junction; some ±30-residue profile windows do span one, because offsets are defined on the modeled construct.

`figure2/wfyl_source_data.csv` has one row per deletion class. `wfyl_free_count` counts deleted regions containing zero W, F, Y, or L; `fraction` and `percent` use `n_windows` as denominator. CIs are stored two-sided gene-cluster bootstrap intervals. Panel D's activation-impairing bar uses LoF only (n=1,341), intentionally differing from panel C. `GoF Q85` is the top-15% activation-enhancing subset (`Delta_minCMV ≥ 2.260200281716645`); `GoF +0.5` uses `Delta_minCMV ≥ 1.6205657816946992`. Historical `robust_z` fields are retained only where present; this deposit does not establish their unavailable upstream scaling constant.

## Extended Data

`extended_data_figure1/interpeak_interval_source_data.csv` has one row per displayed series/interval occurrence. Coordinates and lengths are along the modeled construct and can cross a synthetic junction; they are not necessarily uninterrupted native-protein distances. Nine unique all-interpeak intervals cross such a junction; three of their repeated nested-series occurrences are also present (12 rows total). Series are nested. Identifiers come from the stored upstream interval table rather than length matching. The 196-residue reference interval is retained but lies beyond the plot's 145-residue display limit.

`extended_data_figure1/q85_source_manifest.json` stores the exact activation thresholds.

The five `extended_data_figure2/*.csv` files contain contact-decile estimates and deletion-class estimates/contrasts. Class tables use one row per series/metric; contrast tables use one row per series/reference/metric. Fractions are unitless; CIs are stored upstream bootstrap intervals. `q_BH_across_all_28_deletion_contrasts` is the pLDDT-only 28-test correction. The renderer and manuscript annotations instead recompute BH values from the complete supplied family of 28 pLDDT plus 8 helicity P values (36 tests); no raw P values are discarded.

`extended_data_table1/extended_data_table_1.csv` is one row per published reconstruction. `EMD-31215*` refers to the deposited WT PIC-Mediator map; the listed populations/resolutions are its reported pre/intermediate/holo reconstructions. The 3.6 Å focused WT Tail reconstruction (44,667 particles) preceded state-specific separation and was not assigned to MEDE or MEDB. `n.r.` means not reported separately by state; EMDB and PDB are archive accession systems; ΔMED1-IDR is deletion of the MED1 intrinsically disordered region.

## Figure 3 and RMSD summaries

`figure3/figure3b_nodes.csv` has one row per displayed receptor node. `strict_contact_mass` is the sum of per-pair weights `(1 - (d_CA/12 Å)^2)^2` for strict hydrophobic pairs with Cα distance below 12 Å after masks; it is a weighted contact score, not physical mass, a raw pair count, or the separately normalized ensemble surface map. `strict_pair_count`, `supporting_models`, and `supporting_ADs` are unweighted counts. Node area uses the documented affine display scale.

`figure3/figure3b_edges.csv` contains all 27 nonzero unordered site pairs. `supporting_models` counts retained models with strict engagement of both sites; fractions use 1,139. `figure3b_positions.csv` stores the schematic two-dimensional layout. Positions, edge paths, and curves are not structural distances.

`figure3/figure3c_reported_geometry.csv` is a frozen numerical report. MED23 used 1,334 sequence-matched Cα pairs and a global no-rejection Kabsch fit. MED24 used 897 PyMOL-matched Cα pairs ordered by Extended-object chain/residue, scanned two contiguous blocks with at least 60 pairs each, and minimized their combined squared residual. `X226/X238` is the best boundary between adjacent matched atoms in the Extended/Bent PyMOL objects; `0.00;1.48` are the two block-fit RMSDs at printed precision. The scripts/session were identified, but this deposit does not contain a self-contained coordinate-source manifest or transforms.

`rmsd/af3_vs_8gxs_per_job.csv` contains one row per construct. All 293 selected model 0, which has the maximum recorded ranking score (ties included), and use 3,328 paired Cα atoms. `rmsd/af3_vs_8gxs_summary.json` records the chain selection, global Kabsch fit, no-outlier-rejection policy, median, and IQR. The exact semiglobal sequence-matching implementation and coordinates are not deposited, so atom-level recalculation is outside this repository.
