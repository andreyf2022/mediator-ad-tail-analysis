#!/usr/bin/env python3
"""Validate cross-table identities in the coordinate-free release.

Run from any working directory with ``python tests/validate_release.py``.
This checks deposited values and joins; it does not recalculate contacts,
bootstraps, structural alignments, or AF3 predictions.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def read_csv(relative: str, **kwargs) -> pd.DataFrame:
    path = ROOT / relative
    assert path.is_file(), f"missing {relative}"
    return pd.read_csv(path, **kwargs)


def validate_constructs() -> None:
    manifest = read_csv("data/dataset/construct_model_manifest.csv", dtype=str).fillna("")
    assert len(manifest) == 293
    assert manifest["construct_index"].astype(int).tolist() == list(range(2, 295))
    assert manifest["construct_index"].is_unique and manifest["gene"].is_unique and manifest["fold_dir"].is_unique
    assert (manifest["modeled_ad_sequence"].str.len() == manifest["ad_length"].astype(int)).all()
    assert set(manifest["model_ids"]) == {"0;1;2;3;4"}
    assert set(manifest["af3_platform"]) == {"AlphaFold Server"}
    assert set(manifest["af3_version_recorded"]) == {"AlphaFold-beta-20231127"}
    assert manifest["af3_version_confirmation_status"].str.contains("mmCIF metadata", regex=False).all()
    assert int(manifest["n_retained_contact_models"].astype(int).sum()) == 1139
    assert int(manifest["in_234_contact_ensemble"].eq("yes").sum()) == 234

    statuses = pd.concat([manifest[f"model_{i}_contact_status"] for i in range(5)], ignore_index=True)
    assert len(statuses) == 1465
    assert int(statuses.str.startswith("INITIAL_EXCLUDED:").sum()) == 86
    assert int(statuses.str.startswith("POSTMASK_EXCLUDED:").sum()) == 240
    assert int(statuses.eq("RETAINED").sum()) == 1139
    initial = statuses[statuses.str.startswith("INITIAL_EXCLUDED:")]
    assert initial.str.fullmatch(
        r"INITIAL_EXCLUDED:heavy_atom_overlap_screen;n_AD_heavy_atoms_within_1p2A=\d+"
    ).all()

    components = read_csv("data/dataset/construct_source_components.csv", dtype=str).fillna("")
    assert len(components) == 374
    components["construct_index"] = components["construct_index"].astype(int)
    components["component_order"] = components["component_order"].astype(int)
    assert components[["construct_index", "component_order"]].duplicated().sum() == 0
    reconstructed = (
        components.sort_values(["construct_index", "component_order"])
        .groupby("construct_index", sort=True)["appended_sequence"].sum()
    )
    expected = manifest.set_index(manifest["construct_index"].astype(int))["modeled_ad_sequence"].sort_index()
    assert reconstructed.index.astype(int).tolist() == expected.index.tolist()
    assert reconstructed.tolist() == expected.tolist()
    noncontiguous = components.loc[components["synthetic_junction_after_previous_component"].eq("yes")]
    assert noncontiguous["construct_index"].nunique() == 22

    model_status = {}
    for row in manifest.itertuples(index=False):
        for model_number in range(5):
            model_status[f"{row.fold_dir}__model_{model_number}"] = getattr(row, f"model_{model_number}_contact_status")
    inpaint = read_csv("data/dataset/inpainting_operational_ledger.csv", dtype=str).fillna("")
    masked = 0
    partial = inpaint.loc[inpaint["final_operational_class"].eq("PARTIAL_CONTACT_MASK")]
    for row in partial.itertuples(index=False):
        note = row.masked_region_note
        match = re.search(r"chainH:\s*(\d+)-(\d+)", note)
        assert match, note
        start, end = map(int, match.groups())
        assert start >= 1 and end >= start
        if model_status[row.model_id] == "RETAINED":
            masked += end - start + 1
    assert masked == 913
    threading = read_csv("data/dataset/threading_local_mask_residues.csv")
    assert len(threading) == 61
    assert threading[["model_id", "ad_resnum"]].duplicated().sum() == 0
    assert threading["in_threading_local_mask"].astype(bool).all()


def validate_figure1_2() -> None:
    counts = read_csv("data/figure1/final_counts_used.csv")
    per_ad = counts.loc[counts["row_type"].eq("ad_level_fraction")]
    expected = {"proximity": (986, 219, 1.0), "strict_hydrophobic": (819, 195, 0.8)}
    for definition, (multi_total, recurrent, median) in expected.items():
        frame = per_ad.loc[per_ad["definition"].eq(definition)]
        assert len(frame) == 234 and not frame["fold_dir"].duplicated().any()
        usable = frame["n_models_usable"].astype(int)
        multi = frame["n_models_multisubunit"].astype(int)
        fraction = frame["fraction_models_multisubunit"].astype(float)
        assert usable.sum() == 1139 and multi.sum() == multi_total
        assert (multi >= 2).sum() == recurrent and float(fraction.median()) == median
        assert np.allclose(fraction, multi / usable, rtol=0, atol=1e-15)

    contact = read_csv("data/figure2/contact_profile_source_data.csv")
    assert not contact.duplicated(["series", "offset"]).any()
    for _, frame in contact.groupby("series"):
        assert sorted(frame["offset"].astype(int)) == list(range(-30, 31))
        assert frame["n_windows"].nunique() == 1
    expected_contact = {"GoF Q85": 38, "Neutral": 1838, "LoF + Necessary": 1643}
    assert {key: int(contact.loc[contact.series.eq(key), "n_windows"].iloc[0]) for key in expected_contact} == expected_contact
    assert ((contact.ci95_low <= contact.mean_contact_percentile) & (contact.mean_contact_percentile <= contact.ci95_high)).all()
    assert ((contact.n_valid_at_offset >= 0) & (contact.n_valid_at_offset <= contact.n_windows)).all()

    wfyl = read_csv("data/figure2/wfyl_source_data.csv")
    expected_wfyl = {"Neutral": 1838, "GoF": 234, "GoF +0.5": 101, "GoF Q85": 38, "LoF": 1341}
    assert dict(zip(wfyl.series, wfyl.n_windows.astype(int))) == expected_wfyl
    assert np.allclose(wfyl.fraction, wfyl.wfyl_free_count / wfyl.n_windows, rtol=0, atol=1e-12)
    assert np.allclose(wfyl.percent, 100 * wfyl.fraction, rtol=0, atol=1e-10)


def validate_extended_data() -> None:
    intervals = read_csv("data/extended_data_figure1/interpeak_interval_source_data.csv")
    assert not intervals.duplicated(["series", "interval_id"]).any()
    expected = {
        "All interpeak": (283, 28.0, 15.5, 42.0, 196.0),
        "GoF-targeted": (39, 43.0, 32.0, 57.5, 134.0),
        "GoF +0.5-targeted": (14, 47.0, 38.5, 62.75, 134.0),
        "GoF Q85-targeted": (5, 43.0, 40.0, 51.0, 52.0),
    }
    for series, values in expected.items():
        x = intervals.loc[intervals.series.eq(series), "interpeak_interval_length"].to_numpy(float)
        observed = (len(x), float(np.median(x)), *map(float, np.percentile(x, [25, 75])), float(x.max()))
        assert observed == values
    assert np.array_equal(
        (intervals.construct_end - intervals.construct_start + 1).astype(int),
        intervals.interpeak_interval_length.astype(int),
    )

    deciles = read_csv("data/extended_data_figure2/contact_decile_plddt_helicity.csv")
    assert deciles.contact_decile.astype(int).tolist() == list(range(1, 11))
    assert ((deciles.n_ADs > 0) & (deciles.n_ADs <= 234)).all()
    plddt = read_csv("data/extended_data_figure2/deletion_plddt_contrasts.csv")
    helix = read_csv("data/extended_data_figure2/deletion_helicity_contrasts.csv")
    assert len(plddt) == 28 and len(helix) == 8
    assert not plddt.duplicated(["series", "reference", "metric"]).any()
    assert not helix.duplicated(["series", "reference", "metric"]).any()
    pvalues = np.concatenate([plddt.bootstrap_two_sided_directional_p, helix.two_sided_bootstrap_by_gene_p])
    assert len(pvalues) == 36 and np.isfinite(pvalues).all() and ((pvalues >= 0) & (pvalues <= 1)).all()


def validate_figure3_rmsd() -> None:
    sites = read_csv("data/dataset/receptor_site_definitions.csv", dtype=str).fillna("")
    assert len(sites) == 9 and sites.receptor_site_id.is_unique
    assert sites.included_in_figure3b_network.eq("yes").sum() == 8
    d653 = sites.loc[sites.receptor_site_id.eq("SITE_MED24_22_HAZE")].iloc[0]
    assert d653.included_in_figure3b_network == "yes" and "without downweighting" in d653.methodological_caveat
    tunnel = sites.loc[sites.receptor_site_id.eq("SITE_MED16_MED25_TUNNEL_COMPOSITE")].iloc[0]
    assert ">=0.65" in tunnel.shoulder_support_rule and ">=0.80" in tunnel.shoulder_support_rule

    nodes = read_csv("data/figure3/figure3b_nodes.csv")
    edges = read_csv("data/figure3/figure3b_edges.csv")
    positions = read_csv("data/figure3/figure3b_positions.csv")
    assert len(nodes) == len(positions) == 8 and len(edges) == 27
    assert set(nodes.hotspot_key) == set(positions.hotspot_key)
    assert nodes.hotspot_key.is_unique and positions.hotspot_key.is_unique
    assert nodes.loc[nodes.hotspot_key.eq("MED15_L76"), "residue_label"].item() == "L595"
    pairs = [tuple(sorted(pair)) for pair in zip(edges.hotspot_1, edges.hotspot_2)]
    assert len(set(pairs)) == 27 and all(a != b for a, b in pairs)
    assert np.allclose(edges.fraction_usable_models, edges.supporting_models / 1139, rtol=0, atol=1e-15)
    support = nodes.set_index("hotspot_key").supporting_models.astype(int)
    assert all(row.supporting_models <= min(support[row.hotspot_1], support[row.hotspot_2]) for row in edges.itertuples())

    per_job = read_csv("data/rmsd/af3_vs_8gxs_per_job.csv")
    summary = json.loads((ROOT / "data/rmsd/af3_vs_8gxs_summary.json").read_text())
    manifest = read_csv("data/dataset/construct_model_manifest.csv")
    merged = per_job.merge(manifest, on="fold_dir", validate="one_to_one")
    assert len(merged) == 293 and per_job.model_id.eq(0).all() and per_job.paired_CA_atoms.eq(3328).all()
    score_columns = [f"ranking_score_model_{i}" for i in range(5)]
    assert np.allclose(merged.ranking_score, merged.ranking_score_model_0, rtol=0, atol=1e-15)
    assert (merged.ranking_score_model_0 >= merged[score_columns].max(axis=1)).all()
    rmsd = per_job.global_Kabsch_RMSD_A.to_numpy(float)
    assert np.isclose(np.median(rmsd), summary["median_RMSD_A"], rtol=0, atol=1e-15)
    assert np.allclose(np.percentile(rmsd, [25, 75]), summary["IQR_RMSD_A"], rtol=0, atol=1e-15)
    geometry = read_csv("data/figure3/figure3c_reported_geometry.csv")
    assert geometry.paired_CA_atoms.astype(int).tolist() == [1334, 897]
    assert np.allclose(geometry.rotation_degrees, [46.13, 52.14], rtol=0, atol=1e-12)


def validate_metadata_and_manifest() -> None:
    zenodo = json.loads((ROOT / ".zenodo.json").read_text())
    assert zenodo["version"] == "1.0.1" and zenodo["license"] == "MIT"
    cff = (ROOT / "CITATION.cff").read_text()
    assert "version: 1.0.1" in cff and "license: MIT" in cff
    manifest = read_csv("release_manifest.csv", dtype=str).fillna("")
    assert manifest.staged_path.is_unique
    listed = set(manifest.staged_path)
    actual = {
        str(path.relative_to(ROOT))
        for path in ROOT.rglob("*")
        if path.is_file()
        and ".git" not in path.parts
        and "__pycache__" not in path.parts
        and not str(path.relative_to(ROOT)).startswith("outputs/")
    }
    assert listed == actual - {"release_manifest.csv"}
    for row in manifest.itertuples(index=False):
        digest = hashlib.sha256((ROOT / row.staged_path).read_bytes()).hexdigest()
        assert digest == row.sha256, row.staged_path
    for path in actual:
        if path in {"release_manifest.csv", "tests/validate_release.py"}:
            continue
        raw = (ROOT / path).read_bytes()
        if b"\x00" not in raw:
            text = raw.decode("utf-8")
            assert "/Users/" not in text and "andrey1" not in text, path


def main() -> int:
    validate_constructs()
    validate_figure1_2()
    validate_extended_data()
    validate_figure3_rmsd()
    validate_metadata_and_manifest()
    print("PASS: all coordinate-free release consistency checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
