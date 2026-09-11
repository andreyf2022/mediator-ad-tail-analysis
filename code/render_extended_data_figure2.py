#!/usr/bin/env python3
"""Render the pLDDT/helicity summary in Extended Data Figure 2.

Purpose: redraw frozen class/contrast summaries and calculate only the displayed
Benjamini-Hochberg adjustment across the supplied 36 contrast P values.
Manuscript support: Extended Data Figure 2A/B.
Required inputs: five CSVs under ``data/extended_data_figure2``; each path can
be overridden with the corresponding command-line argument.
Outputs: PNG, PDF and SVG versions in a new output directory.
Example: ``python code/render_extended_data_figure2.py --output-dir outputs/extended_data_figure2``
Fixed assumptions: contact and pLDDT ranks are within AD; Panel A confidence
intervals are AD-bootstrap intervals; Panel B intervals are gene-bootstrap
intervals; top pLDDT means the within-AD top quintile; activation-impairing
means LoF OR necessary. The frozen upstream consensus-helix rule is STRIDE H in
>=40% of usable models, retained in runs of >=3 residues. This script does not
rerun bootstraps, STRIDE or residue mapping and never overwrites an existing
output directory.
"""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib_plddt_helicity"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "xdg_plddt_helicity"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/extended_data_figure2"
DEFAULT_CONTACT_DECILES = DATA / "contact_decile_plddt_helicity.csv"
DEFAULT_HELIX_CLASSES = DATA / "deletion_helicity_class_summary.csv"
DEFAULT_HELIX_CONTRASTS = DATA / "deletion_helicity_contrasts.csv"
DEFAULT_PLDDT_CLASSES = DATA / "deletion_plddt_class_summary.csv"
DEFAULT_PLDDT_CONTRASTS = DATA / "deletion_plddt_contrasts.csv"

BLUE = "#315f8a"
GREEN = "#08733f"
GRAY = "#8a939b"
ORANGE = "#f18e2c"
GRID = "#d8dde2"
ANNOTATION = "#596169"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contact-deciles-csv", type=Path, default=DEFAULT_CONTACT_DECILES)
    parser.add_argument("--helix-classes-csv", type=Path, default=DEFAULT_HELIX_CLASSES)
    parser.add_argument("--helix-contrasts-csv", type=Path, default=DEFAULT_HELIX_CONTRASTS)
    parser.add_argument("--plddt-classes-csv", type=Path, default=DEFAULT_PLDDT_CLASSES)
    parser.add_argument("--plddt-contrasts-csv", type=Path, default=DEFAULT_PLDDT_CONTRASTS)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.5,
            "axes.titlesize": 10.0,
            "axes.labelsize": 9.2,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def save_figure(fig: plt.Figure, stem: str, output_dir: Path) -> None:
    for suffix, kwargs in [(".png", {"dpi": 600}), (".pdf", {}), (".svg", {})]:
        fig.savefig(output_dir / f"{stem}{suffix}", bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(fig)


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def draw_panel_a(ax: plt.Axes, deciles: pd.DataFrame) -> None:
    x = deciles["contact_decile"].to_numpy(float)
    series = [
        ("top_plddt_fraction", "Top pLDDT quintile", BLUE, "o"),
        ("helix_fraction", "Consensus helix", GREEN, "s"),
    ]
    for column, label, color, marker in series:
        y = 100 * deciles[column].to_numpy(float)
        low = 100 * deciles[f"{column}_AD_bootstrap_ci95_low"].to_numpy(float)
        high = 100 * deciles[f"{column}_AD_bootstrap_ci95_high"].to_numpy(float)
        ax.fill_between(x, low, high, color=color, alpha=0.10, linewidth=0)
        ax.plot(
            x,
            y,
            color=color,
            marker=marker,
            markersize=4.1,
            linewidth=1.85,
            label=label,
            zorder=3,
        )
    ax.set_xlim(0.65, 10.35)
    ax.set_xticks(np.arange(1, 11))
    ax.set_xlabel("Tail-contact traffic decile within AD")
    ax.set_ylabel("Fraction of residues (%)")
    ax.set_title("A  Local confidence and helicity across Tail-contact traffic", loc="left", weight="bold", pad=7)
    ax.grid(axis="y", color=GRID, linewidth=0.45, alpha=0.42)
    ax.legend(frameon=False, loc="upper left", fontsize=7.7, handlelength=2.3, labelspacing=0.55)


def draw_panel_b(
    ax: plt.Axes,
    helix_classes: pd.DataFrame,
    helix_contrasts: pd.DataFrame,
    plddt_classes: pd.DataFrame,
    plddt_contrasts: pd.DataFrame,
) -> None:
    plddt_metric = "deletion_interval_fraction_top_relative_plddt_quintile"
    helix_metric = "deletion_interval_fraction_consensus_STRIDE_H"
    groups = [
        (
            "Top pLDDT\nquintile",
            plddt_classes.loc[
                plddt_classes["metric"].eq(plddt_metric)
                & plddt_classes["series"].isin(["Neutral", "GoF Q85", "LoF + Necessary"])
            ],
            "gene_cluster_bootstrap_ci95_low",
            "gene_cluster_bootstrap_ci95_high",
        ),
        (
            "Consensus\nhelix",
            helix_classes.loc[
                helix_classes["metric"].eq(helix_metric)
                & helix_classes["series"].isin(["Neutral", "GoF Q85", "LoF + Necessary"])
            ],
            "bootstrap_by_gene_ci95_low",
            "bootstrap_by_gene_ci95_high",
        ),
    ]
    displayed_series = ["Neutral", "LoF + Necessary"]
    offsets = {"Neutral": -0.11, "LoF + Necessary": 0.11}
    colors = {"Neutral": GRAY, "LoF + Necessary": ORANGE}
    for x_center, (_, frame, low_column, high_column) in enumerate(groups):
        for series in displayed_series:
            row = frame.loc[frame["series"].eq(series)].iloc[0]
            value = 100 * float(row["estimate"])
            low = 100 * float(row[low_column])
            high = 100 * float(row[high_column])
            ax.errorbar(
                x_center + offsets[series],
                value,
                yerr=[[value - low], [high - value]],
                fmt="o",
                markersize=6.8,
                color=colors[series],
                ecolor=colors[series],
                elinewidth=1.55,
                capsize=3.0,
                zorder=3,
            )

    plddt_p = plddt_contrasts["bootstrap_two_sided_directional_p"].to_numpy(float)
    helix_p = helix_contrasts["two_sided_bootstrap_by_gene_p"].to_numpy(float)
    combined_q = bh_adjust(np.concatenate([plddt_p, helix_p]))
    plddt_target_index = plddt_contrasts.index[
        plddt_contrasts["series"].eq("LoF + Necessary")
        & plddt_contrasts["metric"].eq(plddt_metric)
    ][0]
    helix_target_index = helix_contrasts.index[
        helix_contrasts["series"].eq("LoF + Necessary")
        & helix_contrasts["metric"].eq(helix_metric)
    ][0]
    plddt_adjusted_p = combined_q[int(plddt_target_index)]
    helix_adjusted_p = combined_q[len(plddt_contrasts) + int(helix_target_index)]
    labels = [
        f"BH-adjusted $P$ = {plddt_adjusted_p:.4f}",
        f"BH-adjusted $P$ = {helix_adjusted_p:.4f}",
    ]
    for x_center, label, y in zip([0, 1], labels, [34.0, 44.0]):
        left = x_center + offsets["Neutral"]
        right = x_center + offsets["LoF + Necessary"]
        ax.plot([left, left, right, right], [y - 0.7, y, y, y - 0.7], color=ANNOTATION, linewidth=0.8)
        ax.text(x_center, y + 0.6, label, ha="center", va="bottom", fontsize=7.0, color=ANNOTATION)

    ax.set_xlim(-0.48, 1.48)
    ax.set_xticks([0, 1], [groups[0][0], groups[1][0]])
    ax.set_xlabel("")
    ax.set_ylabel("Fraction of deleted residues (%)")
    ax.set_title("B  Local confidence and helicity in deleted residues", loc="left", weight="bold", pad=7)
    ax.grid(axis="y", color=GRID, linewidth=0.45, alpha=0.42)
    ax.scatter([], [], color=GRAY, s=38, label="Neutral")
    ax.scatter([], [], color=ORANGE, s=38, label="Activation-impairing")
    ax.legend(frameon=False, loc="upper left", fontsize=7.4, handletextpad=0.45, labelspacing=0.42)


def main() -> int:
    args = parse_args()
    contact_deciles_path = args.contact_deciles_csv.expanduser().resolve()
    helix_classes_path = args.helix_classes_csv.expanduser().resolve()
    helix_contrasts_path = args.helix_contrasts_csv.expanduser().resolve()
    plddt_classes_path = args.plddt_classes_csv.expanduser().resolve()
    plddt_contrasts_path = args.plddt_contrasts_csv.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    required = [
        contact_deciles_path,
        helix_classes_path,
        helix_contrasts_path,
        plddt_classes_path,
        plddt_contrasts_path,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(missing)
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output directory: {output_dir}")
    output_dir.mkdir(parents=True)
    setup_style()
    deciles = pd.read_csv(contact_deciles_path)
    helix_classes = pd.read_csv(helix_classes_path)
    helix_contrasts = pd.read_csv(helix_contrasts_path)
    plddt_classes = pd.read_csv(plddt_classes_path)
    plddt_contrasts = pd.read_csv(plddt_contrasts_path)
    if deciles["contact_decile"].astype(int).tolist() != list(range(1, 11)):
        raise ValueError("Panel A must contain contact deciles 1 through 10 in order")
    if len(plddt_contrasts) != 28 or len(helix_contrasts) != 8:
        raise ValueError("Expected 28 pLDDT and 8 helicity contrasts for BH adjustment")

    fig, axes = plt.subplots(1, 2, figsize=(8.35, 3.15), gridspec_kw={"width_ratios": [1.18, 0.82]})
    draw_panel_a(axes[0], deciles)
    draw_panel_b(axes[1], helix_classes, helix_contrasts, plddt_classes, plddt_contrasts)
    axes[0].set_ylim(0, 42)
    axes[0].set_yticks(np.arange(0, 41, 10))
    axes[1].set_ylim(0, 50)
    axes[1].set_yticks(np.arange(0, 51, 10))
    for ax in axes:
        ax.tick_params(width=0.8, length=3.5)
        ax.spines["left"].set_linewidth(0.8)
        ax.spines["bottom"].set_linewidth(0.8)
    fig.text(
        0.5,
        0.012,
        "pLDDT and contact-traffic rankings are within each AD. A, 95% CIs from bootstrap resampling by AD; B, 95% CIs from bootstrap resampling by gene.",
        ha="center",
        fontsize=7.1,
        color=ANNOTATION,
    )
    fig.subplots_adjust(left=0.075, right=0.995, top=0.90, bottom=0.27, wspace=0.34)
    save_figure(fig, "plddt_helicity_extended_data_two_panel_20260828", output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
