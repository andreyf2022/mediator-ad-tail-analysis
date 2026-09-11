#!/usr/bin/env python3
"""Render the Figure 1B/C multisubunit-engagement quantitative module.

Purpose: validate and redraw the frozen final model-valency and recurrent-AD
counts; this script does not call contacts or reconstruct receptor occupancy.
Manuscript support: Figure 1B/C and the corresponding 1,139-model/234-AD
multisubunit statistics.
Required input: ``data/figure1/final_counts_used.csv`` (or ``--counts-csv``).
Output: one PNG supplied with ``--output``.
Example: ``python code/render_figure1_multisubunit.py --output outputs/figure1_multisubunit.png``
Fixed assumptions: physical-subunit valency bins are 0, 1, 2 and >=3;
multisubunit means >=2 occupied Tail subunits; recurrent means support in >=2
usable models of the same AD. Expected denominators are 1,139 models and 234
AD constructs. Existing output files are never overwritten.
"""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib_figure1_multisubunit")
)
os.environ.setdefault(
    "XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "xdg_figure1_multisubunit")
)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_COUNTS = PACKAGE_ROOT / "data/figure1/final_counts_used.csv"
EXPECTED_MODEL_COUNTS = {
    "proximity": [15, 138, 161, 825],
    "strict_hydrophobic": [30, 290, 367, 452],
}
EXPECTED_MULTI_COUNTS = {"proximity": 986, "strict_hydrophobic": 819}
EXPECTED_RECURRENT_COUNTS = {"proximity": 219, "strict_hydrophobic": 195}
EXPECTED_MEDIANS = {"proximity": 1.0, "strict_hydrophobic": 0.8}
EXPECTED_DENOMINATORS = {2: 1, 3: 4, 4: 20, 5: 209}
EXPECTED_N_MODELS = 1139
EXPECTED_N_ADS = 234


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-csv", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def require_columns(frame: pd.DataFrame, columns: set[str]) -> None:
    missing = sorted(columns - set(frame.columns))
    if missing:
        raise ValueError(f"Final counts table is missing columns: {missing}")


def load_plot_authority(counts_path: Path) -> dict[str, object]:
    counts = pd.read_csv(counts_path, low_memory=False)
    require_columns(
        counts,
        {
            "row_type",
            "definition",
            "category",
            "fold_dir",
            "fraction_models_multisubunit",
            "count",
            "fraction",
        },
    )
    categories = ["0", "1", "2", ">=3"]
    model_counts: dict[str, list[int]] = {}
    multi_counts: dict[str, int] = {}
    recurrent_counts: dict[str, int] = {}
    medians: dict[str, float] = {}
    for definition in ("proximity", "strict_hydrophobic"):
        summary = counts.loc[
            counts["row_type"].eq("model_valency_summary")
            & counts["definition"].eq(definition)
        ].copy()
        summary["category"] = pd.Categorical(
            summary["category"], categories=categories, ordered=True
        )
        summary = summary.sort_values("category", kind="mergesort")
        model_counts[definition] = summary["count"].astype(int).tolist()
        multi = counts.loc[
            counts["row_type"].eq("model_multisubunit_summary")
            & counts["definition"].eq(definition),
            "count",
        ]
        recurrent = counts.loc[
            counts["row_type"].eq("ad_recurrent_ge2_models_summary")
            & counts["definition"].eq(definition),
            "count",
        ]
        median = counts.loc[
            counts["row_type"].eq("ad_fraction_median_summary")
            & counts["definition"].eq(definition),
            "fraction",
        ]
        if len(multi) != 1 or len(recurrent) != 1 or len(median) != 1:
            raise ValueError(f"Final counts table is incomplete for {definition}")
        multi_counts[definition] = int(multi.iloc[0])
        recurrent_counts[definition] = int(recurrent.iloc[0])
        medians[definition] = float(median.iloc[0])

    denominator_rows = counts.loc[
        counts["row_type"].eq("ad_usable_model_denominator_summary")
    ]
    denominators = {
        int(row.category): int(row.count)
        for row in denominator_rows.itertuples(index=False)
    }
    per_ad = counts.loc[counts["row_type"].eq("ad_level_fraction")]
    if len(per_ad) != 2 * EXPECTED_N_ADS:
        raise ValueError(f"Expected {2 * EXPECTED_N_ADS} per-AD rows; found {len(per_ad)}")
    if per_ad.duplicated(["definition", "fold_dir"]).any():
        raise ValueError("Per-AD fractions are not unique by definition and fold")

    observed = {
        "model_valency_counts": model_counts,
        "model_multisubunit_counts": multi_counts,
        "ad_recurrent_ge2_counts": recurrent_counts,
        "ad_fraction_medians": medians,
        "usable_model_denominator_distribution": denominators,
    }
    expected = {
        "model_valency_counts": EXPECTED_MODEL_COUNTS,
        "model_multisubunit_counts": EXPECTED_MULTI_COUNTS,
        "ad_recurrent_ge2_counts": EXPECTED_RECURRENT_COUNTS,
        "ad_fraction_medians": EXPECTED_MEDIANS,
        "usable_model_denominator_distribution": EXPECTED_DENOMINATORS,
    }
    if observed != expected:
        raise ValueError(f"Final Figure 1B/C values changed: {observed}")
    return observed


def draw_figure(values: dict[str, object], output_path: Path) -> None:
    rc = {
        "font.family": "DejaVu Sans",
        "font.size": 7.5,
        "axes.titlesize": 8.5,
        "axes.labelsize": 7.6,
        "xtick.labelsize": 6.9,
        "ytick.labelsize": 7.2,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
    }
    category_colors = ["#EFEFEF", "#C8C8C8", "#8A8A8A", "#151515"]
    recurrence_colors = ["#8A8A8A", "#303030"]
    categories = ["0", "1", "2", "≥3"]
    definitions = ["proximity", "strict_hydrophobic"]
    row_labels = ["Proximity", "Strict hydrophobic"]
    y_positions = [1, 0]

    with mpl.rc_context(rc):
        fig = plt.figure(figsize=(10.6, 3.20), facecolor="white")
        ax_a = fig.add_axes([0.10, 0.34, 0.62, 0.56])
        for definition, y_position in zip(definitions, y_positions):
            counts = np.asarray(values["model_valency_counts"][definition], dtype=float)
            fractions = counts / EXPECTED_N_MODELS
            left = 0.0
            for category, fraction, color in zip(categories, fractions, category_colors):
                ax_a.barh(
                    y_position,
                    fraction,
                    left=left,
                    height=0.42,
                    color=color,
                    edgecolor="white",
                    linewidth=0.55,
                    zorder=2,
                )
                if fraction >= 0.08:
                    ax_a.text(
                        left + fraction / 2,
                        y_position,
                        f"{fraction:.1%}",
                        ha="center",
                        va="center",
                        fontsize=6.7,
                        color="white" if category == "≥3" else "#22282C",
                        fontweight="bold" if category == "≥3" else "normal",
                    )
                left += fraction

        ax_a.set_title("Occupancy in individual AF3 models", loc="left", fontweight="normal", color="black", pad=7)
        ax_a.text(1.0, 1.025, "n = 1,139 usable AF3 models", transform=ax_a.transAxes, ha="right", va="bottom", fontsize=6.2, color="black")
        ax_a.text(1.025, 1, "≥2: 86.6% (986/1,139)", ha="left", va="center", fontsize=6.7)
        ax_a.text(1.025, 0, "≥2: 71.9% (819/1,139)", ha="left", va="center", fontsize=6.7)
        ax_a.set_yticks(y_positions, row_labels)
        ax_a.set_xlim(0, 1.30)
        ax_a.set_ylim(-0.58, 1.62)
        ax_a.set_xticks([0, 0.25, 0.50, 0.75, 1.0], ["0", "25", "50", "75", "100"])
        ax_a.set_xlabel("")
        ax_a.tick_params(axis="y", length=0)
        ax_a.tick_params(axis="both", colors="black")
        ax_a.spines[["top", "right", "left"]].set_visible(False)
        ax_a.spines["bottom"].set_bounds(0, 1.03)
        legend = ax_a.legend(
            handles=[Patch(facecolor=color, edgecolor="none", label=category) for category, color in zip(categories, category_colors)],
            title="Occupied subunits",
            title_fontsize=6.1,
            frameon=False,
            ncol=4,
            loc="upper left",
            bbox_to_anchor=(0.0, -0.34),
            fontsize=6.1,
            handlelength=0.75,
            handletextpad=0.28,
            columnspacing=0.48,
            borderaxespad=0.0,
            labelspacing=0.20,
        )
        legend.get_title().set_color("black")
        legend.get_title().set_fontweight("normal")
        for legend_text in legend.get_texts():
            legend_text.set_color("black")
            legend_text.set_fontweight("normal")

        ax_b = fig.add_axes([0.78, 0.42, 0.18, 0.32])
        heights = [values["ad_recurrent_ge2_counts"][definition] / EXPECTED_N_ADS for definition in definitions]
        x_positions = np.arange(2)
        ax_b.bar(x_positions, heights, width=0.50, color=recurrence_colors, edgecolor="none", zorder=2)
        for x_position, height, definition in zip(x_positions, heights, definitions):
            count = values["ad_recurrent_ge2_counts"][definition]
            ax_b.text(x_position, height + 0.035, f"{height:.1%}", ha="center", va="bottom", fontsize=6.7, color="black")
            ax_b.text(x_position, -0.055, f"{count}/234", ha="center", va="top", fontsize=6.1, color="black")
        ax_b.text(0.0, 1.52, "Recurrence among models of the same AD", transform=ax_b.transAxes, ha="left", va="bottom", fontsize=8.5, color="black")
        ax_b.text(0.0, 1.25, "≥2 occupied subunits in ≥2 models", transform=ax_b.transAxes, ha="left", va="bottom", fontsize=6.2, color="black")
        ax_b.set_xlim(-0.65, 1.65)
        ax_b.set_ylim(-0.16, 1.08)
        ax_b.set_xticks(x_positions, ["Proximity", "Strict\nhydrophobic"])
        ax_b.tick_params(axis="x", length=0, pad=13, labelsize=6.1)
        for tick_label in ax_b.get_xticklabels():
            tick_label.set_color("black")
            tick_label.set_fontweight("normal")
        ax_b.set_yticks([])
        ax_b.axhline(0, color="#B8C0C4", linewidth=0.5, zorder=1)
        ax_b.spines[["top", "right", "left", "bottom"]].set_visible(False)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=500, facecolor="white", bbox_inches="tight", pad_inches=0.08)
        plt.close(fig)


def main() -> int:
    args = parse_args()
    counts_path = args.counts_csv.expanduser().resolve()
    output_path = args.output.expanduser().resolve()
    if not counts_path.is_file():
        raise FileNotFoundError(f"Missing final counts table: {counts_path}")
    if output_path.suffix.lower() != ".png":
        raise ValueError("--output must end in .png")
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {output_path}")
    values = load_plot_authority(counts_path)
    draw_figure(values, output_path)
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
