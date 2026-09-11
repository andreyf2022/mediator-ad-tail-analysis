#!/usr/bin/env python3
"""Summarize and render Extended Data Figure 1 interpeak lengths.

Purpose: calculate descriptive n/median/IQR values from the frozen interval
table and draw the manuscript-facing distributions.
Manuscript support: Extended Data Figure 1, interpeak interval lengths.
Required inputs: ``data/extended_data_figure1/interpeak_interval_source_data.csv``
and ``q85_source_manifest.json`` (overridable by command-line arguments).
Outputs: PNG/SVG figure, descriptive summary CSV, input-table copy, legend and
run manifest in a new output directory.
Example: ``python code/render_extended_data_figure1.py --output-dir outputs/extended_data_figure1``
Fixed assumptions: Q85 is Delta_minCMV >= 2.260200281716645; the stronger
subset is >=1.6205657816946992; subsets are nested; the displayed x range ends
at 145 residues although the 196-residue reference observation is retained.
The script does not reconstruct peak membership or deletion classes and never
overwrites an existing output directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib_extended_data_interpeak"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "xdg_extended_data_interpeak"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from PIL import Image  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "data/extended_data_figure1/interpeak_interval_source_data.csv"
DEFAULT_Q85_MANIFEST = ROOT / "data/extended_data_figure1/q85_source_manifest.json"

ORDER = ["GoF Q85-targeted", "GoF +0.5-targeted", "GoF-targeted", "All interpeak"]
Y_POSITIONS = {
    "GoF Q85-targeted": 3,
    "GoF +0.5-targeted": 2,
    "GoF-targeted": 1,
    "All interpeak": 0,
}
DISPLAY_LABELS = {
    "GoF Q85-targeted": "Top 15% activation-enhancing",
    "GoF +0.5-targeted": "Stronger activation-enhancing",
    "GoF-targeted": "Activation-enhancing",
    "All interpeak": "All interpeak intervals\n(reference)",
}

# Same manuscript palette as the finalized main panels.
COLORS = {
    "GoF Q85-targeted": "#1B9093",
    "GoF +0.5-targeted": "#56AAA7",
    "GoF-targeted": "#A3D3CF",
    "All interpeak": "#747C84",
}
POINT_EDGES = {
    "GoF Q85-targeted": "#126F72",
    "GoF +0.5-targeted": "#347E7B",
    "GoF-targeted": "#4E9590",
    "All interpeak": "none",
}
INK = "#25282B"
GRID = "#D7DDE0"
Q85_THRESHOLD = 2.260200281716645
STRONGER_THRESHOLD = 1.6205657816946992

EXPECTED_SUMMARY = {
    "All interpeak": {"n": 283, "median": 28.0, "q1": 15.5, "q3": 42.0, "max": 196.0},
    "GoF-targeted": {"n": 39, "median": 43.0, "q1": 32.0, "q3": 57.5, "max": 134.0},
    "GoF +0.5-targeted": {"n": 14, "median": 47.0, "q1": 38.5, "q3": 62.75, "max": 134.0},
    "GoF Q85-targeted": {"n": 5, "median": 43.0, "q1": 40.0, "q3": 51.0, "max": 52.0},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-csv", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--thresholds-json", type=Path, default=DEFAULT_Q85_MANIFEST)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.8,
            "axes.labelsize": 9.6,
            "xtick.labelsize": 8.2,
            "ytick.labelsize": 8.5,
            "axes.linewidth": 0.8,
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def summarize(source: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for series in ORDER:
        values = source.loc[source["series"].eq(series), "interpeak_interval_length"].dropna().to_numpy(float)
        if not len(values):
            raise ValueError(f"Missing required interval series: {series}")
        q1, median, q3 = np.percentile(values, [25, 50, 75])
        rows.append(
            {
                "series": series,
                "display_label": DISPLAY_LABELS[series].replace("\n", " "),
                "n": int(len(values)),
                "median_residues": float(median),
                "q1_residues": float(q1),
                "q3_residues": float(q3),
                "min_residues": float(np.min(values)),
                "max_residues": float(np.max(values)),
            }
        )
    return pd.DataFrame(rows)


def validate_summary(summary: pd.DataFrame) -> None:
    indexed = summary.set_index("series")
    for series, expected in EXPECTED_SUMMARY.items():
        observed = indexed.loc[series]
        checks = {
            "n": int(observed["n"]),
            "median": float(observed["median_residues"]),
            "q1": float(observed["q1_residues"]),
            "q3": float(observed["q3_residues"]),
            "max": float(observed["max_residues"]),
        }
        if checks != expected:
            raise ValueError(f"Unexpected summary for {series}: {checks}; expected {expected}")


def render(source: pd.DataFrame, summary: pd.DataFrame, output_dir: Path) -> list[Path]:
    fig, ax = plt.subplots(figsize=(8.25, 4.65))
    rng = np.random.default_rng(20260901)
    indexed = summary.set_index("series")
    plot_order = ["All interpeak", "GoF-targeted", "GoF +0.5-targeted", "GoF Q85-targeted"]
    positions = np.arange(len(plot_order), dtype=float)
    values_by_series = [
        source.loc[source["series"].eq(series), "interpeak_interval_length"].dropna().to_numpy(float)
        for series in plot_order
    ]
    violins = ax.violinplot(
        values_by_series,
        positions=positions,
        vert=False,
        widths=0.72,
        showextrema=False,
        showmeans=False,
        showmedians=False,
    )
    violin_alpha = {
        "All interpeak": 0.27,
        "GoF-targeted": 0.31,
        "GoF +0.5-targeted": 0.35,
        "GoF Q85-targeted": 0.39,
    }

    for body, series, y, values in zip(violins["bodies"], plot_order, positions, values_by_series):
        body.set_facecolor(COLORS[series])
        body.set_edgecolor("none")
        body.set_alpha(violin_alpha[series])
        row = indexed.loc[series]
        q1 = float(row["q1_residues"])
        median = float(row["median_residues"])
        q3 = float(row["q3_residues"])
        ax.plot(
            [q1, q3],
            [y, y],
            color=COLORS[series],
            linewidth=8.0,
            solid_capstyle="butt",
            alpha=1.0,
            zorder=4,
        )
        ax.scatter(
            [median],
            [y],
            s=62,
            facecolor="white",
            edgecolor=COLORS[series],
            linewidth=1.35,
            zorder=5,
        )
        if len(values) <= 15:
            jitter = rng.uniform(-0.065, 0.065, len(values))
            ax.scatter(
                values,
                y + jitter,
                s=20,
                color=COLORS[series],
                edgecolor=POINT_EDGES[series],
                linewidth=0.4,
                alpha=0.92,
                zorder=4,
            )
        ax.text(
            q3 + 4.0,
            y + 0.39,
            f"median {median:.0f} residues  (n={int(row['n'])})",
            color="#4F575E",
            fontsize=8.2,
            fontweight="bold",
            va="bottom",
            ha="left",
            zorder=6,
        )

    ax.set_xlim(0, 145)
    ax.set_xticks(np.arange(0, 141, 20))
    ax.set_ylim(-0.55, 3.62)
    ax.set_yticks(
        positions,
        [DISPLAY_LABELS[series] for series in plot_order],
    )
    ax.set_xlabel("Interpeak interval length (residues)")
    ax.grid(axis="x", color=GRID, linewidth=0.55, alpha=0.52)
    ax.set_axisbelow(True)
    ax.tick_params(length=3.0, width=0.7)

    fig.text(
        0.275,
        0.955,
        "Interpeak interval lengths associated with activation-enhancing deletions",
        fontsize=11.5,
        fontweight="bold",
        ha="left",
        va="top",
    )
    fig.subplots_adjust(left=0.275, right=0.985, bottom=0.16, top=0.84)

    png = output_dir / "Extended_Data_Fig_1_interpeak_interval_lengths_20260901.png"
    svg = output_dir / "Extended_Data_Fig_1_interpeak_interval_lengths_20260901.svg"
    fig.savefig(png, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    with Image.open(png) as image:
        if image.mode != "RGB":
            image.convert("RGB").save(png, dpi=(600, 600))
    return [png, svg]


def main() -> int:
    args = parse_args()
    source_path = args.source_csv.expanduser().resolve()
    thresholds_path = args.thresholds_json.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    setup_style()
    for path in (source_path, thresholds_path):
        if not path.exists():
            raise FileNotFoundError(path)
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output directory: {output_dir}")
    output_dir.mkdir(parents=True)

    source = pd.read_csv(source_path)
    q85_manifest = json.loads(thresholds_path.read_text())
    recovered_q85 = float(q85_manifest["q85_definition"]["threshold_Delta_minCMV"])
    if recovered_q85 != Q85_THRESHOLD:
        raise ValueError(f"Unexpected Q85 threshold: {recovered_q85}")
    recovered_stronger = float(
        q85_manifest["stronger_activation_enhancing_threshold_Delta_minCMV"]
    )
    if recovered_stronger != STRONGER_THRESHOLD:
        raise ValueError(f"Unexpected stronger-subset threshold: {recovered_stronger}")

    summary = summarize(source)
    validate_summary(summary)
    outputs = render(source, summary, output_dir)

    source.to_csv(output_dir / "Extended_Data_Fig_1_source_data_20260901.csv", index=False)
    summary.to_csv(output_dir / "Extended_Data_Fig_1_summary_20260901.csv", index=False)

    legend = (
        "Extended Data Fig. 1 | Interpeak interval lengths associated with activation-enhancing deletions. "
        "Distributions of interpeak interval lengths are shown for all interpeak intervals (reference; "
        "n=283), intervals associated with activation-enhancing deletions (n=39), a stronger "
        "activation-enhancing subset (n=14), and the top 15% activation-enhancing subset (n=5). "
        "The stronger subset used Delta_minCMV >= 1.6205657816946992, corresponding to 0.5 above "
        "the activation-enhancing threshold; the top-15% subset used Delta_minCMV >= "
        "2.260200281716645 (the original Q85 definition). Violins show interval-length "
        "distributions, thick horizontal segments show interquartile ranges, open circles show "
        "medians, and points show individual observations for the two smallest subsets. "
        "Activation-enhancing subsets are nested and the top-15% subset is small; comparisons are "
        "descriptive and do not establish an optimal linker length or mechanism. The display range "
        "is limited to 145 residues to match the historical panel; one 196-residue reference "
        "interval lies outside the displayed range and is retained in the source-data file.\n"
    )
    legend_path = output_dir / "Extended_Data_Fig_1_legend_20260901.txt"
    legend_path.write_text(legend)

    manifest = {
        "created_date": "2026-09-01",
        "rendering_only": True,
        "source": {"path": str(source_path), "sha256": sha256(source_path)},
        "source_lineage": (
            "Latest validated interval table from the 2026-08-09 r2 Q85 sensitivity replot; "
            "the 2026-08-25 terminology-cleanup renders reused the same interval values."
        ),
        "display_labels": DISPLAY_LABELS,
        "palette": COLORS,
        "subset_definitions": {
            "stronger_activation_enhancing_threshold_Delta_minCMV": STRONGER_THRESHOLD,
            "top_15_percent_threshold_Delta_minCMV": Q85_THRESHOLD,
            "top_15_percent_internal_legacy_name": "Q85",
        },
        "summary": summary.to_dict(orient="records"),
        "rendering_note": (
            "The original horizontal violin encoding, IQR bars, median circles, small-subset "
            "points, and 0-145-residue display range were restored from the historical renderer. "
            "Median and sample-size annotations are positioned above, rather than over, each "
            "distribution, and no extended-data figure number is embedded in the plotting area. "
            "One 196-residue all-interpeak reference interval is retained in source data but lies "
            "outside the historical display range."
        ),
        "outputs": [str(path) for path in outputs] + [str(legend_path)],
        "validation": {
            "all_outputs_exist": all(path.exists() for path in outputs + [legend_path]),
            "all_outputs_nonempty": all(path.stat().st_size > 0 for path in outputs + [legend_path]),
            "source_observations_above_display_domain": int((source["interpeak_interval_length"] > 145).sum()),
            "png_mode": Image.open(outputs[0]).mode,
            "embedded_figure_number": False,
        },
    }
    (output_dir / "manifest_20260901.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
