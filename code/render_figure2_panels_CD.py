#!/usr/bin/env python3
"""Render the quantitative panels in manuscript Figure 2C/D.

Purpose: redraw frozen deletion-centered contact profiles and WFYL-zero
proportions; only a five-residue centered display mean is applied to Panel C.
Manuscript support: Figure 2C (contact profiles) and Figure 2D (WFYL-zero).
Required inputs: ``data/figure2/contact_profile_source_data.csv`` and
``wfyl_source_data.csv`` (overridable by command-line arguments).
Outputs: PNG/SVG panels, source-table copies, legend and run manifest in a new
output directory.
Example: ``python code/render_figure2_panels_CD.py --output-dir outputs/figure2_panels_CD``
Fixed assumptions: Panel C activation-impairing means LoF OR necessary and
uses Q85 Delta_minCMV >= 2.260200281716645; Panel D uses LoF only and displays
the supplied gene-bootstrap 95% intervals. This script does not recalculate
phenotypes, contact profiles, WFYL status or confidence intervals and never
overwrites an existing output directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib_figure2_CD"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "xdg_figure2_CD"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from PIL import Image  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTACT_SOURCE = ROOT / "data/figure2/contact_profile_source_data.csv"
DEFAULT_WFYL_SOURCE = ROOT / "data/figure2/wfyl_source_data.csv"
Q85_THRESHOLD = 2.260200281716645
STRONGER_THRESHOLD = 1.6205657816946992
DISPLAY_SMOOTHING_WINDOW = 5

# Muted, colorblind-compatible publication palette. The teal family identifies
# nested activation-enhancing sets; orange identifies the opposing phenotype.
PALETTE = {
    "neutral": "#747C84",
    "activation_enhancing_light": "#A3D3CF",
    "activation_enhancing_mid": "#56AAA7",
    "activation_enhancing_dark": "#1B9093",
    "activation_impairing": "#E8892C",
    "ink": "#25282B",
    "whisker": "#596168",
    "grid": "#D7DDE0",
    "reference": "#949CA2",
}

CONTACT_ORDER = ["GoF Q85", "Neutral", "LoF + Necessary"]
CONTACT_DISPLAY = {
    "GoF Q85": "Top 15% activation-enhancing deletions",
    "Neutral": "Neutral deletions",
    "LoF + Necessary": "Activation-impairing deletions",
}
CONTACT_COLORS = {
    "GoF Q85": PALETTE["activation_enhancing_dark"],
    "Neutral": PALETTE["neutral"],
    "LoF + Necessary": PALETTE["activation_impairing"],
}
CONTACT_LINEWIDTHS = {"GoF Q85": 2.15, "Neutral": 1.70, "LoF + Necessary": 2.20}

WFYL_ORDER = ["Neutral", "GoF", "GoF +0.5", "GoF Q85", "LoF"]
WFYL_AXIS_LABELS = {
    "Neutral": "Neutral",
    "GoF": "Activation-\nenhancing",
    "GoF +0.5": "Stronger activation-\nenhancing",
    "GoF Q85": "Top 15% activation-\nenhancing",
    "LoF": "Activation-\nimpairing",
}
WFYL_COLORS = {
    "Neutral": PALETTE["neutral"],
    "GoF": PALETTE["activation_enhancing_light"],
    "GoF +0.5": PALETTE["activation_enhancing_mid"],
    "GoF Q85": PALETTE["activation_enhancing_dark"],
    "LoF": PALETTE["activation_impairing"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contact-source-csv", type=Path, default=DEFAULT_CONTACT_SOURCE)
    parser.add_argument("--wfyl-source-csv", type=Path, default=DEFAULT_WFYL_SOURCE)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.5,
            "axes.titlesize": 11.0,
            "axes.labelsize": 9.2,
            "xtick.labelsize": 7.8,
            "ytick.labelsize": 8.0,
            "axes.linewidth": 0.8,
            "axes.edgecolor": PALETTE["ink"],
            "axes.labelcolor": PALETTE["ink"],
            "xtick.color": PALETTE["ink"],
            "ytick.color": PALETTE["ink"],
            "text.color": PALETTE["ink"],
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


def centered_mean(values: np.ndarray, window: int = DISPLAY_SMOOTHING_WINDOW) -> np.ndarray:
    return pd.Series(values).rolling(window=window, center=True, min_periods=1).mean().to_numpy(float)


def add_panel_letter(ax: plt.Axes, letter: str, *, y: float = 1.105) -> None:
    ax.text(
        -0.13,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=25,
        fontweight="bold",
        ha="left",
        va="top",
        clip_on=False,
    )


def save_figure(fig: plt.Figure, stem: str, output_dir: Path) -> list[Path]:
    png = output_dir / f"{stem}.png"
    svg = output_dir / f"{stem}.svg"
    fig.savefig(png, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    with Image.open(png) as image:
        if image.mode != "RGB":
            image.convert("RGB").save(png, dpi=(600, 600))
    return [png, svg]


def render_panel_c(source: pd.DataFrame, output_dir: Path) -> tuple[list[Path], dict[str, int]]:
    fig, ax = plt.subplots(figsize=(7.7, 3.65))
    sample_sizes: dict[str, int] = {}

    for series in CONTACT_ORDER:
        frame = source.loc[source["series"].eq(series)].sort_values("offset")
        if frame.empty:
            raise ValueError(f"Panel C source is missing required series: {series}")
        x = frame["offset"].to_numpy(float)
        mean = centered_mean(frame["mean_contact_percentile"].to_numpy(float))
        sample_sizes[series] = int(frame["n_windows"].iloc[0])
        ax.plot(
            x,
            mean,
            color=CONTACT_COLORS[series],
            linewidth=CONTACT_LINEWIDTHS[series],
            solid_capstyle="round",
            zorder=4,
        )

    ax.axvline(
        0,
        color=PALETTE["reference"],
        linewidth=0.75,
        linestyle=(0, (2.5, 2.5)),
        zorder=2,
    )
    ax.set_xlim(-30, 30)
    ax.set_ylim(0.38, 0.72)
    ax.set_xticks(np.arange(-30, 31, 10))
    ax.set_yticks(np.arange(0.40, 0.71, 0.10))
    ax.set_xlabel("Position relative to deletion center (residues)")
    ax.set_ylabel("Mean contact percentile")
    ax.set_title(
        "Contact profiles around deletion sites",
        loc="left",
        y=1.205,
        fontweight="bold",
        pad=0,
    )
    ax.grid(axis="y", color=PALETTE["grid"], linewidth=0.55, alpha=0.42)
    ax.set_axisbelow(True)
    ax.tick_params(length=3.0, width=0.7)
    add_panel_letter(ax, "C", y=1.265)

    legend_labels = {
        "GoF Q85": "Top 15% activation-enhancing",
        "Neutral": "Neutral",
        "LoF + Necessary": "Activation-impairing",
    }
    handles = [
        Line2D(
            [0],
            [0],
            color=CONTACT_COLORS[series],
            linewidth=CONTACT_LINEWIDTHS[series],
            solid_capstyle="round",
            label=legend_labels[series],
        )
        for series in CONTACT_ORDER
    ]
    legend = ax.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.0, 1.025, 1.0, 0.10),
        mode="expand",
        ncol=3,
        frameon=False,
        borderaxespad=0,
        handlelength=1.75,
        handletextpad=0.55,
        columnspacing=1.1,
        fontsize=7.6,
    )
    for text, series in zip(legend.get_texts(), CONTACT_ORDER):
        text.set_color(CONTACT_COLORS[series])
        text.set_fontweight("bold")

    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.20, top=0.73)
    return save_figure(fig, "Figure2_panel_C_contact_landscape", output_dir), sample_sizes


def render_panel_d(source: pd.DataFrame, output_dir: Path) -> tuple[list[Path], list[dict[str, float | int | str]]]:
    summary = source.set_index("series").loc[WFYL_ORDER].reset_index()
    x = np.arange(len(summary), dtype=float)
    values = summary["percent"].to_numpy(float)
    lower = summary["cluster_bootstrap_ci95_low_percent"].to_numpy(float)
    upper = summary["cluster_bootstrap_ci95_high_percent"].to_numpy(float)

    fig, ax = plt.subplots(figsize=(8.15, 4.55))
    bars = ax.bar(
        x,
        values,
        width=0.64,
        color=[WFYL_COLORS[series] for series in summary["series"]],
        edgecolor="white",
        linewidth=0.8,
        zorder=3,
    )

    errors = np.vstack([values - lower, upper - values])
    ax.errorbar(
        x,
        values,
        yerr=errors,
        fmt="none",
        ecolor=PALETTE["whisker"],
        elinewidth=0.95,
        capsize=3.0,
        capthick=0.95,
        zorder=5,
    )

    labels = [WFYL_AXIS_LABELS[row.series] for row in summary.itertuples(index=False)]
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 75)
    ax.set_yticks(np.arange(0, 71, 10))
    ax.set_ylabel("Deleted regions lacking W/F/Y/L (%)")
    ax.set_title(
        "Activation-enhancing deletions preferentially remove regions lacking W/F/Y/L",
        loc="left",
        fontweight="bold",
        pad=8,
    )
    ax.grid(axis="y", color=PALETTE["grid"], linewidth=0.55, alpha=0.48)
    ax.set_axisbelow(True)
    ax.tick_params(length=3.0, width=0.7)
    ax.tick_params(axis="x", labelsize=8.4)
    add_panel_letter(ax, "D")

    for bar, value, upper_limit in zip(bars, values, upper):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            upper_limit + 1.55,
            f"{value:.1f}%",
            ha="center",
            va="bottom",
            fontsize=8.1,
            fontweight="bold",
        )

    fig.subplots_adjust(left=0.115, right=0.985, bottom=0.21, top=0.84)
    records = [
        {
            "series": str(row.series),
            "display_label": WFYL_AXIS_LABELS[str(row.series)].replace("\n", " "),
            "n_windows": int(row.n_windows),
            "percent": float(row.percent),
            "ci95_low_percent": float(row.cluster_bootstrap_ci95_low_percent),
            "ci95_high_percent": float(row.cluster_bootstrap_ci95_high_percent),
        }
        for row in summary.itertuples(index=False)
    ]
    return save_figure(fig, "Figure2_panel_D_WFYL_poor_deletions", output_dir), records


def validate_sources(contact: pd.DataFrame, wfyl: pd.DataFrame) -> None:
    required_contact = {"series", "offset", "mean_contact_percentile", "n_windows"}
    required_wfyl = {
        "series",
        "n_windows",
        "n_genes",
        "wfyl_free_count",
        "percent",
        "cluster_bootstrap_ci95_low_percent",
        "cluster_bootstrap_ci95_high_percent",
    }
    if missing := sorted(required_contact - set(contact.columns)):
        raise ValueError(f"Panel C source is missing columns: {missing}")
    if missing := sorted(required_wfyl - set(wfyl.columns)):
        raise ValueError(f"Panel D source is missing columns: {missing}")
    expected_contact_n = {"GoF Q85": 38, "Neutral": 1838, "LoF + Necessary": 1643}
    observed_contact_n = {
        series: int(contact.loc[contact["series"].eq(series), "n_windows"].iloc[0])
        for series in CONTACT_ORDER
    }
    if observed_contact_n != expected_contact_n:
        raise ValueError(f"Unexpected Panel C sample sizes: {observed_contact_n}")

    expected_wfyl_n = {"Neutral": 1838, "GoF": 234, "GoF +0.5": 101, "GoF Q85": 38, "LoF": 1341}
    observed_wfyl_n = {
        str(row.series): int(row.n_windows)
        for row in wfyl.set_index("series").loc[WFYL_ORDER].reset_index().itertuples(index=False)
    }
    if observed_wfyl_n != expected_wfyl_n:
        raise ValueError(f"Unexpected Panel D sample sizes: {observed_wfyl_n}")


def main() -> int:
    args = parse_args()
    contact_source = args.contact_source_csv.expanduser().resolve()
    wfyl_source = args.wfyl_source_csv.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    setup_style()
    for path in (contact_source, wfyl_source):
        if not path.exists():
            raise FileNotFoundError(path)
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output directory: {output_dir}")
    output_dir.mkdir(parents=True)

    contact = pd.read_csv(contact_source)
    wfyl = pd.read_csv(wfyl_source)
    validate_sources(contact, wfyl)

    panel_c_outputs, panel_c_n = render_panel_c(contact, output_dir)
    panel_d_outputs, panel_d_values = render_panel_d(wfyl, output_dir)
    contact.to_csv(output_dir / "Figure2_panel_C_source_data.csv", index=False)
    wfyl.set_index("series").loc[WFYL_ORDER].reset_index().to_csv(
        output_dir / "Figure2_panel_D_source_data.csv", index=False
    )

    legend = (
        "Figure 2C. Deletion-centered mean contact-percentile profiles (five-residue centered "
        "display mean) for neutral (n=1,838), activation-impairing (n=1,643), and the top "
        "15% of mapped activation-enhancing deletions (n=38). The top-15% subset corresponds "
        "to the original Q85 definition, Delta_minCMV >= 2.260200281716645.\n\n"
        "Figure 2D. Percentage of mapped/analyzed deleted regions lacking W, F, Y, or L: "
        "neutral (n=1,838), activation-enhancing (n=234), stronger activation-enhancing "
        "(n=101), top 15% activation-enhancing (n=38), and activation-impairing (n=1,341). "
        "The stronger activation-enhancing subset uses Delta_minCMV >= 1.6205657816946992 "
        "(0.5 above the activation-enhancing threshold); the top-15% subset uses "
        "Delta_minCMV >= 2.260200281716645. Whiskers show two-sided gene-cluster bootstrap "
        "95% confidence intervals.\n"
    )
    (output_dir / "Figure2_panels_CD_legend.txt").write_text(legend)

    outputs = panel_c_outputs + panel_d_outputs
    manifest = {
        "rendering_only": True,
        "visible_terminology": "activation-enhancing / activation-impairing deletions",
        "palette": PALETTE,
        "source_files": {
            "panel_C": {"path": str(contact_source), "sha256": sha256(contact_source)},
            "panel_D": {"path": str(wfyl_source), "sha256": sha256(wfyl_source)},
        },
        "panel_C": {
            "series_internal_order": CONTACT_ORDER,
            "display_labels": CONTACT_DISPLAY,
            "sample_sizes": panel_c_n,
            "display_smoothing": f"{DISPLAY_SMOOTHING_WINDOW}-residue centered moving mean, min_periods=1",
            "q85_threshold_Delta_minCMV": Q85_THRESHOLD,
        },
        "panel_D": {
            "series_internal_order": WFYL_ORDER,
            "display_labels": {key: value.replace("\n", " ") for key, value in WFYL_AXIS_LABELS.items()},
            "values": panel_d_values,
            "whisker_display": "full two-sided gene-cluster bootstrap 95% confidence intervals",
            "stronger_threshold_Delta_minCMV": STRONGER_THRESHOLD,
            "q85_threshold_Delta_minCMV": Q85_THRESHOLD,
        },
        "outputs": [str(path) for path in outputs],
        "validation": {
            "all_outputs_exist": all(path.exists() for path in outputs),
            "all_outputs_nonempty": all(path.stat().st_size > 0 for path in outputs),
            "png_modes": {
                path.name: Image.open(path).mode for path in outputs if path.suffix.lower() == ".png"
            },
        },
    }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
