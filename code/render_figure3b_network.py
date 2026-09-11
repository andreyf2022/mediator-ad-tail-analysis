#!/usr/bin/env python3
"""Render the receptor-site co-engagement network in Figure 3B.

Purpose: redraw the accepted network from frozen node, edge and 2-D position
tables; no coordinates or contact pairs are read.
Manuscript support: Figure 3B.
Required inputs: ``data/figure3/figure3b_nodes.csv``, ``figure3b_edges.csv`` and
``figure3b_positions.csv`` (overridable by command-line arguments).
Outputs: ``hotspot_receptor_graph.png`` and ``.svg`` in a new output directory.
Example: ``python code/render_figure3b_network.py --out-dir outputs/figure3b_network``
Fixed assumptions: eight displayed receptor sites (including MED24 D653),
1,139 usable models, node area proportional to strict contact mass, edge width
proportional to same-model strict co-engagement fraction, and the accepted
fixed node positions/three display-only edge curvatures. Existing output
directories are never overwritten.
"""

from __future__ import annotations

import argparse
import math
import os
import tempfile
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mpl_hotspot_receptor_graph_only")
)
os.environ.setdefault(
    "XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "fontconfig_hotspot_receptor_graph_only")
)

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch


EXPECTED_HOTSPOTS = (
    "MED15_L76",
    "MED16_V248",
    "MED16_T601",
    "MED24_V457",
    "MED24_D653",
    "MED23_L343",
    "MED23_F884",
    "MED23_Y1303",
)
EXPECTED_USABLE_MODELS = 1139


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--nodes-csv",
        type=Path,
        default=PACKAGE_ROOT / "data/figure3/figure3b_nodes.csv",
    )
    parser.add_argument(
        "--edges-csv",
        type=Path,
        default=PACKAGE_ROOT / "data/figure3/figure3b_edges.csv",
    )
    parser.add_argument(
        "--positions-csv",
        type=Path,
        default=PACKAGE_ROOT / "data/figure3/figure3b_positions.csv",
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args()


def padded_limits(ext: dict[str, np.ndarray]) -> tuple[tuple[float, float], tuple[float, float]]:
    matrix = np.asarray(list(ext.values()))
    xmin, ymin = matrix.min(axis=0)
    xmax, ymax = matrix.max(axis=0)
    xrange = max(xmax - xmin, 1.0)
    yrange = max(ymax - ymin, 1.0)
    return (
        (xmin - 0.17 * xrange, xmax + 0.17 * xrange),
        (ymin - 0.19 * yrange, ymax + 0.22 * yrange),
    )


def render(
    node: pd.DataFrame,
    edge: pd.DataFrame,
    ext: dict[str, np.ndarray],
    out_dir: Path,
) -> tuple[Path, Path]:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "svg.fonttype": "none",
        }
    )
    red = "#9F2525"
    red_dark = "#7E1717"
    ink = "#151515"

    node = node.set_index("hotspot_key")
    masses = node["strict_contact_mass"].astype(float)
    max_mass = float(masses.max()) or 1.0
    node_scale_factor = 0.78
    node_area = {
        key: node_scale_factor * (95.0 + 910.0 * float(masses.loc[key]) / max_mass)
        for key in masses.index
    }
    max_edge_fraction = float(edge["fraction_usable_models"].max()) if len(edge) else 1.0

    fig, ax = plt.subplots(figsize=(7.4, 5.25))
    xlim, ylim = padded_limits(ext)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")

    fig.text(
        0.065,
        0.955,
        "Extended-state weighted receptor graph",
        fontsize=12,
        fontweight="bold",
        va="top",
    )
    fig.text(
        0.065,
        0.915,
        "strict hydrophobic support • 1,139 artifact-clean usable models",
        fontsize=8.5,
        color="#555555",
        va="top",
    )

    # Edge trajectories are visual connectors, not spatial paths.  Retain the
    # accepted restrained curves used to avoid unrelated projected centroids.
    edge_curvature = {
        frozenset(("MED24_V457", "MED23_Y1303")): 0.06,
        frozenset(("MED16_T601", "MED23_F884")): 0.08,
        frozenset(("MED16_T601", "MED24_V457")): 0.06,
    }
    for row in edge.sort_values("supporting_models").itertuples():
        p1 = ext[row.hotspot_1]
        p2 = ext[row.hotspot_2]
        support_fraction = float(row.fraction_usable_models) / max_edge_fraction
        width = 0.18 + 8.27 * support_fraction
        alpha = 0.16 + 0.58 * math.sqrt(support_fraction)
        curve = edge_curvature.get(frozenset((row.hotspot_1, row.hotspot_2)), 0.0)
        ax.add_patch(
            FancyArrowPatch(
                posA=p1,
                posB=p2,
                arrowstyle="-",
                connectionstyle=f"arc3,rad={curve}",
                color=red_dark,
                linewidth=width,
                alpha=alpha,
                capstyle="round",
                joinstyle="round",
                zorder=1,
            )
        )

    label_offsets = {
        "MED15_L76": (11, 15, "left", "bottom"),
        "MED16_V248": (-19, -5, "right", "top"),
        "MED16_T601": (-10, 15, "right", "bottom"),
        "MED24_V457": (-12, -12, "right", "top"),
        "MED24_D653": (12, -15, "left", "top"),
        "MED23_L343": (18, -4, "left", "top"),
        "MED23_F884": (15, -2, "left", "center"),
        "MED23_Y1303": (8, 15, "left", "bottom"),
    }
    leader_keys = {
        "MED16_V248",
        "MED16_T601",
        "MED24_V457",
        "MED24_D653",
        "MED23_L343",
        "MED15_L76",
    }
    for hotspot in node.reset_index().itertuples():
        point = ext[hotspot.hotspot_key]
        ax.scatter(
            [point[0]],
            [point[1]],
            s=node_area[hotspot.hotspot_key],
            color=red,
            edgecolor="white",
            linewidth=1.25,
            zorder=3,
        )
        dx, dy, ha, va = label_offsets[hotspot.hotspot_key]
        ax.annotate(
            f"{hotspot.subunit}\n{hotspot.residue_label}",
            xy=point,
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha,
            va=va,
            fontsize=9.2,
            color=ink,
            linespacing=1.12,
            arrowprops=(
                {
                    "arrowstyle": "-",
                    "color": "#6E6E6E",
                    "linewidth": 0.55,
                    "shrinkA": 2,
                    "shrinkB": 5,
                }
                if hotspot.hotspot_key in leader_keys
                else None
            ),
            zorder=5,
        )

    example_masses = np.quantile(masses.to_numpy(), [0.5, 1.0])
    example_sizes = [
        node_scale_factor * (95.0 + 910.0 * value / max_mass)
        for value in example_masses
    ]
    edge_values = (
        np.quantile(edge["supporting_models"].to_numpy(), [0.5, 1.0])
        if len(edge)
        else np.array([0, 0])
    )
    edge_legend_styles = []
    for value in edge_values:
        fraction = float(value) / EXPECTED_USABLE_MODELS
        support_fraction = fraction / max_edge_fraction
        width = 0.18 + 8.27 * support_fraction
        alpha = 0.16 + 0.58 * math.sqrt(support_fraction)
        edge_legend_styles.append((width, alpha))

    # Compact, manually spaced key keeps marker areas and labels independent.
    legend_transform = ax.transAxes
    legend_fontsize = 7.6
    ax.text(
        0.02,
        -0.105,
        "Node area = strict mass",
        transform=legend_transform,
        ha="left",
        va="center",
        fontsize=legend_fontsize,
        color=ink,
        clip_on=False,
    )
    for x, size, value, label_x in zip(
        (0.08, 0.31), example_sizes, example_masses, (0.12, 0.38)
    ):
        ax.scatter(
            [x],
            [-0.205],
            s=size,
            transform=legend_transform,
            color=red,
            edgecolor="white",
            linewidth=0.8,
            clip_on=False,
            zorder=6,
        )
        ax.text(
            label_x,
            -0.205,
            f"{value:.0f}",
            transform=legend_transform,
            ha="left",
            va="center",
            fontsize=legend_fontsize,
            color=ink,
            clip_on=False,
        )

    ax.text(
        0.56,
        -0.105,
        "Edge width = model count",
        transform=legend_transform,
        ha="left",
        va="center",
        fontsize=legend_fontsize,
        color=ink,
        clip_on=False,
    )
    for x0, x1, value, (width, alpha), label_x in zip(
        (0.56, 0.79),
        (0.65, 0.88),
        edge_values,
        edge_legend_styles,
        (0.67, 0.90),
    ):
        ax.plot(
            [x0, x1],
            [-0.205, -0.205],
            transform=legend_transform,
            color=red_dark,
            linewidth=width,
            alpha=alpha,
            solid_capstyle="round",
            clip_on=False,
        )
        ax.text(
            label_x,
            -0.205,
            f"{value:.0f}",
            transform=legend_transform,
            ha="left",
            va="center",
            fontsize=legend_fontsize,
            color=ink,
            clip_on=False,
        )

    fig.subplots_adjust(left=0.065, right=0.975, top=0.88, bottom=0.205)
    out_dir.mkdir(parents=True, exist_ok=False)
    png = out_dir / "hotspot_receptor_graph.png"
    svg = out_dir / "hotspot_receptor_graph.svg"
    fig.savefig(png, dpi=600, facecolor="white", bbox_inches="tight")
    fig.savefig(svg, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return png, svg


def main() -> int:
    args = parse_args()
    out_dir = args.out_dir.resolve()
    if out_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output directory: {out_dir}")

    required = [args.nodes_csv, args.edges_csv, args.positions_csv]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))

    node = pd.read_csv(args.nodes_csv, float_precision="round_trip")
    edge = pd.read_csv(args.edges_csv, float_precision="round_trip")
    positions = pd.read_csv(args.positions_csv, float_precision="round_trip")
    required_node = {
        "hotspot_key",
        "subunit",
        "residue_label",
        "strict_contact_mass",
        "usable_model_denominator",
    }
    required_edge = {
        "hotspot_1",
        "hotspot_2",
        "supporting_models",
        "fraction_usable_models",
        "usable_model_denominator",
    }
    required_position = {"hotspot_key", "x", "y"}
    for name, table, required_columns in (
        ("nodes", node, required_node),
        ("edges", edge, required_edge),
        ("positions", positions, required_position),
    ):
        if missing_columns := sorted(required_columns - set(table.columns)):
            raise ValueError(f"{name} table is missing columns: {missing_columns}")
    for name, table in (("nodes", node), ("positions", positions)):
        keys = tuple(table["hotspot_key"].astype(str))
        if keys != EXPECTED_HOTSPOTS:
            raise ValueError(f"Unexpected {name} hotspot set or ordering: {keys}")
    edge_keys = set(edge["hotspot_1"].astype(str)) | set(edge["hotspot_2"].astype(str))
    if not edge_keys.issubset(EXPECTED_HOTSPOTS):
        raise ValueError(f"Unexpected edge endpoint(s): {sorted(edge_keys - set(EXPECTED_HOTSPOTS))}")
    denominators = set(node["usable_model_denominator"].astype(int)) | set(
        edge["usable_model_denominator"].astype(int)
    )
    if denominators != {EXPECTED_USABLE_MODELS}:
        raise ValueError(f"Unexpected usable-model denominator(s): {sorted(denominators)}")
    expected_fractions = edge["supporting_models"].astype(float) / EXPECTED_USABLE_MODELS
    if not np.allclose(edge["fraction_usable_models"], expected_fractions, rtol=0, atol=1e-15):
        raise ValueError("Edge fractions do not equal supporting_models / 1,139")
    ext = {
        str(row.hotspot_key): np.array([float(row.x), float(row.y)])
        for row in positions.itertuples()
    }
    png, svg = render(node, edge, ext, out_dir)

    print(png)
    print(svg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
