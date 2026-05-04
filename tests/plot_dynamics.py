#!/usr/bin/env python3
"""Plot ADAF dynamics radial profiles for smoke-test inspection."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Callable, List

from compare_dynamics import COLUMNS, read_profile


class Panel:
    def __init__(
        self,
        column: str,
        label: str,
        transform: Callable[[float], float] = lambda value: value,
    ) -> None:
        self.column = column
        self.label = label
        self.transform = transform


PANELS = [
    Panel("mach", r"Mach number $v_R/c_s$ [dimensionless]"),
    Panel("log_te", r"$T_e$ [K]", lambda value: 10.0**value),
    Panel("log_ti", r"$T_i$ [K]", lambda value: 10.0**value),
    Panel("q_advi_over_q_vis", r"$q_{adv,i}/q_{vis}$ [dimensionless]"),
    Panel("q_adv_total_over_q_vis", r"$(q_{adv,i}+q_{adv,e})/q_{vis}$ [dimensionless]"),
    Panel("cs_over_c", r"$c_s/c$ [dimensionless]"),
    Panel("height", r"$H$ [$R_S$]"),
    Panel("q_rad", r"$q_{rad}$ [erg s$^{-1}$ cm$^{-3}$]"),
    Panel("tau", r"$\tau$ [dimensionless]"),
    Panel("log_lk", r"$l_K$ [code units]", lambda value: 10.0**value),
    Panel("log_l", r"$l$ [code units]", lambda value: 10.0**value),
    Panel("magnetic_field", r"Magnetic field $B$ [G]"),
    Panel("log_rho", r"$\rho$ [g cm$^{-3}$]", lambda value: 10.0**value),
    Panel("log_ne", r"$n_e$ [cm$^{-3}$]", lambda value: 10.0**value),
]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a gridded dynamics profile plot."
    )
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(args.output.parent / ".matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(args.output.parent / ".cache"))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    reference = rows_by_radius(read_profile(args.reference))
    candidate = rows_by_radius(read_profile(args.candidate))

    fig, axes = plt.subplots(4, 4, figsize=(14, 10.5), sharex=False)
    axes_list = list(axes.flat)

    for axis, panel in zip(axes_list, PANELS):
        plot_panel(axis, panel, reference, "fiducial", "--", "tab:blue")
        plot_panel(axis, panel, candidate, "current", "-", "tab:orange")
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_title(panel.label, fontsize=9)
        axis.set_xlabel(r"Radius $R$ [$R_S$], log scale")
        axis.set_ylabel("Value, log scale")
        axis.grid(True, which="both", alpha=0.25)

    for axis in axes_list[len(PANELS) :]:
        axis.set_visible(False)

    axes_list[0].legend(loc="best", fontsize="small")
    fig.suptitle("ADAF dynamics radial profiles: examples/largeR.dat", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(args.output, dpi=160)
    plt.close(fig)

    print(f"Wrote dynamics profile plot to {args.output}")
    return 0


def rows_by_radius(rows: List[List[float]]) -> List[List[float]]:
    radius_index = COLUMNS.index("radius")
    return sorted(rows, key=lambda row: row[radius_index])


def plot_panel(axis, panel: Panel, rows: List[List[float]], label: str, style: str, color: str) -> None:
    radius_index = COLUMNS.index("radius")
    value_index = COLUMNS.index(panel.column)
    x_values = [row[radius_index] for row in rows]
    y_values = [panel.transform(row[value_index]) for row in rows]

    if any(value <= 0.0 for value in x_values + y_values):
        raise SystemExit(f"{panel.column} contains non-positive values; cannot use log axes.")

    axis.plot(x_values, y_values, style, color=color, linewidth=1.5, label=label)


if __name__ == "__main__":
    raise SystemExit(main())
