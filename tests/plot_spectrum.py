#!/usr/bin/env python3
"""Plot ADAF spectra for smoke-test inspection."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import List, Tuple

from compare_spectrum import read_spectrum


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a log-space ADAF spectrum comparison plot."
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

    reference = read_spectrum(args.reference)
    candidate = read_spectrum(args.candidate)

    fig, (axis, residual_axis) = plt.subplots(
        2,
        1,
        figsize=(8, 5.25),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [4, 1], "hspace": 0.08},
    )
    plot_spectrum(axis, reference, "fiducial", "--", "tab:blue")
    plot_spectrum(axis, candidate, "current", "-", "tab:orange")
    plot_residuals(residual_axis, reference, candidate)

    axis.set_ylabel(r"$\log_{10}(\nu L_\nu / \mathrm{erg}\ \mathrm{s}^{-1})$")
    axis.set_title("ADAF spectrum: examples/largeR.toml")
    axis.grid(True, alpha=0.25)
    axis.legend(loc="best")
    axis.tick_params(labelbottom=False)

    residual_axis.set_xlabel(r"$\log_{10}(\nu / \mathrm{Hz})$")
    residual_axis.set_ylabel("residual")
    residual_axis.grid(True, alpha=0.25)

    fig.savefig(args.output, dpi=160)
    plt.close(fig)

    print(f"Wrote spectrum plot to {args.output}")
    return 0


def plot_spectrum(
    axis,
    rows: List[Tuple[float, float]],
    label: str,
    style: str,
    color: str,
) -> None:
    x_values = [row[0] for row in rows]
    y_values = [row[1] for row in rows]
    axis.plot(x_values, y_values, style, color=color, linewidth=1.5, label=label)


def plot_residuals(
    axis,
    reference: List[Tuple[float, float]],
    candidate: List[Tuple[float, float]],
) -> None:
    if len(reference) != len(candidate):
        raise SystemExit(
            "Cannot plot spectrum residuals for files with different row counts."
        )

    x_values = [row[0] for row in reference]
    residuals = [
        actual_log_luminosity - expected_log_luminosity
        for (_, expected_log_luminosity), (_, actual_log_luminosity)
        in zip(reference, candidate)
    ]

    axis.axhline(0.0, color="0.35", linewidth=0.8)
    axis.plot(x_values, residuals, color="tab:green", linewidth=1.2)


if __name__ == "__main__":
    raise SystemExit(main())
