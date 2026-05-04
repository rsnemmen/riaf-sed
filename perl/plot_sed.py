#!/usr/bin/env python3
"""Plot a single ADAF spectral energy distribution."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import List, Tuple


def main() -> int:
    parser = argparse.ArgumentParser(description="Create an ADAF SED PNG.")
    parser.add_argument("spectrum", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(args.output.parent / ".matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(args.output.parent / ".cache"))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = read_spectrum(args.spectrum)
    x_values = [row[0] for row in rows]
    y_values = [row[1] for row in rows]

    fig, axis = plt.subplots(figsize=(8, 5), constrained_layout=True)
    axis.plot(x_values, y_values, color="tab:blue", linewidth=1.7)
    axis.set_xlabel(r"$\log_{10}(\nu / \mathrm{Hz})$")
    axis.set_ylabel(r"$\log_{10}(\nu L_\nu / \mathrm{erg}\ \mathrm{s}^{-1})$")
    axis.set_title(f"ADAF SED: {args.spectrum.name}")
    axis.grid(True, alpha=0.25)

    fig.savefig(args.output, dpi=160)
    plt.close(fig)

    print(f"Wrote SED plot to {args.output}")
    return 0


def read_spectrum(path: Path) -> List[Tuple[float, float]]:
    rows: List[Tuple[float, float]] = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            fields = stripped.replace("D", "E").replace("d", "e").split()
            if len(fields) < 2:
                continue

            try:
                log_nu = float(fields[0])
                log_nu_lnu = float(fields[1])
            except ValueError:
                continue

            if not (math.isfinite(log_nu) and math.isfinite(log_nu_lnu)):
                raise SystemExit(f"{path} contains non-finite spectrum values.")

            rows.append((log_nu, log_nu_lnu))

    if not rows:
        raise SystemExit(f"{path} does not contain any numeric spectrum rows.")

    previous_log_nu = None
    for log_nu, _ in rows:
        if previous_log_nu is not None and log_nu <= previous_log_nu:
            raise SystemExit(f"{path} has a non-increasing frequency grid.")
        previous_log_nu = log_nu

    return rows


if __name__ == "__main__":
    raise SystemExit(main())
