#!/usr/bin/env python3
"""Compare ADAF dynamics radial profiles against a reference solution."""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


COLUMNS = [
    "radius",
    "mach",
    "log_te",
    "log_ti",
    "q_advi_over_q_vis",
    "q_adv_total_over_q_vis",
    "cs_over_c",
    "height",
    "q_rad",
    "tau",
    "log_lk",
    "log_l",
    "magnetic_field",
    "log_rho",
    "log_ne",
]

TOLERANCES = {
    "radius": (1e-10, 1e-8),
    "q_rad": (1e-5, 1e-20),
    "tau": (1e-5, 1e-20),
}

DEFAULT_RTOL = 1e-6
DEFAULT_ATOL = 1e-10
PROFILE_ATOL = 1e-8


@dataclass
class Mismatch:
    column: str
    radius: float
    expected: float
    actual: float
    abs_error: float
    rel_error: float
    rtol: float
    atol: float


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare an ADAF dynamics profile with a fiducial profile."
    )
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    args = parser.parse_args()

    reference = read_profile(args.reference)
    candidate = read_profile(args.candidate)

    failures = compare_profiles(reference, candidate)
    if not failures:
        print(
            f"Dynamics profile matches {args.reference} "
            f"({len(reference)} shells, {len(COLUMNS)} columns)."
        )
        return 0

    print("Dynamics profile differs from fiducial solution:", file=sys.stderr)
    for failure in failures:
        print(
            f"  {failure.column}: radius={failure.radius:.16g} "
            f"expected={failure.expected:.16g} actual={failure.actual:.16g} "
            f"abs_error={failure.abs_error:.3e} rel_error={failure.rel_error:.3e} "
            f"tolerance=atol {failure.atol:.1e} + rtol {failure.rtol:.1e}",
            file=sys.stderr,
        )
    return 1


def read_profile(path: Path) -> List[List[float]]:
    rows: List[List[float]] = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            fields = stripped.replace("D", "E").replace("d", "e").split()
            if len(fields) < len(COLUMNS):
                continue

            try:
                row = [float(value) for value in fields[: len(COLUMNS)]]
            except ValueError:
                continue

            if all(math.isfinite(value) for value in row):
                rows.append(row)

    if not rows:
        raise SystemExit(f"{path} does not contain any numeric dynamics rows.")

    return rows


def compare_profiles(
    reference: List[List[float]], candidate: List[List[float]]
) -> List[Mismatch]:
    if len(reference) != len(candidate):
        raise SystemExit(
            "Dynamics profile shell count changed: "
            f"expected {len(reference)}, got {len(candidate)}."
        )

    failures: List[Mismatch] = []

    for column_index, column in enumerate(COLUMNS):
        rtol, atol = tolerance_for(column)
        worst: Optional[Mismatch] = None

        for expected_row, actual_row in zip(reference, candidate):
            expected = expected_row[column_index]
            actual = actual_row[column_index]
            abs_error = abs(actual - expected)
            rel_error = abs_error / max(abs(expected), PROFILE_ATOL)

            if abs_error <= atol + rtol * abs(expected):
                continue

            mismatch = Mismatch(
                column=column,
                radius=expected_row[0],
                expected=expected,
                actual=actual,
                abs_error=abs_error,
                rel_error=rel_error,
                rtol=rtol,
                atol=atol,
            )

            if worst is None or mismatch.abs_error > worst.abs_error:
                worst = mismatch

        if worst is not None:
            failures.append(worst)

    return failures


def tolerance_for(column: str) -> Tuple[float, float]:
    return TOLERANCES.get(column, (DEFAULT_RTOL, DEFAULT_ATOL))


if __name__ == "__main__":
    raise SystemExit(main())
