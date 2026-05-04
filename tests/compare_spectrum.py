#!/usr/bin/env python3
"""Compare ADAF spectra against a fiducial log-space spectrum."""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


FREQUENCY_ATOL = 1e-8
LUMINOSITY_ATOL = 5e-5
EXPECTED_LARGE_R_ROWS = 100


@dataclass
class SpectrumMismatch:
    column: str
    log_nu: float
    expected: float
    actual: float
    abs_error: float
    atol: float


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare an ADAF spectrum with a fiducial spectrum."
    )
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    args = parser.parse_args()

    reference = read_spectrum(args.reference)
    candidate = read_spectrum(args.candidate)

    failures = compare_spectra(reference, candidate)
    if not failures:
        print(
            f"Spectrum matches {args.reference} "
            f"({len(reference)} frequency samples)."
        )
        return 0

    print("Spectrum differs from fiducial solution:", file=sys.stderr)
    for failure in failures:
        print(
            f"  {failure.column}: log10(nu)={failure.log_nu:.16g} "
            f"expected={failure.expected:.16g} actual={failure.actual:.16g} "
            f"abs_error={failure.abs_error:.3e} tolerance=atol {failure.atol:.1e}",
            file=sys.stderr,
        )
    return 1


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

    previous_log_nu: Optional[float] = None
    for log_nu, _ in rows:
        if previous_log_nu is not None and log_nu <= previous_log_nu:
            raise SystemExit(f"{path} has a non-increasing frequency grid.")
        previous_log_nu = log_nu

    return rows


def compare_spectra(
    reference: List[Tuple[float, float]],
    candidate: List[Tuple[float, float]],
) -> List[SpectrumMismatch]:
    if len(reference) != EXPECTED_LARGE_R_ROWS:
        raise SystemExit(
            "Fiducial largeR spectrum row count changed: "
            f"expected {EXPECTED_LARGE_R_ROWS}, got {len(reference)}."
        )

    if len(reference) != len(candidate):
        raise SystemExit(
            "Spectrum frequency sample count changed: "
            f"expected {len(reference)}, got {len(candidate)}."
        )

    worst_frequency: Optional[SpectrumMismatch] = None
    worst_luminosity: Optional[SpectrumMismatch] = None

    for expected_row, actual_row in zip(reference, candidate):
        expected_log_nu, expected_log_nu_lnu = expected_row
        actual_log_nu, actual_log_nu_lnu = actual_row

        frequency_error = abs(actual_log_nu - expected_log_nu)
        if frequency_error > FREQUENCY_ATOL:
            mismatch = SpectrumMismatch(
                column="log10(nu)",
                log_nu=expected_log_nu,
                expected=expected_log_nu,
                actual=actual_log_nu,
                abs_error=frequency_error,
                atol=FREQUENCY_ATOL,
            )
            if (
                worst_frequency is None
                or mismatch.abs_error > worst_frequency.abs_error
            ):
                worst_frequency = mismatch

        luminosity_error = abs(actual_log_nu_lnu - expected_log_nu_lnu)
        if luminosity_error > LUMINOSITY_ATOL:
            mismatch = SpectrumMismatch(
                column="log10(nu Lnu)",
                log_nu=expected_log_nu,
                expected=expected_log_nu_lnu,
                actual=actual_log_nu_lnu,
                abs_error=luminosity_error,
                atol=LUMINOSITY_ATOL,
            )
            if (
                worst_luminosity is None
                or mismatch.abs_error > worst_luminosity.abs_error
            ):
                worst_luminosity = mismatch

    return [
        failure for failure in (worst_frequency, worst_luminosity) if failure
    ]


if __name__ == "__main__":
    raise SystemExit(main())
