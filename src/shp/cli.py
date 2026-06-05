"""Command line interface for SHP."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from .core import DEFAULT_THETA0, scan_fasta


FIELDNAMES = [
    "record_id",
    "length",
    "n_windows",
    "n_transitions",
    "mean_h",
    "mean_d",
    "fixed_wit",
    "tail_energy",
    "skew_d",
    "kurt_d",
]


def _add_scan_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    scan = subparsers.add_parser("scan", help="scan a FASTA file with SHP")
    scan.add_argument("--fasta", required=True, help="input FASTA path, plain text or .gz")
    scan.add_argument("--out", help="output TSV path; defaults to stdout")
    scan.add_argument("--ngram", type=int, default=3, help="n-gram size, default: 3")
    scan.add_argument("--dim", type=int, default=64, help="hash projection dimension, default: 64")
    scan.add_argument("--window", type=int, default=128, help="window size, default: 128")
    scan.add_argument("--stride", type=int, default=None, help="window stride, default: window // 5")
    scan.add_argument(
        "--theta0",
        type=float,
        default=DEFAULT_THETA0,
        help=f"fixed event threshold, default: {DEFAULT_THETA0}",
    )
    scan.add_argument("--alphabet", default="ACGT", help="symbols to keep, default: ACGT")
    scan.set_defaults(func=cmd_scan)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="shp",
        description="SHP: calibrated structural readout for symbolic sequences.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_scan_parser(subparsers)
    return parser


def _format_float(value: float) -> str:
    return f"{value:.6g}"


def cmd_scan(args: argparse.Namespace) -> int:
    out_handle = None
    try:
        if args.out:
            out_path = Path(args.out)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_handle = out_path.open("w", newline="", encoding="utf-8")
            handle = out_handle
        else:
            handle = sys.stdout

        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record_id, result in scan_fasta(
            args.fasta,
            ngram=args.ngram,
            dim=args.dim,
            window=args.window,
            stride=args.stride,
            theta0=args.theta0,
            alphabet=args.alphabet,
        ):
            writer.writerow(
                {
                    "record_id": record_id,
                    "length": result.length,
                    "n_windows": result.n_windows,
                    "n_transitions": result.n_transitions,
                    "mean_h": _format_float(result.mean_h),
                    "mean_d": _format_float(result.mean_d),
                    "fixed_wit": _format_float(result.fixed_wit),
                    "tail_energy": _format_float(result.tail_energy),
                    "skew_d": _format_float(result.skew_d),
                    "kurt_d": _format_float(result.kurt_d),
                }
            )
    finally:
        if out_handle is not None:
            out_handle.close()
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
