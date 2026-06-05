"""Command line interface for SHP."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from .core import (
    DEFAULT_THETA0,
    calibrate_theta0,
    compute_shp,
    dinuc_shuffle,
    scan_fasta,
)


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
    scan.add_argument("--window", type=int, default=128, help="window size in symbols, default: 128")
    scan.add_argument("--stride", type=int, default=None, help="window stride, default: window // 5")
    scan.add_argument(
        "--theta0",
        type=float,
        default=DEFAULT_THETA0,
        help=f"fixed event threshold, default: {DEFAULT_THETA0}",
    )
    scan.add_argument("--alphabet", default="ACGT", help="symbols to keep, default: ACGT")
    scan.add_argument("--workers", type=int, default=1, help="parallel workers, default: 1")
    scan.add_argument(
        "--shuffle",
        type=int,
        default=0,
        metavar="N",
        help="also scan N dinucleotide-preserving shuffled copies and report null mean",
    )
    scan.set_defaults(func=cmd_scan)


def _add_calibrate_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    cal = subparsers.add_parser("calibrate", help="calibrate theta0 from fair IID streams")
    cal.add_argument("--ngram", type=int, default=3, help="n-gram size, default: 3")
    cal.add_argument("--dim", type=int, default=64, help="hash projection dimension, default: 64")
    cal.add_argument("--window", type=int, default=128, help="window size, default: 128")
    cal.add_argument("--seeds", type=int, default=20, help="number of fair IID seeds, default: 20")
    cal.add_argument(
        "--stream-length",
        type=int,
        default=10000,
        dest="stream_length",
        help="symbols per fair IID stream, default: 10000",
    )
    cal.set_defaults(func=cmd_calibrate)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="shp",
        description="SHP: calibrated structural readout for symbolic sequences.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_scan_parser(subparsers)
    _add_calibrate_parser(subparsers)
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

        shuffle_columns = []
        if args.shuffle > 0:
            for i in range(args.shuffle):
                shuffle_columns.extend([f"shuffle_{i}_fw", f"shuffle_{i}_te"])
        fieldnames = FIELDNAMES + shuffle_columns

        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        # Use scan_fasta for normal mode; manual iteration for shuffle mode
        if args.shuffle > 0:
            from .core import read_fasta
            import random as _random

            for name, seq in read_fasta(args.fasta):
                result = compute_shp(
                    seq,
                    ngram=args.ngram,
                    dim=args.dim,
                    window=args.window,
                    stride=args.stride,
                    theta0=args.theta0,
                    alphabet=args.alphabet,
                )
                if result.length < args.window:
                    print(f"[shp] warning: {name} length {result.length} < window {args.window}", file=sys.stderr)

                row = _make_row(name, result)
                for i in range(args.shuffle):
                    shuf_seq = dinuc_shuffle(seq, seed=hash(f"{name}:{i}") & 0x7FFFFFFF)
                    shuf_res = compute_shp(
                        shuf_seq,
                        ngram=args.ngram,
                        dim=args.dim,
                        window=args.window,
                        stride=args.stride,
                        theta0=args.theta0,
                        alphabet=args.alphabet,
                    )
                    n_null = args.shuffle
                    row[f"shuffle_{i}_fw"] = _format_float(shuf_res.fixed_wit)
                    row[f"shuffle_{i}_te"] = _format_float(shuf_res.tail_energy)
                writer.writerow(row)
        else:
            for record_id, result in scan_fasta(
                args.fasta,
                ngram=args.ngram,
                dim=args.dim,
                window=args.window,
                stride=args.stride,
                theta0=args.theta0,
                alphabet=args.alphabet,
                workers=args.workers,
            ):
                if result.length < args.window:
                    print(f"[shp] warning: {record_id} length {result.length} < window {args.window}", file=sys.stderr)
                writer.writerow(_make_row(record_id, result))

    finally:
        if out_handle is not None:
            out_handle.close()
    return 0


def _make_row(record_id: str, result: object) -> dict:
    return {
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


def cmd_calibrate(args: argparse.Namespace) -> int:
    theta = calibrate_theta0(
        ngram=args.ngram,
        dim=args.dim,
        window=args.window,
        seeds=args.seeds,
        t=args.stream_length,
    )
    print(f"{theta:.6f}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
