"""Generate publication-style figures for the GeneGrammar/SHP preprint.

The script is intentionally self-contained. It reads the frozen Stage 1 matrix
and local Ensembl CDS FASTA when available, then writes PNG, PDF, and SVG
versions of each figure to paper/figures.
"""

from __future__ import annotations

import gzip
import hashlib
import math
import random
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches
from matplotlib.colors import LinearSegmentedColormap


GENEGRAMMAR = Path(__file__).resolve().parents[2]
EE_ROOT = GENEGRAMMAR.parent
PAPER = GENEGRAMMAR / "paper"
FIG_DIR = PAPER / "figures"
STAGE1 = GENEGRAMMAR / "results" / "stage1"
GENE_MATRIX = STAGE1 / "gene_matrix.csv"
CDS_FASTA = EE_ROOT / "revolution" / "ensembl_release_115" / "cds" / "Homo_sapiens.GRCh38.cds.all.fa.gz"

THETA0 = 0.0999
QUIET = 0.01

COL = {
    "ink": "#19212e",
    "muted": "#677383",
    "grid": "#d9dee6",
    "paper": "#fbfbf8",
    "cds": "#245f73",
    "utr": "#d77a3d",
    "chroma": "#1b9aaa",
    "rhythm": "#d45d79",
    "gold": "#cda349",
    "green": "#4f8f5b",
    "purple": "#6f5aa7",
    "red": "#b94a48",
    "blue": "#3f6fa8",
}


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 360,
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "axes.linewidth": 0.8,
            "axes.edgecolor": COL["ink"],
            "axes.facecolor": "white",
            "xtick.color": COL["ink"],
            "ytick.color": COL["ink"],
            "text.color": COL["ink"],
            "grid.color": COL["grid"],
            "grid.linewidth": 0.6,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def save_all(fig: mpl.figure.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"{stem}.{ext}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def label_panel(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.06,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        fontweight="bold",
        color=COL["ink"],
    )


def load_matrix() -> pd.DataFrame:
    df = pd.read_csv(GENE_MATRIX)
    for c in df.columns:
        if c != "gene":
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def stable_hash(text: str, d: int = 64) -> int:
    return int(hashlib.md5(text.encode("utf-8")).hexdigest(), 16) % d


def shp_trace(seq: str, w: int = 128, n: int = 3, d: int = 64) -> tuple[np.ndarray, np.ndarray]:
    stride = max(1, w // 5)
    hs: list[float] = []
    for start in range(0, max(0, len(seq) - w + 1), stride):
        win = seq[start : start + w]
        chroma: set[int] = set()
        rhythm: set[int] = set()
        kmers = [win[i : i + n] for i in range(0, len(win) - n + 1)]
        for k in kmers:
            chroma.add(stable_hash("C:" + k, d))
        for a, b in zip(kmers[:-1], kmers[1:]):
            rhythm.add(stable_hash("R:" + a + ">" + b, d))
        union = chroma | rhythm
        inter = chroma & rhythm
        hs.append(1.0 - (len(inter) / len(union) if union else 0.0))
    h = np.asarray(hs, dtype=float)
    if len(h) < 2:
        return h, np.asarray([], dtype=float)
    return h, np.abs(np.diff(h))


def random_dna(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def block_dna(length: int, block: int = 64) -> str:
    parts: list[str] = []
    symbols = "ACGT"
    while len("".join(parts)) < length:
        for s in symbols:
            parts.append(s * block)
            if len("".join(parts)) >= length:
                break
    return "".join(parts)[:length]


def markov_dna(length: int, seed: int, stay: float = 0.9) -> str:
    rng = random.Random(seed)
    bases = "ACGT"
    cur = rng.choice(bases)
    out = [cur]
    for _ in range(length - 1):
        if rng.random() > stay:
            cur = rng.choice([b for b in bases if b != cur])
        out.append(cur)
    return "".join(out)


def fixed_wit_for(seq: str) -> tuple[float, float]:
    _, d = shp_trace(seq)
    if len(d) == 0:
        return 0.0, 0.0
    fw = float(np.mean(d > THETA0))
    te = float(np.mean(np.maximum(0.0, d - THETA0)))
    return fw, te


def read_gc_by_gene() -> dict[str, float]:
    gc: dict[str, float] = {}
    length: dict[str, int] = {}
    if not CDS_FASTA.exists():
        return gc
    gene = None
    seq_parts: list[str] = []

    def flush() -> None:
        nonlocal gene, seq_parts
        if not gene or not seq_parts:
            return
        seq = "".join(seq_parts).upper()
        if len(seq) > length.get(gene, 0):
            length[gene] = len(seq)
            gc[gene] = (seq.count("G") + seq.count("C")) / max(1, len(seq))

    with gzip.open(CDS_FASTA, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                flush()
                m = re.search(r"gene_symbol:(\S+)", line)
                gene = m.group(1) if m else None
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        flush()
    return gc


REGIME_ORDER = ["HOX", "MHC", "BRAIN", "AD", "KRTAP"]
BRAIN = {
    "FOXP1",
    "FOXP2",
    "PAX6",
    "NEUROD1",
    "NEUROD2",
    "NEUROD4",
    "NEUROD6",
    "TBR1",
    "DLG4",
    "SYN1",
    "SYN2",
    "SYN3",
    "GABRA1",
    "GABRA2",
    "GABRA3",
    "GABRA4",
    "GABRA5",
    "GABRA6",
    "GRIN1",
    "GRIN2A",
    "GRIN2B",
    "GRIN2C",
    "GRIN2D",
    "GRIN3A",
    "GRIN3B",
}
AD = {
    "MAPT",
    "GSK3B",
    "LAMP1",
    "CD33",
    "BCL2",
    "APP",
    "PSEN1",
    "PSEN2",
    "APOE",
    "TREM2",
    "CLU",
    "CR1",
    "BIN1",
    "PICALM",
    "ABCA7",
    "MS4A6A",
    "EPHA1",
    "CD2AP",
    "SORL1",
    "FERMT2",
}


def regime_for(gene: str) -> str | None:
    if gene.startswith("HOX"):
        return "HOX"
    if gene.startswith("HLA-"):
        return "MHC"
    if gene in BRAIN:
        return "BRAIN"
    if gene in AD:
        return "AD"
    if gene.startswith("KRTAP"):
        return "KRTAP"
    return None


CAT_RULES = {
    "ZincFinger": lambda g: g.startswith("ZNF"),
    "SoluteCarrier": lambda g: g.startswith("SLC"),
    "TranscriptionFactor": lambda g: any(
        g.startswith(p) for p in ["HOX", "FOX", "SOX", "PAX", "TBX", "NKX", "LHX", "GATA", "STAT", "SMAD", "NFKB", "TCF"]
    ),
    "OlfactoryReceptor": lambda g: g.startswith("OR") and len(g) <= 6,
    "Keratin": lambda g: g.startswith("KRT"),
    "Histone": lambda g: any(g.startswith(p) for p in ["HIST", "H2A", "H2B", "H3", "H4", "H1-"]),
    "Ribosomal": lambda g: any(g.startswith(p) for p in ["RPL", "RPS", "MRPL", "MRPS"]),
    "GPCR": lambda g: any(g.startswith(p) for p in ["GPR", "ADGR", "GPRC", "ADRA", "ADRB"]),
    "HLA": lambda g: g.startswith("HLA-"),
}


def category_for(gene: str) -> str | None:
    for cat, fn in CAT_RULES.items():
        if fn(gene):
            return cat
    return None


def draw_bitset(ax: mpl.axes.Axes, x: float, y: float, active: list[int], color: str, label: str) -> None:
    ax.text(x, y + 0.28, label, ha="left", va="center", fontsize=8, color=COL["ink"], fontweight="bold")
    for i in range(16):
        xx = x + (i % 8) * 0.075
        yy = y - (i // 8) * 0.075
        fc = color if i in active else "#eef1f4"
        ec = color if i in active else "#d9dee6"
        ax.add_patch(patches.Rectangle((xx, yy), 0.052, 0.052, fc=fc, ec=ec, lw=0.5))


def fig1_method() -> None:
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    ax.set_axis_off()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)

    def box(x: float, y: float, w: float, h: float, text: str, fc: str, ec: str = "none", fs: int = 9) -> None:
        ax.add_patch(
            patches.FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.02,rounding_size=0.08",
                fc=fc,
                ec=ec if ec != "none" else fc,
                lw=1.0,
            )
        )
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=COL["ink"])

    def arrow(x1: float, y1: float, x2: float, y2: float, color: str = COL["muted"]) -> None:
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", lw=1.4, color=color))

    def readout_matrix(x: float, y: float, w: float, h: float) -> None:
        ax.add_patch(
            patches.FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.02,rounding_size=0.08",
                fc="#edf4ed",
                ec=COL["green"],
                lw=1.0,
            )
        )
        ax.text(x + w / 2, y + h - 0.18, "readout matrix", ha="center", va="center", fontsize=8, fontweight="bold")
        left = x + 0.25
        col1 = x + w - 1.35
        col2 = x + w - 0.55
        top = y + h - 0.42
        ax.text(col1, top, "CDS", ha="center", va="center", fontsize=6.8, color=COL["muted"], fontweight="bold")
        ax.text(col2, top, "UTR", ha="center", va="center", fontsize=6.8, color=COL["muted"], fontweight="bold")
        rows = [
            ("fixed_wit", "0.012", "0.030"),
            ("tail_energy", "0.0001", "0.0007"),
            ("skew_d", "+0.8", "+1.2"),
            ("kurt_d", "+0.1", "+1.7"),
        ]
        row_gap = 0.17
        for i, (name, cds, utr) in enumerate(rows):
            yy = top - (i + 1) * row_gap
            ax.text(left, yy, name, ha="left", va="center", fontsize=6.7)
            ax.text(col1, yy, cds, ha="center", va="center", fontsize=6.7)
            ax.text(col2, yy, utr, ha="center", va="center", fontsize=6.7)
            if i < len(rows) - 1:
                ax.plot([x + 0.18, x + w - 0.18], [yy - 0.105, yy - 0.105], color="#d7e3d8", lw=0.6)

    def axis_panel(
        x: float,
        y: float,
        w: float,
        h: float,
        title: str,
        subtitle: str,
        label: str,
        active: list[int],
        fc: str,
        color: str,
    ) -> None:
        ax.add_patch(
            patches.FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.02,rounding_size=0.08",
                fc=fc,
                ec=color,
                lw=1.0,
            )
        )
        ax.text(x + w / 2, y + h - 0.20, title, ha="center", va="center", fontsize=9, color=COL["ink"])
        ax.text(x + w / 2, y + h - 0.42, subtitle, ha="center", va="center", fontsize=7.2, color=COL["ink"])
        ax.text(x + 0.43, y + 0.23, label, ha="left", va="center", fontsize=7.2, color=COL["ink"], fontweight="bold")
        sx = x + 0.95
        sy = y + 0.25
        for i in range(16):
            xx = sx + (i % 8) * 0.075
            yy = sy - (i // 8) * 0.075
            fill = color if i in active else "#eef1f4"
            edge = color if i in active else "#d9dee6"
            ax.add_patch(patches.Rectangle((xx, yy), 0.052, 0.052, fc=fill, ec=edge, lw=0.5))

    ax.text(0.25, 5.65, "Figure 1. SHP turns one nucleotide stream into two structural views", fontsize=14, fontweight="bold")
    ax.text(
        0.25,
        5.35,
        "Chroma asks what 3-mers are present. Rhythm asks how 3-mers transition. Cross-harm measures their mismatch.",
        fontsize=9,
        color=COL["muted"],
    )

    seq = "ATG CCG TTA GGC AAG TCC GAT TGA".split()
    x0 = 0.55
    for i, kmer in enumerate(seq):
        fc = ["#e7f0f2", "#f7eadf", "#e9edf7"][i % 3]
        box(x0 + i * 0.64, 4.55, 0.52, 0.36, kmer, fc, fs=8)
    box(0.45, 4.05, 5.2, 0.28, "128 nt sliding window", "#f3f5f7", "#d9dee6", fs=8)

    axis_panel(0.45, 2.72, 2.45, 0.83, "Chroma", "3-mer presence", "C_t", [0, 2, 5, 7, 10, 13], "#dff3f5", COL["chroma"])
    axis_panel(3.20, 2.72, 2.45, 0.83, "Rhythm", "3-mer transitions", "R_t", [1, 2, 6, 8, 10, 14], "#f9e4ea", COL["rhythm"])
    arrow(1.68, 4.05, 1.68, 3.55, COL["chroma"])
    arrow(4.43, 4.05, 4.43, 3.55, COL["rhythm"])
    arrow(1.68, 2.72, 2.65, 2.20, COL["chroma"])
    arrow(4.43, 2.72, 4.05, 2.20, COL["rhythm"])
    box(
        1.45,
        1.42,
        3.8,
        0.78,
        "cross-harm trace\nh_t = 1 - J(C_t, R_t)\nd_t > theta_0 -> event",
        "#f4f1e4",
        COL["gold"],
        fs=7,
    )

    readout_matrix(5.85, 1.23, 3.35, 1.18)
    arrow(5.25, 1.81, 5.85, 1.81, COL["green"])
    save_all(fig, "fig1_shp_method")


def fig2_calibration_and_genome(df: pd.DataFrame) -> None:
    rng_ds: list[float] = []
    for seed in range(20):
        _, d = shp_trace(random_dna(10000, seed))
        rng_ds.extend(d.tolist())
    fair_d = np.asarray(rng_ds)

    streams = {
        "Periodic": "ACGT" * 2500,
        "Fair IID": random_dna(10000, 100),
        "Markov 0.9": markov_dna(10000, 101, 0.9),
        "Block 64": block_dna(10000, 64),
    }
    fw = {name: fixed_wit_for(seq)[0] for name, seq in streams.items()}
    te = {name: fixed_wit_for(seq)[1] for name, seq in streams.items()}

    cds = df[df["n_cds_windows"] >= 10]["cds_fw"].dropna()
    utr = df[df["n_utr_windows"] >= 10]["utr_fw"].dropna()

    fig = plt.figure(figsize=(10.6, 7.2))
    gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.28)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    label_panel(ax1, "A")
    ax1.hist(fair_d, bins=44, color="#dfe8ed", edgecolor="white")
    ax1.axvline(THETA0, color=COL["red"], lw=2)
    ax1.text(THETA0 + 0.004, ax1.get_ylim()[1] * 0.86, "theta0 = 0.0999", color=COL["red"], fontsize=9)
    ax1.set_title("Fair-IID calibration")
    ax1.set_xlabel("cross-harm displacement d_t")
    ax1.set_ylabel("windows")
    ax1.grid(axis="y", alpha=0.5)

    label_panel(ax2, "B")
    names = list(fw)
    vals = [fw[n] for n in names]
    colors = ["#b9c2cf", COL["blue"], COL["purple"], COL["gold"]]
    ax2.bar(names, vals, color=colors, width=0.68)
    ax2.set_title("Fixed-threshold control streams")
    ax2.set_ylabel("fixed_wit")
    ax2.tick_params(axis="x", rotation=18)
    ymax = max(vals) * 1.32 if vals else 0.06
    ax2.set_ylim(0, ymax)
    for i, n in enumerate(names):
        ax2.text(i, vals[i] + ymax * 0.04, f"TE={te[n]:.4f}", ha="center", fontsize=7, color=COL["muted"])
    ax2.grid(axis="y", alpha=0.5)

    label_panel(ax3, "C")
    bins = np.linspace(0, 0.11, 44)
    ax3.hist(cds, bins=bins, alpha=0.78, color=COL["cds"], label=f"CDS (N={len(cds):,})")
    ax3.hist(utr, bins=bins, alpha=0.58, color=COL["utr"], label=f"UTR (N={len(utr):,})")
    ax3.axvline(QUIET, color=COL["ink"], lw=1.1, ls="--")
    ax3.set_title("Genome-wide SHP event-rate distribution")
    ax3.set_xlabel("fixed_wit")
    ax3.set_ylabel("genes")
    ax3.legend()
    ax3.grid(axis="y", alpha=0.5)

    label_panel(ax4, "D")
    exact = [float(np.mean(cds == 0)), float(np.mean(utr == 0))]
    quiet = [float(np.mean(cds < QUIET)), float(np.mean(utr < QUIET))]
    x = np.arange(2)
    ax4.bar(x - 0.17, exact, width=0.34, color=["#9ab5c0", "#e1ab83"], label="exact zero")
    ax4.bar(x + 0.17, quiet, width=0.34, color=[COL["cds"], COL["utr"]], label="fw < 0.01")
    ax4.set_xticks(x, ["CDS", "UTR"])
    ax4.set_ylim(0, 0.58)
    ax4.set_ylabel("fraction of genes")
    ax4.set_title("Quiescent sequence regions")
    ax4.legend(loc="upper right")
    for xpos, val in zip(x - 0.17, exact):
        ax4.text(xpos, val + 0.015, f"{100*val:.1f}%", ha="center", fontsize=8)
    for xpos, val in zip(x + 0.17, quiet):
        ax4.text(xpos, val + 0.015, f"{100*val:.1f}%", ha="center", fontsize=8)
    ax4.grid(axis="y", alpha=0.5)

    fig.suptitle("Figure 2. Calibration and genome-wide structural event rates", x=0.02, ha="left", fontsize=14, fontweight="bold")
    save_all(fig, "fig2_calibration_genome")


def regime_summary(df: pd.DataFrame) -> pd.DataFrame:
    temp = df.copy()
    temp["regime"] = temp["gene"].map(regime_for)
    temp = temp[temp["regime"].notna()].copy()
    rows = []
    for reg in REGIME_ORDER:
        sub = temp[temp["regime"] == reg]
        cds = sub[sub["n_cds_windows"] >= 10]["cds_fw"].dropna()
        utr = sub[sub["n_utr_windows"] >= 10]["utr_fw"].dropna()
        rows.append(
            {
                "regime": reg,
                "N": len(sub),
                "cds_fw": float(cds.mean()) if len(cds) else np.nan,
                "utr_fw": float(utr.mean()) if len(utr) else np.nan,
                "cds_quiet": float(np.mean(cds < QUIET)) if len(cds) else np.nan,
                "utr_quiet": float(np.mean(utr < QUIET)) if len(utr) else np.nan,
                "gradient": (float(utr.mean()) if len(utr) else np.nan) - (float(cds.mean()) if len(cds) else np.nan),
            }
        )
    return pd.DataFrame(rows)


def fig3_regime_matrix(df: pd.DataFrame) -> None:
    rs = regime_summary(df)
    fig = plt.figure(figsize=(10.8, 6.9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=0.36, wspace=0.26)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    label_panel(ax1, "A")
    x = np.arange(len(rs))
    ax1.bar(x - 0.18, rs["cds_fw"], width=0.36, color=COL["cds"], label="CDS")
    ax1.bar(x + 0.18, rs["utr_fw"], width=0.36, color=COL["utr"], label="UTR")
    ax1.axhline(QUIET, color=COL["ink"], lw=1, ls="--", alpha=0.75)
    ax1.set_xticks(x, rs["regime"])
    ax1.set_ylabel("mean fixed_wit")
    ax1.set_title("CDS/UTR structural gradient by regime")
    ax1.legend()
    ax1.grid(axis="y", alpha=0.5)

    label_panel(ax2, "B")
    size = np.clip(rs["N"].to_numpy() * 9, 120, 700)
    colors = [COL["green"], COL["red"], COL["blue"], COL["purple"], COL["gold"]]
    ax2.scatter(rs["cds_fw"], rs["utr_fw"], s=size, c=colors, alpha=0.9, edgecolor="white", linewidth=1.5)
    lim = [0, max(rs["cds_fw"].max(), rs["utr_fw"].max()) * 1.22]
    ax2.plot(lim, lim, color=COL["grid"], lw=1.2, ls="--")
    text_pos = {
        "HOX": (0.0067, 0.0206),
        "AD": (0.0130, 0.0215),
        "MHC": (0.0118, 0.0332),
        "BRAIN": (0.0222, 0.0297),
        "KRTAP": (0.0363, 0.0128),
    }
    for _, r in rs.iterrows():
        tx, ty = text_pos.get(r["regime"], (r["cds_fw"] + 0.0008, r["utr_fw"] + 0.0008))
        ax2.annotate(
            r["regime"],
            xy=(r["cds_fw"], r["utr_fw"]),
            xytext=(tx, ty),
            fontsize=9,
            fontweight="bold",
            arrowprops=dict(arrowstyle="-", lw=0.7, color=COL["muted"], alpha=0.65) if r["regime"] in {"HOX", "AD"} else None,
        )
    ax2.set_xlim(lim)
    ax2.set_ylim(lim)
    ax2.set_xlabel("CDS fixed_wit")
    ax2.set_ylabel("UTR fixed_wit")
    ax2.set_title("Only KRTAP is CDS-led")
    ax2.grid(alpha=0.45)

    label_panel(ax3, "C")
    heat = rs[["cds_fw", "utr_fw", "cds_quiet", "utr_quiet", "gradient"]].to_numpy()
    labels = ["CDS fw", "UTR fw", "CDS quiet", "UTR quiet", "UTR-CDS"]
    normed = heat.copy()
    for j in range(normed.shape[1]):
        col = normed[:, j]
        mn, mx = np.nanmin(col), np.nanmax(col)
        normed[:, j] = (col - mn) / (mx - mn + 1e-12)
    cmap = LinearSegmentedColormap.from_list("shp_heat", ["#f7f4ec", "#d7b66d", "#7e496c", "#20364f"])
    ax3.imshow(normed, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax3.set_yticks(np.arange(len(rs)), rs["regime"])
    ax3.set_xticks(np.arange(len(labels)), labels)
    ax3.set_title("Regime signatures in the SHP matrix")
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            v = heat[i, j]
            txt = f"{v:.3f}" if j in (0, 1, 4) else f"{100*v:.0f}%"
            ax3.text(j, i, txt, ha="center", va="center", fontsize=8, color="white" if normed[i, j] > 0.56 else COL["ink"])
    for spine in ax3.spines.values():
        spine.set_visible(False)

    fig.suptitle("Figure 3. Biological regimes occupy distinct CDS/UTR structural states", x=0.02, ha="left", fontsize=14, fontweight="bold")
    save_all(fig, "fig3_regime_signatures")


def functional_summary(df: pd.DataFrame) -> pd.DataFrame:
    gc = read_gc_by_gene()
    temp = df[df["n_cds_windows"] >= 10].copy()
    temp["cat"] = temp["gene"].map(category_for)
    temp = temp[temp["cat"].notna()].copy()
    temp["gc"] = temp["gene"].map(gc)
    rows = []
    for cat in CAT_RULES:
        sub = temp[temp["cat"] == cat]
        if len(sub) < 10:
            continue
        rows.append(
            {
                "cat": cat,
                "N": len(sub),
                "cds_fw": float(sub["cds_fw"].mean()),
                "utr_fw": float(sub["utr_fw"].fillna(0).mean()),
                "gc": float(sub["gc"].mean()) if sub["gc"].notna().any() else np.nan,
            }
        )
    return pd.DataFrame(rows)


def nearest_centroid_accuracy(df: pd.DataFrame) -> pd.DataFrame:
    temp = df[df["n_cds_windows"] >= 10].copy()
    temp["cat"] = temp["gene"].map(category_for)
    temp = temp[temp["cat"].notna()].copy()
    temp["utr_fw"] = temp["utr_fw"].fillna(0)
    cats = [c for c, n in temp["cat"].value_counts().items() if n >= 10]
    temp = temp[temp["cat"].isin(cats)].reset_index(drop=True)
    out = []
    for cat in cats:
        sub = temp[temp["cat"] == cat]
        correct = 0
        for idx, row in sub.iterrows():
            cent = {}
            for c in cats:
                pool = temp[temp["cat"] == c]
                if c == cat:
                    pool = pool[pool.index != idx]
                cent[c] = (pool["cds_fw"].mean(), pool["utr_fw"].mean())
            best = min(cent, key=lambda c: math.hypot(row["cds_fw"] - cent[c][0], row["utr_fw"] - cent[c][1]))
            correct += int(best == cat)
        out.append({"cat": cat, "N": len(sub), "recall": correct / len(sub)})
    return pd.DataFrame(out)


def fig4_functional_orthogonality(df: pd.DataFrame) -> None:
    fs = functional_summary(df)
    acc = nearest_centroid_accuracy(df)
    merged = fs.merge(acc, on=["cat", "N"], how="left")
    merged = merged.sort_values("N", ascending=False)

    fig = plt.figure(figsize=(11.0, 7.0))
    gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.30)
    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])

    label_panel(ax1, "A")
    sc = ax1.scatter(
        merged["cds_fw"],
        merged["utr_fw"],
        s=np.clip(merged["N"] * 1.15, 90, 620),
        c=merged["gc"] * 100,
        cmap="viridis",
        edgecolor="white",
        linewidth=1.4,
        alpha=0.92,
    )
    label_offsets = {
        "HLA": (0.00045, 0.00045),
        "TranscriptionFactor": (0.00045, 0.00045),
        "ZincFinger": (0.00045, 0.00035),
        "GPCR": (-0.00085, 0.00025),
        "SoluteCarrier": (0.00045, -0.00018),
        "Ribosomal": (0.00045, 0.00042),
        "Keratin": (-0.0019, 0.00042),
        "OlfactoryReceptor": (0.00055, 0.00035),
        "Histone": (0.00055, 0.00035),
    }
    for _, r in merged.iterrows():
        dx, dy = label_offsets.get(r["cat"], (0.00045, 0.00045))
        ax1.text(r["cds_fw"] + dx, r["utr_fw"] + dy, r["cat"], fontsize=8)
    ax1.set_xlim(merged["cds_fw"].min() - 0.001, merged["cds_fw"].max() + 0.002)
    ax1.set_ylim(merged["utr_fw"].min() - 0.001, merged["utr_fw"].max() + 0.0014)
    cb = fig.colorbar(sc, ax=ax1, shrink=0.78, pad=0.02)
    cb.set_label("CDS GC%")
    ax1.set_xlabel("CDS fixed_wit")
    ax1.set_ylabel("UTR fixed_wit")
    ax1.set_title("Functional categories in SHP space")
    ax1.grid(alpha=0.45)

    label_panel(ax2, "B")
    pairs = [
        ("ZincFinger", "TranscriptionFactor"),
        ("OlfactoryReceptor", "Histone"),
        ("SoluteCarrier", "GPCR"),
    ]
    y = np.arange(len(pairs))
    shp_d = []
    gc_d = []
    for a, b in pairs:
        ra = fs[fs["cat"] == a].iloc[0]
        rb = fs[fs["cat"] == b].iloc[0]
        shp_d.append(math.hypot((ra["cds_fw"] - rb["cds_fw"]) * 1000, (ra["utr_fw"] - rb["utr_fw"]) * 1000))
        gc_d.append(abs(ra["gc"] - rb["gc"]) * 100)
    ax2.barh(y - 0.18, shp_d, height=0.32, color=COL["blue"], label="SHP distance x1000")
    ax2.barh(y + 0.18, gc_d, height=0.32, color=COL["utr"], label="GC difference (pp)")
    ax2.set_yticks(y, [f"{a}\nvs {b}" for a, b in pairs])
    ax2.invert_yaxis()
    ax2.set_title("SHP is not a GC proxy")
    ax2.legend(fontsize=8)
    ax2.grid(axis="x", alpha=0.45)

    label_panel(ax3, "C")
    acc_plot = merged.sort_values("recall", ascending=True).tail(9)
    y = np.arange(len(acc_plot))
    ax3.barh(y, acc_plot["recall"], color="#9ab5c0")
    ax3.axvline(1 / 9, color=COL["red"], lw=1.4, ls="--", label="chance")
    ax3.set_xlim(0, max(0.48, acc_plot["recall"].max() * 1.18))
    ax3.set_yticks(y, acc_plot["cat"])
    ax3.set_xlabel("LOO recall")
    ax3.set_title("2D nearest-centroid recall")
    ax3.legend()
    for i, (_, r) in enumerate(acc_plot.iterrows()):
        ax3.text(r["recall"] + 0.012, i, f"{100*r['recall']:.0f}%", va="center", fontsize=8)
    ax3.grid(axis="x", alpha=0.45)

    fig.suptitle("Figure 4. SHP carries functional signal orthogonal to simple composition", x=0.02, ha="left", fontsize=14, fontweight="bold")
    save_all(fig, "fig4_functional_orthogonality")


def write_captions() -> None:
    text = """# Figure Captions

## Figure 1. SHP turns one nucleotide stream into two structural views

SHP projects each local nucleotide window into two binary hash activations:
chroma, which records which 3-mers are present, and rhythm, which records which
adjacent 3-mer transitions occur. Cross-harm is the Jaccard distance between
these two views. Consecutive cross-harm displacement defines calibrated
structural events (`fixed_wit`) and excess event intensity (`tail_energy`),
which are assembled into CDS/UTR readout matrices for downstream screening.

## Figure 2. Calibration and genome-wide structural event rates

(A) Fair-IID displacement distribution used to define the fixed threshold
`theta0 = 0.0999`. (B) Fixed-threshold readout separates periodic, fair-IID,
Markov-biased, and block-structured streams. (C) Genome-wide CDS and UTR
`fixed_wit` distributions after the `n_windows >= 10` filter. (D) Exact-zero
and quiescent (`fixed_wit < 0.01`) fractions for CDS and UTR regions.

## Figure 3. Biological regimes occupy distinct CDS/UTR structural states

(A) Mean CDS and UTR `fixed_wit` by regime. (B) Regime centroids in CDS/UTR
structural space; points above the diagonal are UTR-led, points below are
CDS-led. (C) Min-max normalized heat map of regime signatures, with raw values
shown in each cell.

## Figure 4. SHP carries functional signal orthogonal to simple composition

(A) Broad functional gene-symbol categories plotted in CDS/UTR SHP space, with
point size proportional to category size and color indicating CDS GC content.
(B) Selected category pairs show that proximity in SHP space does not reduce to
GC similarity. (C) Leave-one-out nearest-centroid recall using only two SHP
features (`CDS fixed_wit`, `UTR fixed_wit`).
"""
    (PAPER / "figure_captions.md").write_text(text, encoding="utf-8")


def main() -> None:
    setup_style()
    df = load_matrix()
    fig1_method()
    fig2_calibration_and_genome(df)
    fig3_regime_matrix(df)
    fig4_functional_orthogonality(df)
    write_captions()
    print(f"Wrote figures to {FIG_DIR}")


if __name__ == "__main__":
    main()
