#!/usr/bin/env python3
"""
render_hic_tiles.py — Batch‑render C‑Origami / Hi‑C tiles stored as .npy matrices,
optionally with per‑tile tracks from bigWig (e.g., CTCF, ATAC) plotted beneath.

Given a directory of files like chr11_0.npy, chr11_2097152.npy, ..., this script
will iterate, load with numpy (mmap for low RAM), optionally crop to a sub‑window,
optionally log‑transform and clip intensities, optionally downsample, and write a
PNG/PDF/SVG heatmap per matrix. If bigWig tracks are provided, the script extracts
matching genomic windows and renders them as stacked tracks aligned to the heatmap
x‑axis (absolute chromosomal Mb), similar to the C‑Origami figure panels.

Example:
  python render_hic_tiles.py \
    --input /path/to/origami_chr11/npy \
    --output ./hic_png \
    --log1p --clip-percentile 99.5 --downsample 2 --workers 8 \
    --ctcf-bw Ly6CLo_Old_WT4_CTCF.bw \
    --atac-bw Ly6CLo_Old_WT4_IS_slop20_RP20M_minmax01.bw \
    --start 0 --end 500kb --tick-mb 0.25 --scalebar-kb 100

Notes
- origin='lower' so the diagonal runs bottom‑left → top‑right.
- If a matrix is 1‑D and a perfect square, we reshape it; otherwise we error.
- Downsample performs average pooling by an integer factor to keep files small.
- Files are auto‑sorted by chromosome + numeric offset when names follow chr*_*.
- When tracks are supplied, pyBigWig is required (conda/pip: pyBigWig).
"""

import argparse
import logging
import math
from pathlib import Path
import re
from typing import List, Optional, Tuple, Union

import json
import csv
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib
matplotlib.use("Agg")  # headless backend for clusters
import matplotlib.pyplot as plt
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

# ----------------------------- helpers ---------------------------------

CHR_OFFSET_RE = re.compile(r"^(chr[^_]+)_(\d+)\.npy$")

def discover_npy_files(input_dir: Path, pattern: str = "*.npy") -> List[Path]:
    files = sorted(input_dir.glob(pattern))
    if not files:
        logging.warning("No .npy files found under %s", input_dir)
    return files

def parse_chr_offset(p: Path) -> Tuple[Optional[str], int]:
    m = CHR_OFFSET_RE.match(p.name)
    if m:
        return m.group(1), int(m.group(2))
    # fallback: keep original order if pattern doesn't match
    return None, 0

def sort_tiles(files: List[Path]) -> List[Path]:
    try:
        return sorted(files, key=lambda x: (parse_chr_offset(x)[0] or x.name, parse_chr_offset(x)[1]))
    except Exception:
        return sorted(files)

def reshape_if_vector(a: np.ndarray) -> np.ndarray:
    if a.ndim == 2:
        return a
    if a.ndim == 1:
        n = int(math.isqrt(a.size))
        if n * n == a.size:
            return a.reshape(n, n)
        raise ValueError(f"1D array of length {a.size} cannot be reshaped to square")
    raise ValueError(f"Expected 1D or 2D array, got ndim={a.ndim}")

def downsample_mean(a: np.ndarray, factor: int) -> np.ndarray:
    if factor <= 1:
        return a
    h, w = a.shape
    nh, nw = h // factor, w // factor
    a = a[: nh * factor, : nw * factor]
    a = a.reshape(nh, factor, nw, factor).mean(axis=(1, 3))
    return a

def transform_intensity(a: np.ndarray, log1p: bool, clip_pct: Optional[float],
                         vmin: Optional[float], vmax: Optional[float]) -> Tuple[np.ndarray, float, float]:
    # Replace infs safely to a finite maximum so log/percentile behave
    if np.isfinite(a).any():
        finite_max = float(np.nanmax(a[np.isfinite(a)]))
    else:
        finite_max = 0.0
    a = np.nan_to_num(a, copy=False, posinf=finite_max, neginf=np.nanmin(a) if np.isfinite(a).any() else 0.0)

    if log1p:
        a = np.log1p(a)
    # percentile clip overrides vmin/vmax if provided
#    if clip_pct is not None:
#        if not (0.0 < clip_pct <= 100.0):
#            raise ValueError("clip-percentile must be in (0, 100]")
#        hi = float(np.percentile(a[np.isfinite(a)], clip_pct))
#        a = np.clip(a, a_min=None, a_max=hi)
#        vmin_calc, vmax_calc = float(np.nanmin(a)), float(np.nanmax(a))
#        return a, vmin_calc, vmax_calc
    if clip_pct is not None:
            if not (0.0 < clip_pct <= 100.0):
                raise ValueError("clip-percentile must be in (0, 100]")
            # --- symmetric limits around 0 ---
            hi = float(np.percentile(np.abs(a[np.isfinite(a)]), clip_pct))
            a = np.clip(a, -hi, hi)
            return a, -hi, hi

# else honor explicit vmin/vmax or compute from data
    vmin_calc = float(np.nanmin(a)) if vmin is None else vmin
    vmax_calc = float(np.nanmax(a)) if vmax is None else vmax
    return a, vmin_calc, vmax_calc

def add_scalebar_mb(ax: plt.Axes, scalebar_kb: float,
                    frac_from_edge: float = 0.06, color: str = "#CC5500"):
    if not scalebar_kb or scalebar_kb <= 0:
        return
    L_mb = scalebar_kb / 1000.0
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xs = x1 - (x1 - x0) * frac_from_edge - L_mb
    ys = y0 + (y1 - y0) * frac_from_edge
    ax.plot([xs, xs + L_mb], [ys, ys], linewidth=3, color=color, solid_capstyle='butt')
    ax.text(xs + L_mb/2, ys + (y1 - y0) * 0.03, f"{int(scalebar_kb)} kb",
            ha='center', va='bottom', fontsize=9, color=color)

# ----------------------------- stats helpers ----------------------------

def compute_tile_stats(arr: np.ndarray, npy_path: Path, bp_per_bin: float,
                        vmin_final: float, vmax_final: float) -> dict:
    """Compute per‑tile stats (on the *rendered* array) for sidecar/CSV."""
    vals = arr[np.isfinite(arr)]
    stats = {
        "file": npy_path.name,
        "shape": [int(arr.shape[0]), int(arr.shape[1])],
        "dtype": str(arr.dtype),
        "tile_width_bp": int(round(bp_per_bin * arr.shape[0])),
        "bp_per_bin": float(bp_per_bin),
        "vmin": float(vmin_final) if vmin_final is not None else None,
        "vmax": float(vmax_final) if vmax_final is not None else None,
        "min": float(vals.min()) if vals.size else None,
        "max": float(vals.max()) if vals.size else None,
        "median": float(np.median(vals)) if vals.size else None,
        "mean": float(np.mean(vals)) if vals.size else None,
        "neg_fraction": float((vals < 0).mean()) if vals.size else None,
    }
    if vals.size:
        abs_p99 = float(np.percentile(np.abs(vals), 99.0))
        stats.update({
            "abs_p99": abs_p99,
            "recommended_vmin": -abs_p99,
            "recommended_vmax": abs_p99,
        })
    return stats

def write_sidecar_json(out_path: Path, stats: dict) -> Path:
    json_path = out_path.with_suffix('.json')
    with open(json_path, 'w') as fh:
        json.dump(stats, fh, indent=2)
    return json_path

def collect_stats_csv_from_sidecars(out_dir: Path, csv_path: Path) -> Optional[Path]:
    json_paths = sorted(out_dir.glob("*.json"))
    if not json_paths:
        logging.warning("No sidecar JSONs found in %s; CSV not written", out_dir)
        return None
    rows = []
    for jp in json_paths:
        try:
            with open(jp) as fh:
                d = json.load(fh)
            shape = d.get("shape")
            d["shape_str"] = f"{shape[0]}x{shape[1]}" if isinstance(shape, list) and len(shape) == 2 else None
            rows.append(d)
        except Exception as e:
            logging.error("Failed to read %s: %s", jp.name, e)
    if not rows:
        logging.warning("No readable sidecars; CSV not written")
        return None
    fieldnames = [
        "file", "shape_str", "dtype", "tile_width_bp", "bp_per_bin",
        "vmin", "vmax", "min", "max", "median", "mean",
        "neg_fraction", "abs_p99", "recommended_vmin", "recommended_vmax"
    ]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fieldnames})
    logging.info("Wrote stats CSV: %s (%d rows)", csv_path, len(rows))
    return csv_path

# ----------------------------- bigWig helpers ---------------------------

def parse_bp(value: Optional[Union[str, int, float]]) -> Optional[int]:
    """Parse a genomic bp string with optional suffix (kb/mb) into bp."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return int(value)
    s = str(value).strip().replace(",", "").lower()
    mult = 1
    if s.endswith("bp"):
        s = s[:-2]
    elif s.endswith("kb") or s.endswith("k"):
        mult = 1_000
        s = s[:-2] if s.endswith("kb") else s[:-1]
    elif s.endswith("mb") or s.endswith("m"):
        mult = 1_000_000
        s = s[:-2] if s.endswith("mb") else s[:-1]
    try:
        return int(float(s) * mult)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Could not parse genomic value: {value!r}")

def _open_bw(path: Optional[Path]):
    if not path:
        return None
    try:
        import pyBigWig  # local import so base rendering works without it
    except ImportError as e:
        raise RuntimeError("pyBigWig is required when providing bigWig tracks.\n"
                           "Install with: conda install -c bioconda pybigwig  (or pip install pybigwig)") from e
    bw = pyBigWig.open(str(path))
    if bw is None:
        raise RuntimeError(f"Failed to open bigWig: {path}")
    return bw

def bw_summary(bw, chrom: str, start: int, end: int, n_bins: int) -> np.ndarray:
    """Mean signal in equal‑width bins across [start,end). Uses pyBigWig.stats."""
    vals = bw.stats(chrom, int(start), int(end), nBins=int(max(1, n_bins)), type='mean')
    arr = np.array([np.nan if v is None else float(v) for v in vals], dtype=float)
    return arr

def smooth_1d(arr: np.ndarray, win: int) -> np.ndarray:
    if win is None or win <= 1:
        return arr
    k = np.ones(win, dtype=float)
    valid = np.isfinite(arr)
    arr2 = np.copy(arr)
    arr2[~valid] = 0.0 #arr2[\u223cvalid] = 0.0
    num = np.convolve(arr2, k, mode='same')
    den = np.convolve(valid.astype(float), k, mode='same')
    out = np.divide(num, den, out=np.full_like(arr, np.nan, dtype=float), where=den>0)
    return out

# ----------------------------- renderer --------------------------------

def render_tile(
    npy_a: Path, npy_b: Path, out_dir: Path, fmt: str, dpi: int, figsize: float,
    log1p: bool, clip_pct: Optional[float], vmin: Optional[float], vmax: Optional[float],
    cmap: str, interpolation: str, downsample: int, with_colorbar: bool,
    tight: bool, origin: str,
    tile_width_bp: int, show_axes: bool, tick_mb: float,
    show_title: bool, scalebar_kb: Optional[float],
    start_str: Optional[str], end_str: Optional[str],
    write_sidecar: bool,
    # tracks
    #ctcf_bw_path: Optional[Path], atac_bw_path: Optional[Path],
    ctcf_bw_a: Optional[Path], atac_bw_a: Optional[Path],
    ctcf_bw_b: Optional[Path], atac_bw_b: Optional[Path],
    track_height: float, track_smooth_bins: int,
    #ctcf_label: str, atac_label: str,
    label_a: str, label_b: str,
    ctcf_color: str, atac_color: str,
) -> Tuple[Path, Optional[str]]:
    try:
        arrA = np.load(npy_a, mmap_mode='r')
        arrA = reshape_if_vector(arrA).astype(np.float64, copy=False)
        
        arrB = np.load(npy_b, mmap_mode='r')
        arrB = reshape_if_vector(arrB).astype(np.float64, copy=False)


        # Full tile geometry (assumes both matrices represent same tile width)
        N_full = arrA.shape[0]
        if arrB.shape[0] != N_full:
            raise ValueError(f"Raw shape mismatch: {arrA.shape} vs {arrB.shape} for {npy_a.name} / {npy_b.name}")
        
        bp_per_bin_full = tile_width_bp / float(N_full)
        
        # Sub-window crop before downsampling
        start_bp = parse_bp(start_str) if start_str is not None else 0
        end_bp   = parse_bp(end_str)   if end_str   is not None else tile_width_bp
        if start_bp < 0 or end_bp <= start_bp or end_bp > tile_width_bp:
            raise ValueError(f"Invalid --start/--end for tile of width {tile_width_bp} bp: start={start_bp}, end={end_bp}")
        
        i0 = int(math.floor(start_bp / bp_per_bin_full))
        i1 = int(math.ceil(end_bp  / bp_per_bin_full))
        i0 = max(0, min(i0, N_full))
        i1 = max(0, min(i1, N_full))
        if i1 - i0 < 2:
            raise ValueError("Requested window too small after rounding (fewer than 2 bins)")
        
        # crop BOTH
        if i0 != 0 or i1 != N_full:
            arrA = arrA[i0:i1, i0:i1]
            arrB = arrB[i0:i1, i0:i1]
        
        # Downsample BOTH (keep bp/bin consistent)
        bp_per_bin = bp_per_bin_full
        arrA = downsample_mean(arrA, downsample)
        arrB = downsample_mean(arrB, downsample)
        bp_per_bin *= downsample
        
        # Symmetrize BOTH (optional)
        if arrA.shape[0] == arrA.shape[1]:
            arrA = 0.5 * (arrA + arrA.T)
        if arrB.shape[0] == arrB.shape[1]:
            arrB = 0.5 * (arrB + arrB.T)
        
        # After crop/downsample they must match
        if arrA.shape != arrB.shape:
            raise ValueError(f"Processed shape mismatch: {arrA.shape} vs {arrB.shape} for {npy_a.name} / {npy_b.name}")

        # apply log if requested
        if log1p:
            arrA = np.log1p(arrA)
            arrB = np.log1p(arrB)
        
        # shared clipping (symmetric around 0 like you do)
        if clip_pct is not None:
            vals = np.concatenate([
                np.abs(arrA[np.isfinite(arrA)]).ravel(),
                np.abs(arrB[np.isfinite(arrB)]).ravel()
            ])
            hi = float(np.percentile(vals, clip_pct))
            arrA = np.clip(arrA, -hi, hi)
            arrB = np.clip(arrB, -hi, hi)
            vmin_final, vmax_final = -hi, hi
        else:
            # shared min/max if user didn’t force vmin/vmax
            vmin_final = vmin if vmin is not None else float(min(np.nanmin(arrA), np.nanmin(arrB)))
            vmax_final = vmax if vmax is not None else float(max(np.nanmax(arrA), np.nanmax(arrB)))


        N = arrA.shape[0]
        comp = np.empty((N, N), dtype=np.float64)
        
        iu = np.triu_indices(N, k=0)   # include diagonal
        il = np.tril_indices(N, k=-1)  # below diagonal only
        
        comp[iu] = arrA[iu]
        comp[il] = arrB[il]

        # Output path (annotate crop)
        out_dir.mkdir(parents=True, exist_ok=True)
        crop_tag = f"_{int(start_bp)}-{int(end_bp)}bp" if (start_bp != 0 or end_bp != tile_width_bp) else ""
        out_path = out_dir / f"{npy_a.stem}{crop_tag}.{fmt}"

        # Absolute coordinates
        chrom, tile_offset = parse_chr_offset(npy_a)
        abs_start = (tile_offset or 0) + start_bp
        abs_end   = (tile_offset or 0) + end_bp

        # ---------------- figure layout ----------------
        
        tracks = []
        if ctcf_bw_a: tracks.append((ctcf_bw_a, f"{label_a} CTCF", ctcf_color))
        if atac_bw_a: tracks.append((atac_bw_a, f"{label_a} ATAC", atac_color))
        if ctcf_bw_b: tracks.append((ctcf_bw_b, f"{label_b} CTCF", ctcf_color))
        if atac_bw_b: tracks.append((atac_bw_b, f"{label_b} ATAC", atac_color))
        n_tracks = len(tracks)

        fig_height = figsize + n_tracks * track_height + (0.25 if n_tracks else 0.0)

        # 2-column grid: left = heatmap + tracks, right = vertical colorbar
        fig = plt.figure(figsize=(figsize, fig_height))
        
        gs = fig.add_gridspec(
            nrows=1 + n_tracks,
            ncols=2,
            width_ratios=[1.0, 0.04],        # small column for colorbar
            height_ratios=[5] + [1] * n_tracks,
            wspace=0.05, hspace=0.12
        )

        # Heatmap (absolute Mb on both axes via extent)
        ax = fig.add_subplot(gs[0, 0])
        x0_mb, x1_mb = abs_start / 1e6, abs_end / 1e6
        extent = [x0_mb, x1_mb, x0_mb, x1_mb]

        # keep it square so the matrix doesn’t collapse
        #ax.set_aspect('equal', adjustable='box')
        im = ax.imshow(
            comp, origin=origin, interpolation=interpolation,
            vmin=vmin_final, vmax=vmax_final, cmap=cmap,
            extent=extent,
            aspect='auto'        # <-- key: no internal padding
        )
        #ax.set_box_aspect(1) 
        ax.set_aspect('auto')
        #ax.set_xlim(x0_mb, x1_mb)
        #ax.set_ylim(x0_mb, x1_mb)
        #ax.margins(x=0)
        ax.margins(0, 0)           # remove the 5% default padding on both axes
        ax.set_anchor('SW')
        ax.plot([x0_mb, x1_mb], [x0_mb, x1_mb], linewidth=0.8, alpha=0.6)
        ax.text(0.02, 0.98, f"{label_a} (upper)", transform=ax.transAxes, va='top', ha='left', fontsize=20)
        ax.text(0.8, 0.06, f"{label_b} (lower)", transform=ax.transAxes, va='top', ha='left', fontsize=20)

        # Axes in absolute Mb
#        if show_axes:
#            if tick_mb and tick_mb > 0:
#                start_tick = math.floor(x0_mb / tick_mb) * tick_mb
#                xticks = np.arange(start_tick, x1_mb + 1e-9, tick_mb)
#                ax.set_xticks(xticks)
#                ax.set_yticks(xticks)
#            label = f"{chrom} position (Mb)" if chrom else "Genomic position (Mb)"
#            #ax.set_xlabel("")
#            #ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
#            ax.set_xlabel(label)
#            ax.set_ylabel(label)
#        else:
#            ax.set_axis_off()
        if show_axes:
            if tick_mb and tick_mb > 0:
                start_tick = math.floor(x0_mb / tick_mb) * tick_mb
                xticks = np.arange(start_tick, x1_mb + 1e-9, tick_mb)
                ax.set_xticks(xticks)
                ax.set_yticks(xticks)
        
            x_label = f"{chrom} position (Mb)" if chrom else "Genomic position (Mb)"
        
            # Keep y label on heatmap
            ax.set_ylabel(x_label)
        
            # If we have tracks, move xlabel to bottom track later
            if n_tracks > 0:
                ax.set_xlabel("")
                ax.tick_params(axis='x', labelbottom=False)
            else:
                ax.set_xlabel(x_label)
        else:
            ax.set_axis_off()

        # Colorbar in its own column
        if with_colorbar:
            cax = fig.add_subplot(gs[0, 1])
            cb = fig.colorbar(im, cax=cax)
            cb.ax.tick_params(labelsize=8)

        # Title & scalebar
        if show_title:
            bpbin = int(round(bp_per_bin))
            ax.set_title(f"{chrom or npy_a.stem}:{abs_start:,}-{abs_end:,}  (≈{bpbin} bp/bin)", fontsize=11)
        if scalebar_kb:
            add_scalebar_mb(ax, scalebar_kb)

        # ---------------- tracks (left column, share x with heatmap) ----------------

        def _plot_track(ax_tr, bw_path: Optional[Path], label: str, color: str):
            if not bw_path:
                return
            bw = _open_bw(bw_path)
            try:
                n_bins = comp.shape[1]              # match heatmap width
                vals = bw_summary(bw, chrom or (parse_chr_offset(npy_a)[0] or ''), 
                                  abs_start, abs_end, n_bins)
            finally:
                bw.close()
        
            vals = smooth_1d(vals, track_smooth_bins)
        
            # --- plot by edges so the track ends exactly at x1 ---
            x0_mb, x1_mb = abs_start / 1e6, abs_end / 1e6
            edges = np.linspace(x0_mb, x1_mb, len(vals) + 1)     # N+1 edges for N bins
            ax_tr.fill_between(edges, np.r_[vals, vals[-1]], step='post',
                               color=color, linewidth=0)
        
            # align with Hi-C axes exactly
            ax_tr.set_xlim(x0_mb, x1_mb)
            #ax_tr.margins(x=0)                   # remove default 5% padding
            ax_tr.margins(0, 0)
            ax_tr.set_ylim(bottom=0)
            ax_tr.set_anchor('SW')
        
            ax_tr.set_ylabel(label, rotation=0, ha='right', va='center', labelpad=12)
            ax_tr.spines['top'].set_visible(False)
            ax_tr.spines['right'].set_visible(False)
            ax_tr.tick_params(axis='x', labelbottom=True)

#        row = 1
#        if ctcf_bw_path:
#            ax_ctcf = fig.add_subplot(gs[row, 0], sharex=ax)
#            _plot_track(ax_ctcf, ctcf_bw_path, ctcf_label, ctcf_color)
#            row += 1
#        if atac_bw_path:
#            ax_atac = fig.add_subplot(gs[row, 0], sharex=ax)
#            _plot_track(ax_atac, atac_bw_path, atac_label, atac_color)

#        row = 1
#        for j, (bw_path, lab, col) in enumerate(tracks):
#            ax_tr = fig.add_subplot(gs[row, 0], sharex=ax)
#            _plot_track(ax_tr, bw_path, lab, col)
#            if j == 1 and len(tracks) > 2:  # after A's second track
#                ax_tr.spines['bottom'].set_linewidth(2.0)
#                ax_tr.spines['bottom'].set_alpha(0.4)
#            row += 1

        # v3 tracks
        row = 1
        last_ax_tr = None
        for j, (bw_path, lab, col) in enumerate(tracks):
            ax_tr = fig.add_subplot(gs[row, 0], sharex=ax)
            _plot_track(ax_tr, bw_path, lab, col)
        
            # Hide x tick labels for all but the last track
            if j != len(tracks) - 1:
                ax_tr.tick_params(axis='x', labelbottom=False)
            else:
                if show_axes:
                    ax_tr.set_xlabel(x_label)
                    ax_tr.tick_params(axis='x', labelbottom=True)
        
            if j == 1 and len(tracks) > 2:
                ax_tr.spines['bottom'].set_linewidth(2.0)
                ax_tr.spines['bottom'].set_alpha(0.4)
        
            last_ax_tr = ax_tr
            row += 1


#        if len(tracks) >= 3:
#        # after plotting A tracks, before B tracks
#            ax_tr.axhline(ax_tr.get_ylim()[1], linewidth=1.2, alpha=0.3)
    
        if tight:
            plt.tight_layout()
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight', pad_inches=0.15)
        plt.close(fig)

        # --- sidecar stats ---
        if write_sidecar:
            stats = compute_tile_stats(comp, npy_a, bp_per_bin, vmin_final, vmax_final)
            
            stats["file_a"] = npy_a.name
            stats["file_b"] = npy_b.name
            
#            stats.update({
#                "source_tile_width_bp": int(tile_width_bp),
#                "window_start_bp": int(start_bp),
#                "window_end_bp": int(end_bp),
#                "abs_start_bp": int(abs_start),
#                "abs_end_bp": int(abs_end),
#                "downsample": int(downsample),
#                "ctcf_bw": str(ctcf_bw_path) if ctcf_bw_path else None,
#                "atac_bw": str(atac_bw_path) if atac_bw_path else None,
#            })
            stats.update({
                "source_tile_width_bp": int(tile_width_bp),
                "window_start_bp": int(start_bp),
                "window_end_bp": int(end_bp),
                "abs_start_bp": int(abs_start),
                "abs_end_bp": int(abs_end),
                "downsample": int(downsample),
                "ctcf_bw_a": str(ctcf_bw_a) if ctcf_bw_a else None,
                "atac_bw_a": str(atac_bw_a) if atac_bw_a else None,
                "ctcf_bw_b": str(ctcf_bw_b) if ctcf_bw_b else None,
                "atac_bw_b": str(atac_bw_b) if atac_bw_b else None,
                "label_a": label_a,
                "label_b": label_b,
            })

            write_sidecar_json(out_path, stats)

        return out_path, None
    except Exception as e:
        return npy_a, f"ERROR: {type(e).__name__}: {e}"

def index_tiles(files: List[Path]) -> dict:
    idx = {}
    for f in files:
        chrom, off = parse_chr_offset(f)
        idx[(chrom, off)] = f
    return idx

def pair_tiles(input_a: Path, pat_a: str, input_b: Path, pat_b: str) -> List[Tuple[Path, Path]]:
    fa = sort_tiles(discover_npy_files(input_a, pat_a))
    fb = sort_tiles(discover_npy_files(input_b, pat_b))
    ia, ib = index_tiles(fa), index_tiles(fb)
    keys = sorted(set(ia.keys()) & set(ib.keys()), key=lambda k: ((k[0] or ""), k[1]))
    pairs = [(ia[k], ib[k]) for k in keys]
    missing_a = set(ib.keys()) - set(ia.keys())
    missing_b = set(ia.keys()) - set(ib.keys())
    if missing_a:
        logging.warning("Missing %d tiles in A", len(missing_a))
    if missing_b:
        logging.warning("Missing %d tiles in B", len(missing_b))
    return pairs

# ----------------------------- CLI ------------------------------------- -------------------------------------

def main():
    p = argparse.ArgumentParser(description="Render .npy Hi‑C tiles to images with optional stacked bigWig tracks")
    p.add_argument('--input', type=Path, required=True, help='Directory containing .npy files')
    p.add_argument('--output', type=Path, required=True, help='Directory to write images')
    p.add_argument('--pattern', type=str, default='*.npy', help='Glob pattern for inputs (default: *.npy)')
    p.add_argument('--fmt', type=str, default='png', choices=['png', 'pdf', 'svg'], help='Image format (png|pdf|svg)')
    p.add_argument('--dpi', type=int, default=300, help='DPI for raster outputs (png)')
    p.add_argument('--figsize', type=float, default=16.0, help='Heatmap side length (inches)')
    p.add_argument('--interpolation', type=str, default='nearest', help='imshow interpolation (nearest, bilinear, etc.)')
    p.add_argument('--origin', type=str, default='lower', choices=['lower','upper'], help='Matrix origin for imshow')

    p.add_argument('--log1p', action='store_true', help='Apply log1p to intensities')
    p.add_argument('--clip-percentile', type=float, default=None, help='Clip intensities at this upper percentile (e.g., 99.5). Overrides vmin/vmax.')
    p.add_argument('--vmin', type=float, default=None, help='Explicit vmin (ignored if clip-percentile set)')
    p.add_argument('--vmax', type=float, default=None, help='Explicit vmax (ignored if clip-percentile set)')
    p.add_argument('--cmap', type=str, default='Reds', help='Matplotlib colormap (e.g., Reds, viridis, magma)')

    p.add_argument('--downsample', type=int, default=1, help='Integer mean‑pooling factor for downsampling (>=1)')
    p.add_argument('--workers', type=int, default=1, help='Parallel workers (use 1 if you hit mpl fork issues)')
    p.add_argument('--limit', type=int, default=None, help='Process only the first N files (for testing)')
    p.add_argument('--colorbar', dest='with_colorbar', action='store_true', help='Add a colorbar to each figure')
    p.add_argument('--no-colorbar', dest='with_colorbar', action='store_false')
    p.set_defaults(with_colorbar=False)
    p.add_argument('--tight', action='store_true', help='Use tight_layout and small padding')
    p.add_argument('--log-level', type=str, default='INFO', help='Logging level (DEBUG, INFO, WARNING, ERROR)')

    # Coordinate‑aware options
    p.add_argument('--tile-width-bp', type=int, default=2_097_152, help='Biological width covered by each tile (bp)')
    p.add_argument('--start', type=str, default=None, help='Start of subwindow within the tile (e.g., 250kb, 1.2mb). Default: 0')
    p.add_argument('--end', type=str, default=None, help='End of subwindow within the tile. Default: full tile width')
    p.add_argument('--axes', dest='show_axes', action='store_true', help='Show Mb axes and labels (absolute coordinates)')
    p.add_argument('--no-axes', dest='show_axes', action='store_false')
    p.set_defaults(show_axes=True)
    p.add_argument('--tick-mb', type=float, default=0.5, help='Tick interval in Mb when axes are shown')
    p.add_argument('--title', dest='show_title', action='store_true', help='Show a title with coordinates and bp/bin')
    p.add_argument('--no-title', dest='show_title', action='store_false')
    p.set_defaults(show_title=True)
    p.add_argument('--scalebar-kb', type=float, default=100.0, help='Draw a scalebar of this length (kb)')

    # Stats options
    p.add_argument('--no-sidecar', dest='write_sidecar', action='store_false', help='Disable writing <image>.json sidecars')
    p.add_argument('--sidecar', dest='write_sidecar', action='store_true')
    p.set_defaults(write_sidecar=True)
    p.add_argument('--stats-csv', type=Path, default=None, help='If set, aggregate all sidecar JSONs into this CSV after rendering')

    # Track options (two common tracks for convenience)
    p.add_argument('--ctcf-bw', type=Path, default=None, help='bigWig path for predicted/observed CTCF track')
    p.add_argument('--atac-bw', type=Path, default=None, help='bigWig path for ATAC‑seq (normalized)')
    p.add_argument('--ctcf-color', type=str, default='#111111')
    p.add_argument('--atac-color', type=str, default='#222222')
    p.add_argument('--track-height', type=float, default=1.2, help='Height (inches) per track panel')
    p.add_argument('--track-smooth-bins', type=int, default=3, help='Moving‑average smoothing window (bins) for tracks')

    # second condition (B)
    p.add_argument('--input-b', type=Path, required=True, help='Directory containing condition B .npy files')
    p.add_argument('--pattern-b', type=str, default=None, help='Glob for condition B (default: same as --pattern)')
    
    p.add_argument('--ctcf-bw-b', type=Path, default=None)
    p.add_argument('--atac-bw-b', type=Path, default=None)
    p.add_argument('--label-a', type=str, default='A')
    p.add_argument('--label-b', type=str, default='B')
    args = p.parse_args()
    
    if args.pattern_b is None:
        args.pattern_b = args.pattern

    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO),
                        format='[%(levelname)s] %(message)s')

    # 2 numpy matrices are iterated through to compare
    pairs = pair_tiles(args.input, args.pattern, args.input_b, args.pattern_b)

    if args.limit is not None:
        pairs = pairs[:args.limit]
    if not pairs:
        return
    
    logging.info("Found %d matched tile pairs", len(pairs))
    out_dir: Path = args.output
    if args.workers <= 1:
        errors = []
        for i, (fa, fb) in enumerate(pairs, 1):
            logging.info("[%d/%d] %s vs %s", i, len(pairs), fa.name, fb.name)
            out, err = render_tile(
                fa, fb, out_dir, args.fmt, args.dpi, args.figsize,
                args.log1p, args.clip_percentile, args.vmin, args.vmax,
                args.cmap, args.interpolation, args.downsample,
                args.with_colorbar, args.tight, args.origin,
                args.tile_width_bp, args.show_axes, args.tick_mb,
                args.show_title, args.scalebar_kb,
                args.start, args.end,
                args.write_sidecar,
                args.ctcf_bw, args.atac_bw,
                args.ctcf_bw_b, args.atac_bw_b,
                args.track_height, args.track_smooth_bins,
                args.label_a, args.label_b,
                args.ctcf_color, args.atac_color,
            )
            if err:
                logging.error("%s", err)
                errors.append(((fa,fb), err))
        if errors:
            logging.warning("Completed with %d errors", len(errors))
#    else:
#        errors = []
#        with ProcessPoolExecutor(max_workers=args.workers) as ex:
#            futs = {
#                ex.submit(
#                    render_tile,
#                    fa, fb, out_dir, args.fmt, args.dpi, args.figsize,
#                    args.log1p, args.clip_percentile, args.vmin, args.vmax,
#                    args.cmap, args.interpolation, args.downsample,
#                    args.with_colorbar, args.tight, args.origin,
#                    args.tile_width_bp, args.show_axes, args.tick_mb,
#                    args.show_title, args.scalebar_kb,
#                    args.start, args.end,
#                    args.write_sidecar,
#                    # tracks A
#                    args.ctcf_bw, args.atac_bw,
#                    # tracks B
#                    args.ctcf_bw_b, args.atac_bw_b,
#                    args.track_height, args.track_smooth_bins,
#                    args.label_a, args.label_b,
#                    args.ctcf_color, args.atac_color,
#                ): (fa, fb)
#                for (fa, fb) in pairs
#            }
#
#        for i, fut in enumerate(as_completed(futs), 1):
#            fa, fb = futs[fut]
#            try:
#                out, err = fut.result()
#                logging.info(
#                    "[%d/%d] %s vs %s -> %s",
#                    i, len(pairs),
#                    fa.name, fb.name,
#                    out.name if isinstance(out, Path) else out
#                )
#                if err:
#                    logging.error("%s", err)
#                    errors.append(((fa, fb), err))
#            except Exception as e:
#                logging.error("Worker error on %s vs %s: %s", fa.name, fb.name, e)
#                errors.append(((fa, fb), str(e)))
#
#    if errors:
#        logging.warning("Completed with %d errors", len(errors))

    else:
        errors = []
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = {
                ex.submit(
                    render_tile,
                    fa, fb, out_dir, args.fmt, args.dpi, args.figsize,
                    args.log1p, args.clip_percentile, args.vmin, args.vmax,
                    args.cmap, args.interpolation, args.downsample,
                    args.with_colorbar, args.tight, args.origin,
                    args.tile_width_bp, args.show_axes, args.tick_mb,
                    args.show_title, args.scalebar_kb,
                    args.start, args.end,
                    args.write_sidecar,
                    # tracks A
                    args.ctcf_bw, args.atac_bw,
                    # tracks B
                    args.ctcf_bw_b, args.atac_bw_b,
                    args.track_height, args.track_smooth_bins,
                    args.label_a, args.label_b,
                    args.ctcf_color, args.atac_color,
                ): (fa, fb)
                for (fa, fb) in pairs
            }
    
            for i, fut in enumerate(as_completed(futs), 1):
                fa, fb = futs[fut]
                try:
                    out, err = fut.result()
                    logging.info(
                        "[%d/%d] %s vs %s -> %s",
                        i, len(pairs),
                        fa.name, fb.name,
                        out.name if isinstance(out, Path) else out
                    )
                    if err:
                        logging.error("%s", err)
                        errors.append(((fa, fb), err))
                except Exception as e:
                    logging.error("Worker error on %s vs %s: %s", fa.name, fb.name, e)
                    errors.append(((fa, fb), str(e)))
    
        if errors:
            logging.warning("Completed with %d errors", len(errors))


    # Post‑pass: aggregate CSV if requested
    if args.stats_csv:
        try:
            collect_stats_csv_from_sidecars(out_dir, args.stats_csv)
        except Exception as e:
            logging.error("Failed to write stats CSV: %s", e)

if __name__ == "__main__":
    main()
