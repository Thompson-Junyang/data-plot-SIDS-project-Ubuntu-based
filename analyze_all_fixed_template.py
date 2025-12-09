#!/usr/bin/env python3
import argparse
import csv
import glob
import math
import os
import re

import numpy as np


def load_csv_elapsed_and_dist(path):
    """
    Load one CSV of the form:
        timestamp_ms, elapsed_ms, pix_dist, m_dist
    Returns:
        t_sec  : np.array of elapsed time in seconds (starting at 0)
        d_m    : np.array of distances in meters (m_dist)
    Skips rows with empty or non-numeric m_dist.
    """
    t_ms = []
    d_m = []

    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        if header is None:
            raise ValueError(f"No header found in {path}")
        required = {"elapsed_ms", "m_dist"}
        if not required.issubset(set(header)):
            raise ValueError(
                f"{path} missing required columns {required}, header={header}"
            )

        for row in reader:
            e = row.get("elapsed_ms", "").strip()
            m = row.get("m_dist", "").strip()
            if e == "" or m == "":
                continue
            try:
                e_val = float(e)
                m_val = float(m)
            except ValueError:
                continue
            t_ms.append(e_val)
            d_m.append(m_val)

    if len(t_ms) == 0:
        raise ValueError(f"No usable rows in {path}")

    t_ms = np.array(t_ms, dtype=np.float64)
    d_m = np.array(d_m, dtype=np.float64)

    # Convert to seconds and set first time to 0
    t_sec = (t_ms - t_ms[0]) / 1000.0
    return t_sec, d_m


def resample_uniform(t, x):
    """
    Interpolate (t, x) to a uniform grid for FFT.

    Returns:
        t_u : uniform time grid (seconds)
        x_u : resampled signal
        fs  : sampling rate (Hz)
    """
    t = np.asarray(t, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)

    # Use median dt as robust estimate
    dt_est = np.median(np.diff(t))
    if dt_est <= 0:
        raise ValueError("Non-increasing time stamps detected.")
    fs = 1.0 / dt_est

    t0 = t[0]
    t1 = t[-1]
    n_samples = int(np.floor((t1 - t0) * fs)) + 1
    t_u = t0 + np.arange(n_samples) / fs

    x_u = np.interp(t_u, t, x)
    return t_u, x_u, fs


def bandpass_fft(x, fs, fmin, fmax):
    """
    Simple band-pass filter in the frequency domain.
    Zeroes out frequencies outside [fmin, fmax].
    """
    x = np.asarray(x, dtype=np.float64)
    n = len(x)
    X = np.fft.rfft(x)
    freqs = np.fft.rfftfreq(n, d=1.0 / fs)

    mask = (freqs >= fmin) & (freqs <= fmax)
    X_bp = X * mask
    x_bp = np.fft.irfft(X_bp, n=n)
    return x_bp, freqs, X_bp


def analyze_one_file(path, fmin, fmax, bpm_ref):
    """
    Analyze a single CSV file:
      - load data
      - resample
      - band-pass filter
      - compute FFT-based BPM
      - compute correlation with a fixed-amplitude sine template
      - compute SNR around f_ref
    Returns a dict of metrics.
    """
    t_sec, d_m = load_csv_elapsed_and_dist(path)
    duration = t_sec[-1] - t_sec[0]

    # Resample to uniform grid
    t_u, d_u, fs = resample_uniform(t_sec, d_m)

    # Remove overall mean before filtering (optional)
    d_u = d_u - np.mean(d_u)

    # Band-pass in [fmin, fmax]
    d_bp, freqs, X_bp = bandpass_fft(d_u, fs, fmin, fmax)

    # Power spectrum (only band-limited is left in X_bp)
    power = np.abs(X_bp) ** 2

    # Dominant frequency from the band-limited spectrum
    # (ignore 0 Hz if present)
    band_mask = (freqs >= fmin) & (freqs <= fmax)
    band_mask &= (freqs > 0.0)
    if not np.any(band_mask):
        f_peak = np.nan
        bpm_fft = np.nan
    else:
        idx_peak = np.argmax(power[band_mask])
        peak_indices = np.where(band_mask)[0]
        idx_global = peak_indices[idx_peak]
        f_peak = freqs[idx_global]
        bpm_fft = 60.0 * f_peak

    # --- Fixed template sine: same frequency & amplitude for all files ---
    f_ref = bpm_ref / 60.0
    # Use same time grid t_u and a unit-amplitude sine
    s = np.sin(2.0 * np.pi * f_ref * t_u)

    # Normalize both signals to zero mean + unit std for correlation
    def zscore(z):
        z = np.asarray(z, dtype=np.float64)
        z = z - np.mean(z)
        std = np.std(z)
        if std <= 0:
            return np.zeros_like(z)
        return z / std

    s_z = zscore(s)
    d_z = zscore(d_bp)

    corr = float(np.mean(s_z * d_z))

    # --- SNR around f_ref ---
    # Use power in a narrow band around f_ref as "signal",
    # remaining band power as "noise".
    df = freqs[1] - freqs[0] if len(freqs) > 1 else 0.0
    # ~ ±0.05 Hz window around f_ref
    width_hz = 0.05
    signal_mask = (freqs >= (f_ref - width_hz)) & (freqs <= (f_ref + width_hz))
    signal_mask &= band_mask

    noise_mask = band_mask & (~signal_mask)

    P_signal = float(np.sum(power[signal_mask])) if np.any(signal_mask) else 0.0
    P_noise = float(np.sum(power[noise_mask])) if np.any(noise_mask) else 0.0

    if P_signal > 0.0 and P_noise > 0.0:
        snr_db = 10.0 * math.log10(P_signal / P_noise)
    elif P_signal > 0.0 and P_noise == 0.0:
        snr_db = float("inf")
    else:
        snr_db = float("nan")

    return {
        "file": os.path.basename(path),
        "duration_s": duration,
        "fs": fs,
        "bpm_fft": bpm_fft,
        "corr_sine": corr,
        "snr_db": snr_db,
        "f_peak": f_peak,
    }


def extract_distance_cm(filename):
    """
    Try to extract something like '50' from 'distance_log_50cm.csv'.
    Returns float or None.
    """
    m = re.search(r"(\d+(?:\.\d+)?)\s*cm", filename)
    if not m:
        # fallback: any number in the name
        m = re.search(r"(\d+(?:\.\d+)?)", filename)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            return None
    return None


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Analyze multiple distance_log_*.csv files with a fixed 40 BPM "
            "unit-amplitude sine template and band-pass filtering."
        )
    )
    parser.add_argument(
        "--dir",
        type=str,
        default=".",
        help="Directory containing distance_log_*.csv (default: current dir).",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="distance_log_*cm.csv",
        help="Glob pattern for CSVs (default: distance_log_*cm.csv).",
    )
    parser.add_argument(
        "--fmin",
        type=float,
        default=0.4,
        help="Band-pass lower cutoff in Hz (default: 0.4).",
    )
    parser.add_argument(
        "--fmax",
        type=float,
        default=1.0,
        help="Band-pass upper cutoff in Hz (default: 1.0).",
    )
    parser.add_argument(
        "--bpm_ref",
        type=float,
        default=40.0,
        help="Reference BPM for fixed sine template (default: 40).",
    )

    args = parser.parse_args()

    pattern = os.path.join(args.dir, args.pattern)
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"No files matched pattern: {pattern}")
        return

    print(f"Found {len(files)} file(s) matching {pattern}")
    print(
        f"Using band-pass [{args.fmin:.2f}, {args.fmax:.2f}] Hz "
        f"and fixed template BPM = {args.bpm_ref:.2f}"
    )
    print()

    results = []
    for path in files:
        try:
            res = analyze_one_file(path, args.fmin, args.fmax, args.bpm_ref)
            res["distance_cm"] = extract_distance_cm(res["file"])
            results.append(res)
        except Exception as e:
            print(f"[WARN] Failed to analyze {path}: {e}")

    # Sort by distance if available
    results.sort(
        key=lambda r: float("inf")
        if r["distance_cm"] is None
        else r["distance_cm"]
    )

    # Print table header
    header = (
        "File".ljust(20),
        "Dist(cm)".rjust(8),
        "Dur(s)".rjust(8),
        "BPM_FFT".rjust(10),
        "Corr(sine)".rjust(12),
        "SNR(dB)".rjust(10),
    )
    print(" ".join(header))
    print("-" * 72)

    for r in results:
        fname = r["file"]
        dcm = r["distance_cm"]
        dur = r["duration_s"]
        bpm = r["bpm_fft"]
        corr = r["corr_sine"]
        snr = r["snr_db"]

        dist_str = f"{dcm:.0f}" if dcm is not None else "NA"
        dur_str = f"{dur:7.1f}"
        bpm_str = "   NaN" if math.isnan(bpm) else f"{bpm:9.2f}"
        corr_str = "     NaN" if math.isnan(corr) else f"{corr:11.3f}"
        if math.isinf(snr):
            snr_str = "     inf"
        elif math.isnan(snr):
            snr_str = "     NaN"
        else:
            snr_str = f"{snr:9.2f}"

        print(
            f"{fname.ljust(20)} {dist_str.rjust(8)} "
            f"{dur_str} {bpm_str} {corr_str} {snr_str}"
        )


if __name__ == "__main__":
    main()
