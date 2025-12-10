import argparse
import csv
import os
import re

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

# ---------- nominal targets ----------
TARGET_A = 3e-3   # 0.003 m
TARGET_BPM = 40.0 # 40 BPM
F_TARGET_HZ = TARGET_BPM / 60.0  # 40 BPM -> 0.666... Hz


# ---------- basic IO ----------

def load_distance_csv(path):
    """Return elapsed time (s) and distance (m) from one CSV."""
    t_ms = []
    d_m = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        if "elapsed_ms" not in reader.fieldnames or "m_dist" not in reader.fieldnames:
            raise ValueError(
                f"{path}: header must contain 'elapsed_ms' and 'm_dist', "
                f"got {reader.fieldnames}"
            )
        for row in reader:
            try:
                t_ms.append(float(row["elapsed_ms"]))
                d_m.append(float(row["m_dist"]))
            except (ValueError, TypeError):
                continue
    t_ms = np.asarray(t_ms)
    d_m = np.asarray(d_m)
    idx = np.argsort(t_ms)
    return t_ms[idx] / 1000.0, d_m[idx]


def crop_first_5min(t, d, duration_limit=300.0):
    """Keep only the first 5 minutes (or less if the recording is shorter)."""
    t0 = t[0]
    t_end = t0 + duration_limit
    mask = t <= t_end
    return t[mask], d[mask]


def resample_uniform(t, x):
    """Resample x(t) to a uniform grid, returning t_u, x_u, fs."""
    dt_med = np.median(np.diff(t))
    fs = 1.0 / dt_med
    t_u = np.arange(t[0], t[-1], dt_med)
    x_u = np.interp(t_u, t, x)
    return t_u, x_u, fs


# ---------- band-pass + rescale ----------

def bandpass_filter(x, fs, fmin, fmax, order=4):
    """Zero-phase Butterworth band-pass."""
    nyq = 0.5 * fs
    low = fmin / nyq
    high = fmax / nyq
    b, a = signal.butter(order, [low, high], btype="bandpass")
    return signal.filtfilt(b, a, x)


def filter_and_match_peak_to_peak(t, d, fmin, fmax):
    """
    1) Crop to first 5 minutes.
    2) Resample to uniform grid.
    3) Band-pass filter.
    4) Rescale filtered signal so its peak-to-peak equals the original
       (over the same 5-min window).
    Returns t_u, d_proc, fs, ptp_orig.
    """
    # crop to 5 minutes in the original irregularly spaced data
    t5, d5 = crop_first_5min(t, d, duration_limit=300.0)

    # peak-to-peak of original CSV data (5-min window)
    ptp_orig = float(d5.max() - d5.min())

    # resample and band-pass
    t_u, d_u, fs = resample_uniform(t5, d5)
    d_bp = bandpass_filter(d_u - np.mean(d_u), fs, fmin, fmax)

    # peak-to-peak of filtered signal
    ptp_filt = float(d_bp.max() - d_bp.min())
    if ptp_filt > 0:
        scale = ptp_orig / ptp_filt
    else:
        scale = 1.0

    d_proc = d_bp * scale
    return t_u, d_proc, fs, ptp_orig


# ---------- BPM estimate from FFT ----------

def estimate_bpm_fft(x, fs, fmin, fmax):
    """Estimate BPM by FFT peak within [fmin, fmax]."""
    N = len(x)
    freqs = np.fft.rfftfreq(N, d=1.0 / fs)
    X = np.fft.rfft(x)
    power = np.abs(X) ** 2

    band = (freqs >= fmin) & (freqs <= fmax)
    if not np.any(band):
        return None, None

    idx_peak = np.argmax(power[band])
    f_peak = freqs[band][idx_peak]
    bpm = 60.0 * f_peak
    return f_peak, bpm


# ---------- sine fitting ----------

def fit_sine_same_freq(t_u, x, f0):
    """
    Fit x(t) ≈ A * sin(2π f0 t + phi) via least-squares
    on the FULL 5-min processed signal.
    Returns A, phi.
    """
    t0 = t_u - t_u[0]
    s_sin = np.sin(2 * np.pi * f0 * t0)
    s_cos = np.cos(2 * np.pi * f0 * t0)

    M = np.vstack([s_sin, s_cos]).T
    params, *_ = np.linalg.lstsq(M, x, rcond=None)
    a, b = params

    A = np.sqrt(a ** 2 + b ** 2)
    phi = np.arctan2(b, a)
    return A, phi


# ---------- metadata parsing from filename (video quality) ----------

def parse_video_metadata_from_path(path):
    """
    Parse resolution (p) and nominal fps from filename.

    Supports patterns like:
      distance_log_720p_5fps.csv
      distance_log_720p_60fps.csv
      distance_log_1080p_30fps.csv
    """
    base = os.path.basename(path).lower()

    resolution_p = None
    fps = None

    # resolution: _720p_, _1080p_, etc.
    m = re.search(r'_(\d{3,4})p', base)
    if m:
        try:
            resolution_p = int(m.group(1))
        except ValueError:
            resolution_p = None

    # fps: _5fps, _30fps, etc.
    m = re.search(r'_(\d+)fps', base)
    if m:
        try:
            fps = int(m.group(1))
        except ValueError:
            fps = None

    return resolution_p, fps


# ---------- main per-file pipeline ----------

def analyze_file(path, fmin, fmax, t_plot):
    # 1) load and filter+rescale (5-min window)
    t, d = load_distance_csv(path)
    t_u, d_proc, fs, ptp_orig = filter_and_match_peak_to_peak(t, d, fmin, fmax)
    duration = t_u[-1] - t_u[0]
    print(f"\n=== {path} ===")
    print(f"  5-min window duration actually used: {duration:.1f} s")
    print(f"  Peak-to-peak original (5 min): {ptp_orig:.4e} m")

    # parse resolution and fps from filename
    resolution_p, fps_nominal = parse_video_metadata_from_path(path)

    # 2) BPM from FFT on full processed signal
    f0, bpm_est = estimate_bpm_fft(d_proc, fs, fmin, fmax)
    if f0 is None:
        print("  [WARN] No frequency in the given band.")
        return None

    bpm_err_pct = (bpm_est - TARGET_BPM) / TARGET_BPM * 100.0
    print(f"  BPM_est ≈ {bpm_est:.2f} (f_peak={f0:.4f} Hz)")
    print(f"  BPM error ≈ {bpm_err_pct:.2f}% vs {TARGET_BPM:.1f} BPM")

    # time offset starting at 0
    t0 = t_u - t_u[0]

    # -------- Figure 1: processed data vs fixed 40 BPM, 0.003 m template --------
    mask_plot = t0 <= t_plot
    t_plot_arr = t0[mask_plot]
    x_plot = d_proc[mask_plot]

    # fixed template: 40 BPM, 0.003 m, zero phase
    template_fixed = TARGET_A * np.sin(2 * np.pi * F_TARGET_HZ * t_plot_arr + 0.0)

    plt.figure(figsize=(10, 4))
    title1 = (f"{os.path.basename(path)} | BPM_est≈{bpm_est:.2f} "
              f"(err={bpm_err_pct:.2f}% vs {TARGET_BPM:.1f})")
    plt.title(title1)
    plt.plot(t_plot_arr, x_plot, label="Processed distance (band-pass + rescaled)")
    plt.plot(t_plot_arr, template_fixed, "--",
             label="Template: 40 BPM, 0.003 m, phase=0")
    plt.xlabel("Time [s]")
    plt.ylabel("Distance [m]")
    plt.legend(loc="upper right")
    plt.grid(True)

    # -------- Fit sine with same freq as data, get amplitude & phase --------
    A_data, phi_data = fit_sine_same_freq(t_u, d_proc, f0)
    amp_err_pct = (A_data - TARGET_A) / TARGET_A * 100.0
    print(f"  Data amplitude A_est ≈ {A_data:.4e} m")
    print(f"  Amplitude error ≈ {amp_err_pct:.2f}% vs {TARGET_A:.4f} m")
    print(f"  Phase phi_est ≈ {phi_data:.3f} rad")

    # -------- Figure 2: same data vs aligned template (same f, same phase) --------
    template_aligned = TARGET_A * np.sin(2 * np.pi * f0 * t_plot_arr + phi_data)

    plt.figure(figsize=(10, 4))
    title2 = (f"{os.path.basename(path)} | Aligned template (f=f_peak, phase=phi_est)\n"
              f"A_est={A_data:.4e} m (err={amp_err_pct:.2f}% vs 0.003 m)")
    plt.title(title2)
    plt.plot(t_plot_arr, x_plot, label="Processed distance (same as Fig.1)")
    plt.plot(t_plot_arr, template_aligned, "--",
             label="Template: 0.003 m, f=f_peak, phase aligned")
    plt.xlabel("Time [s]")
    plt.ylabel("Distance [m]")
    plt.legend(loc="upper right")
    plt.grid(True)

    # summary row for CSV
    row = {
        "file": os.path.basename(path),
        "resolution_p": resolution_p,
        "fps_nominal": fps_nominal,
        "duration_used_s": duration,
        "fs_Hz": fs,

        # BPM values
        "BPM_est": bpm_est,
        "BPM_err_pct": f"{bpm_err_pct}%",

        # Amplitude values
        "A_est_m": A_data,
        "A_err_pct": f"{amp_err_pct}%"
    }

    return row


# ---------- CLI ----------

def main():
    parser = argparse.ArgumentParser(
        description="Breathing analysis for video-quality test: BPM error + amplitude error with aligned template"
    )
    parser.add_argument(
        "--csv", nargs="+", required=True,
        help="One or more CSV files (video-quality test: distance_log_720p_5fps.csv, etc.)"
    )
    parser.add_argument(
        "--fmin", type=float, default=0.4,
        help="Lower cutoff / search frequency in Hz (default 0.4)"
    )
    parser.add_argument(
        "--fmax", type=float, default=1.0,
        help="Upper cutoff / search frequency in Hz (default 1.0)"
    )
    parser.add_argument(
        "--t-plot", type=float, default=30.0,
        help="Seconds to display from start (default 30); "
             "analysis always uses first 5 minutes."
    )
    parser.add_argument(
        "--out-csv", type=str, default=None,
        help="Optional path to write summary CSV (for Excel)"
    )
    args = parser.parse_args()

    rows = []
    for path in args.csv:
        row = analyze_file(path, args.fmin, args.fmax, args.t_plot)
        if row is not None:
            rows.append(row)

    if args.out_csv is not None and rows:
        fieldnames = [
            "file",
            "resolution_p",
            "fps_nominal",
            "duration_used_s",
            "fs_Hz",
            "BPM_est",
            "BPM_err_pct",
            "A_est_m",
            "A_err_pct",
        ]
        with open(args.out_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in rows:
                writer.writerow(r)
        print(f"\n[INFO] Wrote summary CSV to: {args.out_csv}")

    plt.show()


if __name__ == "__main__":
    main()
