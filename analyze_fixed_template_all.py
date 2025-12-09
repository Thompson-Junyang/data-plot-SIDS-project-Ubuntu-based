# python3 analyze_fixed_template_all.py     --csv distance_log_50cm.csv distance_log_70cm.csv distance_log_90cm.csv distance_log_110cm.csv distance_log_130cm.csv distance_log_150cm.csv distance_log_170cm.csv distance_log_190cm.csv distance_log_210cm.csv distance_log_230cm.csv distance_log_250cm.csv distance_log_270cm.csv distance_log_290cm.csv    --fmin 0.4 --fmax 1.0     --t-plot 999999





# #!/usr/bin/env python3
# import argparse
# import csv
# import os

# import numpy as np
# import matplotlib.pyplot as plt


# def load_csv_elapsed_mdist(path):
#     """
#     Load elapsed_ms (ms) and m_dist (meters) from your CSV.

#     Expects header with at least:
#         elapsed_ms, m_dist
#     """
#     t_ms = []
#     d = []

#     with open(path, "r", newline="", encoding="utf-8") as f:
#         reader = csv.DictReader(f)
#         header = reader.fieldnames
#         print(f"[INFO] {os.path.basename(path)} header: {header}")

#         if "elapsed_ms" not in header or "m_dist" not in header:
#             raise ValueError(f"{path}: CSV does not have required columns "
#                              f"'elapsed_ms' and 'm_dist'.")

#         for row in reader:
#             if row["m_dist"] is None or row["m_dist"] == "":
#                 continue
#             try:
#                 t_ms.append(float(row["elapsed_ms"]))
#                 d.append(float(row["m_dist"]))
#             except ValueError:
#                 # skip malformed rows
#                 continue

#     t_ms = np.array(t_ms, dtype=np.float64)
#     d = np.array(d, dtype=np.float64)

#     # Convert ms to seconds; make time start at 0
#     t_s = (t_ms - t_ms[0]) / 1000.0
#     return t_s, d


# def resample_uniform(t, x):
#     """
#     Interpolate onto a uniform time grid using the median dt
#     (good enough for your ~31 Hz data).
#     """
#     if len(t) < 2:
#         raise ValueError("Not enough samples to resample.")

#     dt_est = np.median(np.diff(t))
#     t_u = np.arange(t[0], t[-1], dt_est)
#     x_u = np.interp(t_u, t, x)
#     fs = 1.0 / dt_est
#     return t_u, x_u, fs


# def bandpass_fft(x, fs, fmin, fmax):
#     """
#     Simple FFT band-pass filter:
#         - x: 1D signal
#         - fs: sampling rate (Hz)
#         - [fmin, fmax]: pass band (Hz)
#     """
#     N = len(x)
#     X = np.fft.rfft(x)
#     freqs = np.fft.rfftfreq(N, d=1.0/fs)

#     mask = (freqs >= fmin) & (freqs <= fmax)
#     X_bp = X * mask

#     x_bp = np.fft.irfft(X_bp, n=N)
#     power = (np.abs(X_bp) ** 2) / N
#     return x_bp, freqs, power


# def build_fixed_template(t, bpm, A_ref=1.0):
#     """
#     Fixed sine template: amplitude A_ref, frequency from BPM.
#     """
#     f0 = bpm / 60.0  # Hz
#     phase0 = 0.0
#     s = A_ref * np.sin(2.0 * np.pi * f0 * t + phase0)
#     return s, f0


# def compute_corr_and_snr(x, s):
#     """
#     x: filtered data (1D)
#     s: template (1D, same length)

#     Returns:
#         corr: Pearson correlation between x and s
#         snr_db: SNR in dB using projection of x onto s as "signal"
#     """
#     # --- Correlation (on normalized signals) ---
#     x_norm = x - np.mean(x)
#     s_norm = s - np.mean(s)

#     sx = np.std(x_norm)
#     ss = np.std(s_norm)
#     if sx == 0 or ss == 0:
#         corr = np.nan
#     else:
#         corr = np.mean((x_norm / sx) * (s_norm / ss))

#     # --- SNR via projection ---
#     # Remove mean before projection
#     x0 = x_norm
#     s0 = s_norm

#     denom = np.dot(s0, s0)
#     if denom == 0:
#         snr_db = np.nan
#     else:
#         alpha = np.dot(x0, s0) / denom     # best scaling for template
#         x_hat = alpha * s0                 # signal estimate
#         noise = x0 - x_hat

#         Ps = np.mean(x_hat ** 2)
#         Pn = np.mean(noise ** 2)
#         if Pn <= 0:
#             snr_db = np.inf
#         else:
#             snr_db = 10.0 * np.log10(Ps / Pn)

#     return corr, snr_db


# def plot_time_with_template(t, x_bp, s, title, t_plot=None):
#     """
#     Plot band-passed signal x_bp and template s in time domain.
#     If t_plot is given (seconds), only show [0, t_plot].
#     """
#     if t_plot is not None:
#         mask = t <= t_plot
#         t = t[mask]
#         x_bp = x_bp[mask]
#         s = s[mask]

#     plt.figure(figsize=(10, 4))
#     plt.plot(t, x_bp, label="Band-passed distance", linewidth=1.0)
#     plt.plot(t, s, label="Fixed sine template", linestyle="--", linewidth=1.0)
#     plt.xlabel("Time [s]")
#     plt.ylabel("Distance / Template (a.u.)")
#     plt.title(title)
#     plt.legend()
#     plt.grid(True)


# def main():
#     parser = argparse.ArgumentParser(
#         description="Compare CSV distance signals to a fixed sine template "
#                     "with band-pass filtering.")
#     parser.add_argument(
#         "--csv", nargs="+", required=True,
#         help="List of CSV files (each with elapsed_ms, m_dist).")
#     parser.add_argument(
#         "--fmin", type=float, default=0.4,
#         help="Lower band-pass cutoff in Hz (default: 0.4).")
#     parser.add_argument(
#         "--fmax", type=float, default=1.0,
#         help="Upper band-pass cutoff in Hz (default: 1.0).")
#     parser.add_argument(
#         "--template-bpm", type=float, default=40.0,
#         help="BPM for fixed sine template (same for all CSVs).")
#     parser.add_argument(
#         "--template-A", type=float, default=1.0,
#         help="Amplitude for fixed sine template (same for all CSVs).")
#     parser.add_argument(
#         "--t-plot", type=float, default=30.0,
#         help="Seconds of data to show in plots (default: 30). "
#              "Use a large number (e.g., 9999) if you want full length.")
#     args = parser.parse_args()

#     print(f"[INFO] Using fixed template: BPM={args.template_bpm}, "
#           f"A={args.template_A}, band=[{args.fmin}, {args.fmax}] Hz\n")

#     for path in args.csv:
#         print(f"=== {path} ===")
#         # 1) Load
#         t, d = load_csv_elapsed_mdist(path)
#         print(f"  Loaded {len(d)} samples; duration ≈ {t[-1] - t[0]:.1f} s")

#         # Remove overall mean first
#         d_centered = d - np.mean(d)

#         # 2) Resample to uniform grid
#         t_u, d_u, fs = resample_uniform(t, d_centered)
#         print(f"  Resampled: fs ≈ {fs:.2f} Hz, N = {len(d_u)}")

#         # 3) Band-pass filter
#         d_bp, freqs, power = bandpass_fft(d_u, fs, args.fmin, args.fmax)

#         # 4) Build fixed template (same amplitude & BPM for ALL files)
#         s, f0 = build_fixed_template(t_u, args.template_bpm, A_ref=args.template_A)

#         # 5) Compute correlation & SNR
#         corr, snr_db = compute_corr_and_snr(d_bp, s)
#         print(f"  Corr with template: {corr:.3f}")
#         print(f"  SNR vs template:    {snr_db:.2f} dB\n")

#         # 6) Plot time waveforms (filtered signal + template)
#         title = (f"{os.path.basename(path)} | corr={corr:.3f}, "
#                  f"SNR={snr_db:.2f} dB, template {args.template_bpm:.2f} BPM")
#         plot_time_with_template(t_u, d_bp, s, title, t_plot=args.t_plot)

#     plt.tight_layout()
#     plt.show()


# if __name__ == "__main__":
#     main()




















# #!/usr/bin/env python3
# import argparse
# import csv
# import math
# import os

# import numpy as np
# from scipy.signal import butter, filtfilt, periodogram
# import matplotlib.pyplot as plt


# # ------------------------------------------------------------
# # 1. CSV loading
# # ------------------------------------------------------------

# def load_csv(csv_path):
#     """
#     Load a distance CSV: expects columns 'elapsed_ms' and 'm_dist'.
#     Returns:
#         t_sec_raw: np.ndarray, time in seconds (not normalized)
#         d_raw: np.ndarray, distance in meters
#     """
#     with open(csv_path, "r", newline="") as f:
#         reader = csv.DictReader(f)
#         t_list = []
#         d_list = []
#         for row in reader:
#             try:
#                 t_ms = float(row["elapsed_ms"])
#                 d = row["m_dist"]
#                 if d == "" or d is None:
#                     continue
#                 d_val = float(d)
#                 t_list.append(t_ms * 1e-3)  # ms -> s
#                 d_list.append(d_val)
#             except (KeyError, ValueError):
#                 # skip malformed rows
#                 continue

#     t_sec_raw = np.asarray(t_list, dtype=np.float64)
#     d_raw = np.asarray(d_list, dtype=np.float64)

#     # Normalize time starting from zero
#     if len(t_sec_raw) > 0:
#         t_sec_raw = t_sec_raw - t_sec_raw[0]

#     return t_sec_raw, d_raw


# # ------------------------------------------------------------
# # 2. Resample to uniform grid
# # ------------------------------------------------------------

# def resample_uniform(t_raw, d_raw, fs_target=None):
#     """
#     Resample (t_raw, d_raw) to a uniform grid using linear interpolation.
#     If fs_target is None, use 1 / median(dt).
#     Returns:
#         t_u: uniform time [s]
#         d_u: resampled distance [m]
#         fs: sampling rate [Hz]
#     """
#     if len(t_raw) < 2:
#         raise ValueError("Not enough samples to resample.")

#     dt_raw = np.diff(t_raw)
#     dt_med = np.median(dt_raw)
#     if fs_target is None:
#         fs = 1.0 / dt_med
#     else:
#         fs = fs_target

#     duration = t_raw[-1] - t_raw[0]
#     n_samples = int(round(duration * fs)) + 1

#     t_u = np.linspace(t_raw[0], t_raw[-1], n_samples)
#     d_u = np.interp(t_u, t_raw, d_raw)

#     return t_u, d_u, fs


# # ------------------------------------------------------------
# # 3. Band-pass filtering + *amplitude rescaling*
# # ------------------------------------------------------------

# def butter_bandpass(lowcut, highcut, fs, order=4):
#     nyq = 0.5 * fs
#     low = lowcut / nyq
#     high = highcut / nyq
#     if low <= 0:
#         low = 1e-6
#     if high >= 1:
#         high = 0.999999
#     b, a = butter(order, [low, high], btype="band")
#     return b, a


# def bandpass_and_rescale(d_u, fs, fmin, fmax):
#     """
#     1. Detrend (remove mean).
#     2. Band-pass filter between [fmin, fmax].
#     3. Rescale filtered signal so its *peak-to-peak* matches
#        the original detrended signal's peak-to-peak.
#        (This approximately restores the original amplitude.)
#     Returns:
#         d_bp_rescaled: filtered + rescaled signal (meters, approx)
#     """
#     # Detrend
#     d0 = d_u - np.mean(d_u)

#     # Original robust peak-to-peak amplitude (ignore outliers with percentiles)
#     orig_min = np.percentile(d0, 5.0)
#     orig_max = np.percentile(d0, 95.0)
#     orig_ptp = orig_max - orig_min

#     # Band-pass filter
#     b, a = butter_bandpass(fmin, fmax, fs, order=4)
#     d_bp = filtfilt(b, a, d0)

#     # Filtered robust peak-to-peak
#     bp_min = np.percentile(d_bp, 5.0)
#     bp_max = np.percentile(d_bp, 95.0)
#     bp_ptp = bp_max - bp_min

#     # Rescale to match original ptp
#     if bp_ptp > 1e-12 and orig_ptp > 0:
#         scale = orig_ptp / bp_ptp
#     else:
#         scale = 1.0

#     d_bp_rescaled = d_bp * scale

#     return d_bp_rescaled


# # ------------------------------------------------------------
# # 4. BPM estimation from FFT (after band-pass)
# # ------------------------------------------------------------

# def estimate_bpm_from_fft(d_bp, fs, fmin, fmax):
#     """
#     Use FFT (periodogram) on band-passed signal to find dominant frequency.
#     Returns:
#         f_peak: peak frequency [Hz]
#         bpm: f_peak * 60 [bpm]
#     """
#     freqs, power = periodogram(d_bp, fs=fs, scaling="spectrum")

#     # Restrict to [fmin, fmax]
#     mask = (freqs >= fmin) & (freqs <= fmax)
#     if not np.any(mask):
#         return np.nan, np.nan

#     idx = np.argmax(power[mask])
#     f_peak = freqs[mask][idx]
#     bpm = f_peak * 60.0
#     return f_peak, bpm


# # ------------------------------------------------------------
# # 5. Build fixed sine template and compare
# # ------------------------------------------------------------

# def build_fixed_template(t_u, bpm_template, A_template):
#     """
#     Build a fixed amplitude sine template:
#         s(t) = A * sin(2π f t) where f = bpm/60
#     We do NOT adapt amplitude or bpm per file.
#     """
#     f = bpm_template / 60.0
#     # zero phase is fine because both template and signal are band-passed around same f
#     template = A_template * np.sin(2.0 * math.pi * f * t_u)
#     return template


# def compare_to_template(d_bp_rescaled, template):
#     """
#     Compare filtered+rescaled signal to fixed template:
#     - correlation coefficient
#     - SNR in dB with respect to template
#     - normalized RMSE (RMSE / template amplitude)
#     """
#     # Truncate to common length if needed
#     n = min(len(d_bp_rescaled), len(template))
#     x = d_bp_rescaled[:n]
#     s = template[:n]

#     # Mean-center both (for correlation), but do NOT rescale amplitudes
#     x0 = x - np.mean(x)
#     s0 = s - np.mean(s)

#     # Correlation coefficient
#     denom = np.std(x0) * np.std(s0)
#     if denom < 1e-12:
#         corr = np.nan
#     else:
#         corr = np.mean(x0 * s0) / denom

#     # Error signal & powers
#     error = x - s
#     signal_power = np.mean(s ** 2)
#     noise_power = np.mean(error ** 2)

#     if noise_power <= 0 or signal_power <= 0:
#         snr_db = np.nan
#     else:
#         snr_db = 10.0 * np.log10(signal_power / noise_power)

#     # Normalized RMSE: divide by template amplitude (max abs)
#     A = np.max(np.abs(s)) if np.max(np.abs(s)) > 0 else 1.0
#     rmse = math.sqrt(np.mean(error ** 2))
#     nrmse = rmse / A

#     return corr, snr_db, nrmse, x, s


# # ------------------------------------------------------------
# # 6. Plotting helper
# # ------------------------------------------------------------

# def plot_signal_and_template(csv_name, t_u, x, s, corr, snr_db, nrmse, bpm_est, bpm_template, t_plot):
#     """
#     Plot first t_plot seconds of filtered+rescaled signal vs fixed template.
#     """
#     if t_plot is not None:
#         mask = t_u <= t_plot
#         if not np.any(mask):
#             mask = slice(None)
#     else:
#         mask = slice(None)

#     t_plot_vals = t_u[mask]
#     x_plot = x[mask]
#     s_plot = s[mask]

#     plt.figure(figsize=(10, 4))
#     title = (
#         f"{csv_name} | BPM_est≈{bpm_est:.2f}, template={bpm_template:.2f} "
#         f"| corr={corr:.3f}, SNR={snr_db:.2f} dB, nRMSE={nrmse:.3f}"
#     )
#     plt.title(title)
#     plt.plot(t_plot_vals, x_plot, label="Filtered + rescaled distance (m)")
#     plt.plot(t_plot_vals, s_plot, "--", label="Fixed sine template (m)")
#     plt.xlabel("Time [s]")
#     plt.ylabel("Distance [m]")
#     plt.legend()
#     plt.tight_layout()


# # ------------------------------------------------------------
# # 7. Main
# # ------------------------------------------------------------

# def main():
#     parser = argparse.ArgumentParser(
#         description="Analyze distance CSVs: band-pass + rescale, BPM, fixed-template SNR."
#     )
#     parser.add_argument(
#         "--csv",
#         nargs="+",
#         required=True,
#         help="List of CSV files to analyze.",
#     )
#     parser.add_argument(
#         "--fmin",
#         type=float,
#         default=0.4,
#         help="Lower cutoff frequency for band-pass [Hz].",
#     )
#     parser.add_argument(
#         "--fmax",
#         type=float,
#         default=1.0,
#         help="Upper cutoff frequency for band-pass [Hz].",
#     )
#     parser.add_argument(
#         "--template-bpm",
#         type=float,
#         default=40.0,
#         help="BPM of the fixed sine template used for ALL files.",
#     )
#     parser.add_argument(
#         "--template-A",
#         type=float,
#         default=0.003,
#         help="Amplitude (meters) of the fixed sine template used for ALL files.",
#     )
#     parser.add_argument(
#         "--t-plot",
#         type=float,
#         default=30.0,
#         help="Duration in seconds to show in plots (per file).",
#     )
#     args = parser.parse_args()

#     print(
#         f"[INFO] Using fixed template for ALL files: "
#         f"BPM={args.template_bpm}, A={args.template_A}, "
#         f"band=[{args.fmin}, {args.fmax}] Hz"
#     )

#     for csv_path in args.csv:
#         print(f"\n=== {csv_path} ===")
#         if not os.path.isfile(csv_path):
#             print(f"  [ERROR] File not found.")
#             continue

#         # 1) Load
#         t_raw, d_raw = load_csv(csv_path)
#         if len(t_raw) < 10:
#             print("  [ERROR] Not enough samples.")
#             continue

#         duration = t_raw[-1] - t_raw[0]
#         print(
#             f"  Loaded {len(d_raw)} samples; duration ≈ {duration:.1f} s"
#         )

#         # 2) Resample
#         t_u, d_u, fs = resample_uniform(t_raw, d_raw, fs_target=None)
#         print(f"  Resampled: fs ≈ {fs:.2f} Hz, N = {len(d_u)}")

#         # 3) Band-pass + rescale
#         d_bp_rescaled = bandpass_and_rescale(d_u, fs, args.fmin, args.fmax)

#         # 4) BPM estimation
#         f_peak, bpm_est = estimate_bpm_from_fft(d_bp_rescaled, fs, args.fmin, args.fmax)
#         print(f"  BPM estimate from FFT: {bpm_est:.2f} (f_peak={f_peak:.4f} Hz)")

#         # 5) Fixed template
#         template = build_fixed_template(t_u, args.template_bpm, args.template_A)

#         # 6) Compare
#         corr, snr_db, nrmse, x, s = compare_to_template(d_bp_rescaled, template)
#         print(f"  Corr with template: {corr:.3f}")
#         print(f"  SNR vs template:    {snr_db:.2f} dB")
#         print(f"  Normalized RMSE:    {nrmse:.3f}")

#         # 7) Plot
#         csv_name = os.path.basename(csv_path)
#         plot_signal_and_template(csv_name, t_u, x, s, corr, snr_db, nrmse, bpm_est, args.template_bpm, args.t_plot)

#     plt.show()


# if __name__ == "__main__":
#     main()















#!/usr/bin/env python3
"""
Phase-aligned breathing analysis:

1. Load m_dist vs elapsed_ms from one or more CSV files.
2. Resample to a uniform time grid.
3. Band-pass filter around breathing band [fmin, fmax].
4. Rescale filtered signal back to original amplitude (same std as raw).
5. Estimate BPM from FFT (dominant frequency in band).
6. Build a **phase-aligned sine template** at that frequency
   using least-squares (best-fit sin+cos).
7. Compute:
   - Pearson correlation between signal and template
   - 1D SSIM between signal and template
8. Plot first T seconds overlay (filtered+rescaled vs template).
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy import signal


# ---------- basic IO ----------

def load_distance_csv(path):
    """Return elapsed time (s) and distance (m) from one CSV."""
    t_ms = []
    d_m = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        if "elapsed_ms" not in reader.fieldnames or "m_dist" not in reader.fieldnames:
            raise ValueError(f"{path}: header must contain 'elapsed_ms' and 'm_dist', got {reader.fieldnames}")
        for row in reader:
            try:
                t_ms.append(float(row["elapsed_ms"]))
                d_m.append(float(row["m_dist"]))
            except (ValueError, TypeError):
                continue
    t_ms = np.asarray(t_ms)
    d_m = np.asarray(d_m)
    # sort by time just in case
    idx = np.argsort(t_ms)
    return t_ms[idx] / 1000.0, d_m[idx]


def resample_uniform(t, x):
    """Resample x(t) to uniform grid, returning t_u, x_u, fs."""
    # estimate sampling interval from median dt
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


def filter_and_rescale(t_u, d_u, fmin, fmax):
    """
    Band-pass d_u and rescale so that
    std(filtered) == std(original zero-mean signal).
    Returns filtered_rescaled (same units as d_u) and fs.
    """
    # remove DC from original for a fair std comparison
    d_zero = d_u - np.mean(d_u)
    t_u, d_zero, fs = resample_uniform(t_u, d_zero)  # t_u is already uniform typically, but safe

    d_bp = bandpass_filter(d_zero, fs, fmin, fmax)
    std_orig = np.std(d_zero)
    std_filt = np.std(d_bp)
    if std_filt > 0:
        scale = std_orig / std_filt
    else:
        scale = 1.0
    d_bp_scaled = d_bp * scale
    return t_u, d_bp_scaled, fs


# ---------- BPM estimate from FFT ----------

def estimate_bpm_fft(x, fs, fmin, fmax):
    """Estimate BPM by FFT peak within [fmin, fmax]."""
    N = len(x)
    # real FFT
    freqs = np.fft.rfftfreq(N, d=1.0/fs)
    X = np.fft.rfft(x)
    power = np.abs(X) ** 2

    # restrict to band
    band = (freqs >= fmin) & (freqs <= fmax)
    if not np.any(band):
        return None, None

    idx_peak = np.argmax(power[band])
    f_peak = freqs[band][idx_peak]
    bpm = 60.0 * f_peak
    return f_peak, bpm


# ---------- phase-aligned template ----------

def build_phase_aligned_template(t_u, x, f0):
    """
    Given time grid t_u and signal x (already band-passed & rescaled),
    fit x(t) ≈ A * sin(2π f0 t + phi) via least-squares using
    sin & cos basis, returning template y and (A, phi).
    """
    t0 = t_u - t_u[0]  # start from zero
    s_sin = np.sin(2 * np.pi * f0 * t0)
    s_cos = np.cos(2 * np.pi * f0 * t0)

    # x ≈ a * sin + b * cos
    M = np.vstack([s_sin, s_cos]).T  # N x 2
    params, *_ = np.linalg.lstsq(M, x, rcond=None)
    a, b = params
    A = np.sqrt(a**2 + b**2)
    phi = np.arctan2(b, a)

    y = A * np.sin(2 * np.pi * f0 * t0 + phi)
    return y, A, phi


# ---------- correlation + SSIM ----------

def pearson_corr(x, y):
    """Pearson correlation for 1D arrays."""
    x = np.asarray(x)
    y = np.asarray(y)
    if x.size != y.size:
        raise ValueError("pearson_corr: x and y must have same length")
    if np.allclose(np.std(x), 0) or np.allclose(np.std(y), 0):
        return 0.0
    return np.corrcoef(x, y)[0, 1]


def ssim_1d(x, y):
    """
    Simple 1D SSIM implementation (single window over full signal).
    x, y are 1D arrays.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size != y.size:
        raise ValueError("ssim_1d: x and y must have same length")

    mu_x = x.mean()
    mu_y = y.mean()
    sigma_x2 = x.var()
    sigma_y2 = y.var()
    sigma_xy = ((x - mu_x) * (y - mu_y)).mean()

    # dynamic range for constants
    L = max(x.max() - x.min(), y.max() - y.min())
    if L == 0:
        L = 1.0
    C1 = (0.01 * L) ** 2
    C2 = (0.03 * L) ** 2

    num = (2 * mu_x * mu_y + C1) * (2 * sigma_xy + C2)
    den = (mu_x**2 + mu_y**2 + C1) * (sigma_x2 + sigma_y2 + C2)
    return float(num / den)


# ---------- main per-file pipeline ----------

def analyze_file(path, fmin, fmax, t_plot):
    t, d = load_distance_csv(path)
    duration = t[-1] - t[0]
    print(f"\n=== {path} ===")
    print(f"  Loaded {len(t)} samples; duration ≈ {duration:.1f} s")

    # band-pass + rescale
    t_u, d_bp_scaled, fs = filter_and_rescale(t, d, fmin, fmax)
    print(f"  Resampled fs ≈ {fs:.2f} Hz, N = {len(t_u)}")

    # BPM from FFT
    f0, bpm = estimate_bpm_fft(d_bp_scaled, fs, fmin, fmax)
    if f0 is None:
        print("  [WARN] No frequency in the given band.")
        return
    print(f"  BPM_est ≈ {bpm:.2f} (f_peak={f0:.4f} Hz)")

    # phase-aligned template with same frequency
    template, A, phi = build_phase_aligned_template(t_u, d_bp_scaled, f0)
    print(f"  Fitted template amplitude A ≈ {A:.4e} m, phase φ ≈ {phi:.3f} rad")

    # metrics
    corr = pearson_corr(d_bp_scaled, template)
    ssim_val = ssim_1d(d_bp_scaled, template)
    print(f"  Pearson corr: {corr:.3f}")
    print(f"  SSIM (1D):    {ssim_val:.3f}")

    # plot first t_plot seconds
    mask = t_u - t_u[0] <= t_plot
    t_plot_arr = t_u[mask] - t_u[0]
    sig_plot = d_bp_scaled[mask]
    tmpl_plot = template[mask]

    plt.figure(figsize=(10, 4))
    title = (f"{path} | BPM_est≈{bpm:.2f}, "
             f"corr={corr:.3f}, SSIM={ssim_val:.3f}")
    plt.title(title)
    plt.plot(t_plot_arr, sig_plot, label="Filtered + rescaled distance (m)")
    plt.plot(t_plot_arr, tmpl_plot, "--", label="Phase-aligned sine template (m)")
    plt.xlabel("Time [s]")
    plt.ylabel("Distance [m]")
    plt.legend(loc="upper right")
    plt.grid(True)


# ---------- CLI ----------

def main():
    parser = argparse.ArgumentParser(
        description="Phase-aligned breathing analysis: BPM + correlation + SSIM"
    )
    parser.add_argument(
        "--csv", nargs="+", required=True,
        help="One or more CSV files (distance_log_*.csv)"
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
        help="Seconds to plot from start (default 30)"
    )

    args = parser.parse_args()

    for path in args.csv:
        analyze_file(path, args.fmin, args.fmax, args.t_plot)

    plt.show()


if __name__ == "__main__":
    main()
