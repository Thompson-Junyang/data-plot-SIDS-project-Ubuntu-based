#!/usr/bin/env python3
"""
Analyze breathing-like distance signal from AprilTag CSV.

Functions:
- Load CSV (elapsed_ms, m_dist)
- Resample to uniform time grid
- FFT-based BPM estimation
- Template matching (sine & triangle)
- SNR estimation by sine fit
- Plots time series + spectrum

Usage:
  python analyze_fft_template_snr.py --csv distance_log_50cm.csv --fmin 0.4 --fmax 1.2
"""

import argparse
import csv
import math
import os

import numpy as np
import matplotlib.pyplot as plt


# ------------------------- CSV LOADING ------------------------- #

def load_csv(path):
    """
    Load elapsed_ms and m_dist from CSV.
    Assumes header has 'elapsed_ms' and 'm_dist' columns.
    Returns:
        t_s : np.array of shape (N,) - time in seconds (elapsed)
        d_m : np.array of shape (N,) - distance in meters
    """
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)

        # Find column indices
        header_lower = [h.strip().lower() for h in header]
        try:
            idx_t = header_lower.index("elapsed_ms")
            idx_d = header_lower.index("m_dist")
        except ValueError:
            raise RuntimeError(
                f"Header must contain 'elapsed_ms' and 'm_dist'. Got: {header}"
            )

        t_list = []
        d_list = []
        for row in reader:
            if not row or len(row) <= max(idx_t, idx_d):
                continue
            t_str = row[idx_t].strip()
            d_str = row[idx_d].strip()
            if d_str == "" or d_str.lower() == "nan":
                continue
            try:
                t_ms = float(t_str)
                d = float(d_str)
            except ValueError:
                # skip malformed line
                continue
            t_list.append(t_ms * 1e-3)  # ms -> s
            d_list.append(d)

    if len(t_list) < 10:
        raise RuntimeError("Not enough valid samples in CSV.")

    t = np.asarray(t_list, dtype=np.float64)
    d = np.asarray(d_list, dtype=np.float64)
    return t, d


# ------------------------- RESAMPLING ------------------------- #

def resample_uniform(t, d):
    """
    Resample (t, d) to a uniform time grid using linear interpolation.
    Returns:
        t_u   : uniform time
        d_u   : resampled distance
        fs    : sampling frequency (Hz)
    """
    # Estimate typical dt from differences
    dt = np.median(np.diff(t))
    if dt <= 0:
        raise RuntimeError("Non-increasing timestamps in data.")
    fs = 1.0 / dt

    t_start = t[0]
    t_end = t[-1]
    t_u = np.arange(t_start, t_end, dt)
    d_u = np.interp(t_u, t, d)

    return t_u, d_u, fs


# ------------------------- FFT-BASED BPM ------------------------- #

def estimate_bpm_fft(t, d, fmin=0.4, fmax=1.2):
    """
    Estimate breathing frequency via FFT on detrended signal.
    Args:
        t    : uniform time (s)
        d    : distance (m)
        fmin : min frequency (Hz) to search, default ~24 bpm
        fmax : max frequency (Hz) to search, default ~72 bpm
    Returns:
        f_peak   : dominant frequency in that band (Hz)
        bpm      : f_peak * 60
        freqs    : array of FFT frequencies
        power    : power spectrum (magnitude^2)
    """
    dt = np.median(np.diff(t))
    y = d - np.mean(d)  # remove DC

    # Window to reduce leakage
    win = np.hanning(len(y))
    y_win = y * win

    spec = np.fft.rfft(y_win)
    freqs = np.fft.rfftfreq(len(y_win), d=dt)
    power = np.abs(spec) ** 2

    # Restrict to breathing band
    band = (freqs >= fmin) & (freqs <= fmax)
    if not np.any(band):
        raise RuntimeError("No frequencies in specified band. Check fmin/fmax.")

    idx_band = np.where(band)[0]
    idx_peak_in_band = idx_band[np.argmax(power[band])]
    f_peak = freqs[idx_peak_in_band]
    bpm = f_peak * 60.0

    return f_peak, bpm, freqs, power


# ------------------------- TEMPLATE MATCHING ------------------------- #

def triangle_wave(f, t):
    """
    Generate a normalized triangle wave at frequency f (Hz)
    using arcsin(sin()) identity. Range approx [-1, 1].
    """
    return (2.0 / np.pi) * np.arcsin(np.sin(2.0 * np.pi * f * t))


def normalized_correlation(x, y):
    """
    Pearson correlation between two 1D signals, ignoring DC.
    """
    x = x - np.mean(x)
    y = y - np.mean(y)
    denom = np.linalg.norm(x) * np.linalg.norm(y)
    if denom == 0:
        return 0.0
    return float(np.dot(x, y) / denom)


def template_matching(t, d, f_peak):
    """
    Compare detrended distance with sine and triangle templates
    at the same dominant frequency.
    Returns:
        corr_sin : correlation with sine template
        corr_tri : correlation with triangle template
    """
    y = d - np.mean(d)

    # Generate templates with unit amplitude
    sin_ref = np.sin(2.0 * np.pi * f_peak * t)
    tri_ref = triangle_wave(f_peak, t)

    corr_sin = normalized_correlation(y, sin_ref)
    corr_tri = normalized_correlation(y, tri_ref)

    return corr_sin, corr_tri


# ------------------------- SNR ESTIMATION ------------------------- #

def estimate_snr(t, d, f_peak):
    """
    Estimate SNR by fitting a single-frequency sine:
        d(t) ≈ A*sin(2πft) + B*cos(2πft)  (+ DC already removed)
    Then:
        signal_hat = fitted sine
        noise = residual
        SNR = 20 * log10(rms(signal_hat) / rms(noise))

    Returns:
        snr_db       : SNR in dB
        signal_hat   : fitted periodic component
    """
    y = d - np.mean(d)

    # Design matrix for sin & cos at f_peak
    sin_col = np.sin(2.0 * np.pi * f_peak * t)
    cos_col = np.cos(2.0 * np.pi * f_peak * t)
    X = np.column_stack([sin_col, cos_col])

    # Least-squares fit
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    signal_hat = X @ beta

    noise = y - signal_hat
    rms_signal = math.sqrt(np.mean(signal_hat ** 2))
    rms_noise = math.sqrt(np.mean(noise ** 2))

    if rms_noise == 0:
        snr_db = float("inf")
    else:
        snr_db = 20.0 * math.log10(rms_signal / rms_noise)

    return snr_db, signal_hat


# ------------------------- PLOTTING ------------------------- #

def plot_results(t, d, signal_hat, freqs, power, f_peak, bpm):
    """
    Plot:
      1) Distance vs time (full) with fitted sine overlaid
      2) Power spectrum with breathing peak marked
    """
    # Time-domain
    plt.figure(figsize=(12, 4))
    plt.plot(t, d, label="Measured distance (m)", linewidth=1)
    plt.plot(t, signal_hat + np.mean(d), label="Fitted sine (shifted)", linewidth=2)
    plt.xlabel("Time (s)")
    plt.ylabel("Distance (m)")
    plt.title(f"Distance vs Time (BPM ≈ {bpm:.1f})")
    plt.legend()
    plt.tight_layout()

    # Frequency-domain
    plt.figure(figsize=(8, 4))
    plt.plot(freqs, power)
    plt.axvline(f_peak, color="r", linestyle="--", label=f"Peak {f_peak:.3f} Hz")
    plt.xlim(0, 3.0)  # show up to 3 Hz
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power")
    plt.title("Power Spectrum")
    plt.legend()
    plt.tight_layout()

    plt.show()


# ------------------------- MAIN ------------------------- #

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, help="Path to distance_log_*.csv")
    parser.add_argument("--fmin", type=float, default=0.4,
                        help="Min frequency to search (Hz), default 0.4 ~ 24 bpm")
    parser.add_argument("--fmax", type=float, default=1.2,
                        help="Max frequency to search (Hz), default 1.2 ~ 72 bpm")
    args = parser.parse_args()

    print(f"Loading CSV: {args.csv}")
    t, d = load_csv(args.csv)
    print(f"Loaded {len(t)} samples; duration ≈ {t[-1] - t[0]:.1f} s")

    t_u, d_u, fs = resample_uniform(t, d)
    print(f"Resampled to uniform grid: fs ≈ {fs:.2f} Hz, N = {len(t_u)}")

    f_peak, bpm, freqs, power = estimate_bpm_fft(t_u, d_u, args.fmin, args.fmax)
    print(f"\n[FFT] Dominant frequency: {f_peak:.4f} Hz  ->  BPM ≈ {bpm:.2f}")

    corr_sin, corr_tri = template_matching(t_u, d_u, f_peak)
    print(f"[Template] Correlation with sine:     {corr_sin:.3f}")
    print(f"[Template] Correlation with triangle: {corr_tri:.3f}")

    snr_db, signal_hat = estimate_snr(t_u, d_u, f_peak)
    print(f"[SNR] Estimated SNR: {snr_db:.2f} dB")

    # Plot time-domain and spectrum
    plot_results(t_u, d_u, signal_hat, freqs, power, f_peak, bpm)


if __name__ == "__main__":
    main()
