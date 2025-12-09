#!/usr/bin/env python3
import numpy as np
import csv
import argparse
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
from scipy.interpolate import interp1d


# ============================================================
# BAND-PASS FILTER (0.3–1.5 Hz)
# ============================================================

def bandpass_filter(x, fs, f_low=0.3, f_high=1.5, order=4):
    """
    Zero-phase Butterworth band-pass filter.
    x      : 1D signal
    fs     : sampling rate (Hz)
    f_low  : lower cutoff (Hz)
    f_high : upper cutoff (Hz)
    order  : filter order
    """
    nyq = 0.5 * fs
    low = f_low / nyq
    high = f_high / nyq
    b, a = butter(order, [low, high], btype="band")
    y = filtfilt(b, a, x)
    return y


# ============================================================
# RESAMPLING TO UNIFORM GRID
# ============================================================

def resample_uniform(t, x, fs_desired=31.25):
    """
    Resample irregular time-stamps to a uniform grid (required for FFT).
    """
    t0, t1 = t[0], t[-1]
    N = int((t1 - t0) * fs_desired)
    t_uniform = np.linspace(t0, t1, N)
    interp = interp1d(t, x, kind='linear', fill_value="extrapolate")
    x_u = interp(t_uniform)
    return t_uniform, x_u, fs_desired


# ============================================================
# FFT + DOMINANT FREQUENCY
# ============================================================

def compute_fft_peak(x, fs, fmin=0.4, fmax=1.0):
    N = len(x)
    X = np.fft.rfft(x)
    freqs = np.fft.rfftfreq(N, d=1/fs)
    power = np.abs(X)**2

    mask = (freqs >= fmin) & (freqs <= fmax)
    idx = np.argmax(power[mask])
    f_peak = freqs[mask][idx]
    return f_peak, freqs, power


# ============================================================
# TEMPLATE MATCHING
# ============================================================

def template_correlation(signal, fs, f_peak):
    t = np.arange(len(signal)) / fs

    sine_template = np.sin(2*np.pi*f_peak*t)
    triangle_template = 2*np.abs(2*((t*f_peak) % 1) - 1) - 1

    def corr(a, b):
        a = (a - np.mean(a)) / (np.std(a) + 1e-8)
        b = (b - np.mean(b)) / (np.std(b) + 1e-8)
        return np.mean(a * b)

    return corr(signal, sine_template), corr(signal, triangle_template)


# ============================================================
# SNR ESTIMATION
# ============================================================

def estimate_snr(freqs, power, f_peak, bandwidth=0.05):
    sig_mask = (freqs >= f_peak - bandwidth) & (freqs <= f_peak + bandwidth)
    noise_mask = ~sig_mask

    signal_power = np.sum(power[sig_mask])
    noise_power = np.sum(power[noise_mask])

    if noise_power <= 0:
        return np.inf

    snr_db = 10 * np.log10(signal_power / noise_power)
    return snr_db


# ============================================================
# PLOTTING
# ============================================================

def plot_results(t, raw, filtered, freqs, power, f_peak, bpm):
    fig, axs = plt.subplots(2, 1, figsize=(14, 10))

    # ---- RAW vs FILTERED ----
    axs[0].plot(t, raw, alpha=0.45, label="Raw signal", linewidth=0.8)
    axs[0].plot(t, filtered, alpha=0.90, label="Filtered (band-pass)", linewidth=1.2)
    axs[0].legend()
    axs[0].set_title("Raw vs Filtered Signal (Band-pass 0.3–1.5 Hz)")
    axs[0].set_xlabel("Time (s)")
    axs[0].set_ylabel("Distance (m)")

    # ---- FFT ----
    axs[1].plot(freqs, power, label="FFT Power Spectrum")
    axs[1].axvline(f_peak, color="red", linestyle="--",
                   label=f"Peak = {f_peak:.3f} Hz ({bpm:.2f} BPM)")
    axs[1].set_xlim(0, 3)
    axs[1].set_title("Power Spectrum")
    axs[1].set_xlabel("Frequency (Hz)")
    axs[1].set_ylabel("Power")
    axs[1].legend()

    plt.tight_layout()
    plt.show()


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=str, required=True)
    parser.add_argument("--fmin", type=float, default=0.4)
    parser.add_argument("--fmax", type=float, default=1.0)
    args = parser.parse_args()

    # --------------------------------------------------------
    # LOAD CSV
    # --------------------------------------------------------
    print("Loading CSV:", args.csv)
    t, d = [], []
    with open(args.csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            t.append(float(row["elapsed_ms"]) / 1000.0)
            d.append(float(row["m_dist"]))

    t = np.array(t)
    d = np.array(d)
    print(f"Loaded {len(d)} samples; duration ≈ {t[-1] - t[0]:.1f} s")

    # --------------------------------------------------------
    # RESAMPLE
    # --------------------------------------------------------
    t_u, d_u, fs = ressample = resample_uniform(t, d)
    d_u = d_u - np.mean(d_u)

    # --------------------------------------------------------
    # NEW: BAND-PASS FILTER
    # --------------------------------------------------------
    filtered = bandpass_filter(d_u, fs, f_low=0.3, f_high=1.5, order=4)

    # --------------------------------------------------------
    # FFT ON FILTERED SIGNAL
    # --------------------------------------------------------
    f_peak, freqs, power = compute_fft_peak(filtered, fs, args.fmin, args.fmax)
    bpm = f_peak * 60
    print(f"\n[FFT] Dominant frequency: {f_peak:.4f} Hz -> BPM ≈ {bpm:.2f}")

    # --------------------------------------------------------
    # TEMPLATE MATCHING
    # --------------------------------------------------------
    corr_sin, corr_tri = template_correlation(filtered, fs, f_peak)
    print(f"[Template] Correlation with sine:     {corr_sin:.3f}")
    print(f"[Template] Correlation with triangle: {corr_tri:.3f}")

    # --------------------------------------------------------
    # SNR
    # --------------------------------------------------------
    snr_db = estimate_snr(freqs, power, f_peak)
    print(f"[SNR] Estimated SNR: {snr_db:.2f} dB")

    # --------------------------------------------------------
    # PLOTS
    # --------------------------------------------------------
    plot_results(t_u, d_u, filtered, freqs, power, f_peak, bpm)


if __name__ == "__main__":
    main()
