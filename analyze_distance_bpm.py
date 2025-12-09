#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analyze distance-vs-time CSV from the AprilTag system, estimate BPM,
and compare to a reference sine wave (e.g., 40 or 41 BPM).

Behavior:
  - 使用 CSV 里的 elapsed_ms 和 m_dist
  - 只在波形的「上半波」找峰值
  - 用“最小峰间距”保证每个周期只算一个峰
"""

import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt


def load_csv(path):
    """
    Load CSV with columns:
        timestamp_ms, elapsed_ms, pix_dist, m_dist
    Return:
        t_s     : np.ndarray, time in seconds (from elapsed_ms)
        dist_m  : np.ndarray, distance in meters (m_dist)
    """
    t_list = []
    d_list = []

    # 尝试 utf-8 读取；如果失败可以手动改成 encoding="gbk"
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        print("Header:", header)

        if ("elapsed_ms" not in header) or ("m_dist" not in header):
            raise ValueError(f"CSV must contain 'elapsed_ms' and 'm_dist', got: {header}")

        for row in reader:
            if not row["m_dist"]:
                continue  # 跳过没有距离的行
            try:
                t_ms = float(row["elapsed_ms"])
                d_m = float(row["m_dist"])
            except ValueError:
                continue
            t_list.append(t_ms)
            d_list.append(d_m)

    if len(t_list) == 0:
        raise ValueError("No valid samples found in CSV (check columns & data).")

    t_s = np.asarray(t_list, dtype=float) / 1000.0   # ms -> s
    dist_m = np.asarray(d_list, dtype=float)

    return t_s, dist_m


def detect_upper_peaks(t_s, dist_m, max_bpm=80.0, percentile_thresh=75.0):
    """
    只找上半波的峰值（local maxima），并且加“最小峰间距”约束。

    Parameters
    ----------
    t_s : np.ndarray
        Time in seconds.
    dist_m : np.ndarray
        Distance in meters.
    max_bpm : float
        你认为最大可能的 BPM，用来换算“最小峰间距”。
        如果预期在 30~60 bpm，可以设成 80。
    percentile_thresh : float
        只找高于这个百分位数的峰。75 或 80 一般比较安全。

    Returns
    -------
    peak_times : np.ndarray
        峰值所对应的时间（秒）
    peak_values : np.ndarray
        峰值距离（米）
    """

    N = len(dist_m)
    if N < 3:
        return np.array([]), np.array([])

    # 估算采样间隔 dt（秒）
    dt = np.median(np.diff(t_s))
    print(f"Estimated dt ≈ {dt:.4f} s  (~{1.0/dt:.1f} Hz)")

    # 以 75 / 80 百分位作为「上半波」阈值，自动适配波形
    thr = np.percentile(dist_m, percentile_thresh)
    print(f"Upper-peak threshold (percentile {percentile_thresh}) ≈ {thr:.6f} m")

    # 初步：找“局部最大值”且在上半波
    candidate_idx = []
    for i in range(1, N - 1):
        if dist_m[i] > thr and dist_m[i] > dist_m[i - 1] and dist_m[i] >= dist_m[i + 1]:
            candidate_idx.append(i)

    if not candidate_idx:
        print("No candidates found above threshold, try lowering percentile_thresh.")
        return np.array([]), np.array([])

    candidate_idx = np.asarray(candidate_idx, dtype=int)

    # 最小峰间距：由 max_bpm 决定。
    # 最大 BPM -> 最小周期（秒） -> 最小样本间隔
    min_period_s = 60.0 / max_bpm       # e.g. max_bpm=80 -> 0.75 s
    min_dist_samples = max(int(min_period_s / dt), 1)
    print(f"Min peak distance: {min_dist_samples} samples (for max_bpm={max_bpm})")

    # 按顺序筛选：相邻峰至少隔 min_dist_samples
    filtered_idx = []
    last_idx = -10 ** 9
    for i in candidate_idx:
        if i - last_idx >= min_dist_samples:
            filtered_idx.append(i)
            last_idx = i

    filtered_idx = np.asarray(filtered_idx, dtype=int)

    peak_times = t_s[filtered_idx]
    peak_values = dist_m[filtered_idx]

    print(f"Detected {len(peak_times)} upper peaks.")
    return peak_times, peak_values


def estimate_bpm_from_peaks(peak_times):
    """由峰值时间序列估算 BPM。"""
    if len(peak_times) < 2:
        return float("nan")

    duration_s = peak_times[-1] - peak_times[0]
    n_cycles = len(peak_times)
    bpm = n_cycles * 60.0 / duration_s
    return bpm


def main():
    parser = argparse.ArgumentParser(
        description="Estimate BPM from distance CSV (upper peaks only) and compare to reference sine."
    )
    parser.add_argument("--csv", type=str, required=True, help="CSV file, e.g. distance_log_50cm.csv")
    parser.add_argument("--bpm_ref", type=float, default=40.0, help="Reference BPM for sine wave (e.g. 40 or 41)")
    parser.add_argument("--max_bpm", type=float, default=80.0, help="Expected upper BPM bound for min-peak-distance")
    parser.add_argument("--percentile", type=float, default=75.0,
                        help="Percentile threshold for upper peaks (default 75)")
    args = parser.parse_args()

    # 1. 读数据
    t_s, dist_m = load_csv(args.csv)
    print(f"Loaded {len(t_s)} samples from {args.csv}")
    duration = t_s[-1] - t_s[0]
    print(f"Duration ≈ {duration:.1f} s")

    # 2. 找上半波的峰值
    peak_t, peak_y = detect_upper_peaks(
        t_s, dist_m,
        max_bpm=args.max_bpm,
        percentile_thresh=args.percentile
    )

    bpm_est = estimate_bpm_from_peaks(peak_t)
    print(f"Estimated BPM (upper peaks): {bpm_est:.2f}")

    # 3. 生成参考正弦波（同一时间轴）
    bpm_ref = args.bpm_ref
    freq_ref = bpm_ref / 60.0  # Hz
    omega = 2.0 * np.pi * freq_ref

    # 用真实波形的均值和振幅来塑形正弦波（只是为了好看，可随意调）
    mean_y = dist_m.mean()
    amp = (np.percentile(dist_m, 95) - np.percentile(dist_m, 5)) / 4.0
    y_ref = mean_y + amp * np.sin(omega * (t_s - t_s[0]))

    # 4. 画图：真实数据 + 峰值 + 参考正弦
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(t_s, dist_m, label="Measured distance (m)", linewidth=0.8)
    if len(peak_t) > 0:
        ax.scatter(peak_t, peak_y, color="red", s=10, label="Detected upper peaks")

    ax.plot(t_s, y_ref, color="orange", alpha=0.7,
            label=f"Reference sine ({bpm_ref:.1f} BPM)")

    title = f"Distance vs Time (BPM ≈ {bpm_est:.1f})"
    ax.set_title(title)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Distance (m)")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
