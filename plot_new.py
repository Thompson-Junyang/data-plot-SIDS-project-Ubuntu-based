#!/usr/bin/env python3

# python plot_new.py distance_log_WIN_20251102_16_49_00_Pro_20251103_003808.csv --ms 120000 --dedup

# -*- coding: utf-8 -*-
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser(description="Plot distance vs time for first 200s")
    ap.add_argument("csv", help="input CSV with columns: timestamp_ms,elapsed_ms,pix_dist,m_dist")
    ap.add_argument("--out", default="distance_vs_time_0_200s.png", help="output figure filename")
    ap.add_argument("--ms", type=int, default=200_000, help="range in milliseconds (default: 200000 = 200s)")
    ap.add_argument("--dedup", action="store_true", help="drop duplicate elapsed_ms (keep first)")
    args = ap.parse_args()

    # 读取
    df = pd.read_csv(args.csv)

    # 基本字段检查
    needed = {"elapsed_ms", "m_dist"}
    missing = needed - set(df.columns)
    if missing:
        raise SystemExit(f"Missing columns: {missing}. CSV must contain {needed}")

    # 仅取前 200s
    df = df[df["elapsed_ms"] <= args.ms].copy()

    # 可选：按 elapsed_ms 去重（你给的样本末尾有重复时间戳）
    if args.dedup:
        df = df.sort_values("elapsed_ms", kind="stable").drop_duplicates(subset=["elapsed_ms"], keep="first")

    # 毫秒→秒
    df["t_sec"] = df["elapsed_ms"] / 1000.0
    y = df["m_dist"].astype(float)

    # 作图
    plt.figure(figsize=(10, 4.5))
    plt.plot(df["t_sec"], y)
    plt.title("Distance vs Time (first 200 s)")
    plt.xlabel("Time (s)")
    plt.ylabel("Distance (m)")
    plt.tight_layout()
    plt.savefig(args.out, dpi=180)
    print(f"Saved: {args.out}  (points={len(df)})")

if __name__ == "__main__":
    main()
