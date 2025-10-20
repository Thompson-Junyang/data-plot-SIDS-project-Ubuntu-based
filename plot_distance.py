#!/usr/bin/env python3
"""
Plot elapsed time vs. m_dist from the AprilTag CSV.

Usage:
  python plot_distance.py --csv distance_log_YYYYMMDD_HHMMSS.csv \
                          --out distance_plot.png \
                          --ms   # (optional) keep x-axis in milliseconds (default is seconds)
"""

import argparse
import csv
import math
import os
import matplotlib.pyplot as plt  # install: pip install matplotlib

def load_csv(path):
    xs_elapsed_ms = []
    ys_m = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        # Expect columns: timestamp_ms, elapsed_ms, pix_dist, m_dist
        for row in reader:
            m_dist_str = (row.get("m_dist") or "").strip()
            if m_dist_str == "":
                continue
            try:
                elapsed_ms = float(row["elapsed_ms"])
                m_dist = float(m_dist_str)
            except (KeyError, ValueError):
                continue
            xs_elapsed_ms.append(elapsed_ms)
            ys_m.append(m_dist)
    return xs_elapsed_ms, ys_m

def main():
    ap = argparse.ArgumentParser(description="Plot elapsed time vs. m_dist from CSV.")
    ap.add_argument("--csv", required=True, help="Path to CSV file (from your detector).")
    ap.add_argument("--out", default=None, help="Output image path (e.g., plot.png). If omitted, will show only.")
    ap.add_argument("--ms", action="store_true", help="Use milliseconds on x-axis (default: seconds).")
    args = ap.parse_args()

    xs_ms, ys = load_csv(args.csv)
    if not xs_ms or not ys:
        raise SystemExit("No valid (elapsed_ms, m_dist) pairs found in CSV.")

    # Convert x-axis to seconds by default (cleaner to read)
    if args.ms:
        xs = xs_ms
        x_label = "Elapsed time (ms)"
    else:
        xs = [v / 1000.0 for v in xs_ms]
        x_label = "Elapsed time (s)"

    plt.figure(figsize=(10, 5))
    plt.plot(xs, ys, linewidth=1.0)
    plt.xlabel(x_label)
    plt.ylabel("Tag center distance (m)")
    plt.title("Elapsed time vs. m_dist")
    plt.grid(True)

    if args.out:
        plt.savefig(args.out, dpi=150, bbox_inches="tight")
        print(f"[INFO] Saved plot -> {os.path.abspath(args.out)}")
    else:
        plt.show()

if __name__ == "__main__":
    main()
