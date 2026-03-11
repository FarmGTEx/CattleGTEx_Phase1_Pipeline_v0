#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Multi-set SNP enrichment in functional regions (E1–E15), no plotting.
Efficient streaming version for large SNP lists.
"""

import sys, math, argparse
from collections import defaultdict
import pandas as pd


# ---------- 1. Load annotation ----------
def load_annotation(bed_file):
    ann = defaultdict(list)
    with open(bed_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, region = parts[0], parts[1], parts[2], parts[3]
            if not chrom.startswith("chr"):
                chrom = "chr" + str(chrom)
            try:
                s, e = int(start), int(end)
            except ValueError:
                continue
            ann[chrom].append((s, e, region))
    for chrom in ann:
        ann[chrom].sort(key=lambda x: x[0])
    return ann


# ---------- 2. Find region ----------
def find_region(chrom, pos, ann):
    if chrom not in ann:
        return None
    intervals = ann[chrom]
    lo, hi = 0, len(intervals) - 1
    idx = -1
    while lo <= hi:
        mid = (lo + hi) // 2
        s, e, _ = intervals[mid]
        if s <= pos:
            idx = mid
            lo = mid + 1
        else:
            hi = mid - 1
    if idx == -1:
        return None
    for i in range(idx, max(idx - 5, -1), -1):
        s, e, region = intervals[i]
        if s <= pos < e:
            return region
    return None


# ---------- 3. Stream SNPs ----------
def process_snp_file(snp_file, ann):
    region_counts = defaultdict(int)
    total = 0
    with open(snp_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("_")
            if len(parts) != 2:
                continue
            chrom = parts[0]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            total += 1
            region = find_region(chrom, pos, ann)
            if region is not None:
                region_counts[region] += 1
    return total, region_counts


# ---------- 4. Fisher Exact Test ----------
def fisher_exact_twosided(a, b, c, d):
    if a == 0 or b == 0 or c == 0 or d == 0:
        odds_ratio = (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5))
    else:
        odds_ratio = (a * d) / (b * c)
    n1, n2 = a + b, c + d
    k, n = a + c, n1 + n2
    max_n = n
    log_fact = [0.0] * (max_n + 1)
    for i in range(1, max_n + 1):
        log_fact[i] = log_fact[i - 1] + math.log(i)
    def log_comb(n_, k_):
        if k_ < 0 or k_ > n_:
            return float("-inf")
        return log_fact[n_] - log_fact[k_] - log_fact[n_ - k_]
    def hypergeo_prob(x):
        return math.exp(log_comb(k, x) + log_comb(n - k, n1 - x) - log_comb(n, n1))
    p_obs = hypergeo_prob(a)
    x_min, x_max = max(0, n1 - (n - k)), min(n1, k)
    p_val = sum(hypergeo_prob(x) for x in range(x_min, x_max + 1) if hypergeo_prob(x) <= p_obs + 1e-12)
    return odds_ratio, min(p_val, 1.0)


# ---------- 5. Enrichment ----------
def calculate_enrichment(name, total_set, counts_set, total_bg, counts_bg):
    regions = sorted(set(counts_bg.keys()) | set(counts_set.keys()))
    rows = []
    for r in regions:
        a = counts_set.get(r, 0)
        c = counts_bg.get(r, 0)
        b, d = total_set - a, total_bg - c
        if (a + b) == 0 or (c + d) == 0:
            odds, p = float("nan"), 1.0
        else:
            odds, p = fisher_exact_twosided(a, b, c, d)
        rows.append({
            "set": name, "region": r,
            "set_count": a, "bg_count": c,
            "odds_ratio": odds, "p_value": p
        })
    df = pd.DataFrame(rows)
    df = df.sort_values("p_value")
    df.to_csv(f"enrichment_{name}.csv", index=False)
    return df


# ---------- 6. Main ----------
def main():
    parser = argparse.ArgumentParser(description="Multi-set SNP enrichment in functional regions (no plotting).")
    parser.add_argument("--annotation", required=True, help="Functional annotation BED file (no header)")
    parser.add_argument("--background", required=True, help="Background SNP list file (chr_pos)")
    parser.add_argument("--sets", nargs="+", required=True,
                        help="SNP sets in format Label:File, e.g. eQTL:eQTL_SNP.txt GWAS:GWAS_SNP.txt")
    args = parser.parse_args()

    print(">>> Loading annotation...")
    ann = load_annotation(args.annotation)

    print(">>> Processing background SNPs...")
    total_bg, counts_bg = process_snp_file(args.background, ann)
    print(f"Background total: {total_bg}, mapped: {sum(counts_bg.values())}")

    all_results = []

    for item in args.sets:
        name, path = item.split(":", 1)
        print(f"\n>>> Processing set: {name}")
        total, counts = process_snp_file(path, ann)
        print(f"  Total SNPs: {total}, mapped: {sum(counts.values())}")
        df = calculate_enrichment(name, total, counts, total_bg, counts_bg)
        all_results.append(df)

    df_all = pd.concat(all_results)
    df_all.to_csv("region_enrichment_all.csv", index=False)
    print(">>> Saved combined results: region_enrichment_all.csv")
    print("Done.")

if __name__ == "__main__":
    main()
