import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
import os
os.chdir(Path(__file__).parent.parent)

import numpy as np
import pandas as pd
from dge_analysis import benjamini_hochberg, run_dge


def test_bh_correction_never_decreases_order():
    pvals = np.array([0.001, 0.2, 0.03, 0.5, 0.0001])
    padj = benjamini_hochberg(pvals)
    assert np.all(padj >= pvals)  # BH-adjusted p is always >= raw p
    assert np.all(padj <= 1.0)


def test_bh_correction_monotonic_with_rank():
    pvals = np.sort(np.array([0.001, 0.01, 0.02, 0.03, 0.04, 0.05]))
    padj = benjamini_hochberg(pvals)
    # adjusted p-values should be non-decreasing when sorted by raw p-value
    assert np.all(np.diff(padj) >= -1e-10)


def test_run_dge_recovers_known_genes():
    res_df, log_cpm, coldata = run_dge()
    truth = pd.read_csv("data/ground_truth.csv")
    curated = truth[truth.is_curated_real_gene]
    merged = res_df.merge(curated, on="gene")
    recovery_rate = merged.significant.sum() / len(curated)
    # pipeline sanity check: most strongly-designed curated genes
    # should be recovered as significant - this validates the
    # pipeline logic, not real biology
    assert recovery_rate > 0.7


def test_direction_matches_ground_truth():
    res_df, log_cpm, coldata = run_dge()
    truth = pd.read_csv("data/ground_truth.csv")
    merged = res_df.merge(truth[truth.is_curated_real_gene], on="gene")
    sig = merged[merged.significant]
    # sign of detected log2FC should match the sign it was generated with
    same_sign = np.sign(sig.log2FoldChange) == np.sign(sig.true_log2fc)
    assert same_sign.mean() > 0.95


if __name__ == "__main__":
    tests = [v for k, v in list(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS: {t.__name__}")
    print(f"\n{len(tests)} tests passed.")
