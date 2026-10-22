"""Simulate paired RNA-seq (expression) and ATAC-seq (accessibility)
profiles for 220 tumor samples across 8 true molecular subtypes, each
subtype driven by its own latent factor affecting both omics layers --
the ground truth a multi-omics integration pipeline must recover."""
from __future__ import annotations

import numpy as np
import pandas as pd

N_SAMPLES = 220
N_SUBTYPES = 8
N_RNA_FEATURES = 60
N_ATAC_FEATURES = 60
SEED = 521


def simulate(n_samples: int = N_SAMPLES, n_subtypes: int = N_SUBTYPES, seed: int = SEED):
    rng = np.random.default_rng(seed)

    subtype = rng.integers(0, n_subtypes, n_samples)
    subtype_factor_loadings_rna = rng.normal(0, 2.0, (n_subtypes, N_RNA_FEATURES))
    subtype_factor_loadings_atac = rng.normal(0, 2.0, (n_subtypes, N_ATAC_FEATURES))

    rna = np.zeros((n_samples, N_RNA_FEATURES))
    atac = np.zeros((n_samples, N_ATAC_FEATURES))

    for i in range(n_samples):
        st = subtype[i]
        rna[i] = subtype_factor_loadings_rna[st] + rng.normal(0, 1.0, N_RNA_FEATURES)
        atac[i] = subtype_factor_loadings_atac[st] + rng.normal(0, 1.0, N_ATAC_FEATURES)

    rna_cols = [f"RNA_{i:03d}" for i in range(N_RNA_FEATURES)]
    atac_cols = [f"ATAC_{i:03d}" for i in range(N_ATAC_FEATURES)]

    rna_df = pd.DataFrame(rna, columns=rna_cols)
    atac_df = pd.DataFrame(atac, columns=atac_cols)
    sample_ids = [f"S-{i:04d}" for i in range(n_samples)]
    rna_df.insert(0, "sample_id", sample_ids)
    atac_df.insert(0, "sample_id", sample_ids)

    true_subtypes = pd.Series(subtype, index=sample_ids, name="true_subtype")
    return rna_df, atac_df, true_subtypes


if __name__ == "__main__":
    rna_df, atac_df, true_subtypes = simulate()
    print(rna_df.shape, atac_df.shape, true_subtypes.value_counts().to_dict())
