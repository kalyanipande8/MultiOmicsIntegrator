"""Real MOFA-style multi-view joint dimensionality reduction: z-scores
each omics view independently, concatenates them, and fits PCA across
the joint feature space to extract shared latent factors -- the same
core idea as MOFA (Multi-Omics Factor Analysis), implemented in
scikit-learn since the R MOFA2 package is not available in this
environment (see docs/mofa_reference.R for the authored-but-not-
executed real MOFA2 reference script)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

N_FACTORS = 10


def joint_factor_analysis(rna_df: pd.DataFrame, atac_df: pd.DataFrame, n_factors: int = N_FACTORS):
    rna_features = [c for c in rna_df.columns if c != "sample_id"]
    atac_features = [c for c in atac_df.columns if c != "sample_id"]

    rna_scaled = StandardScaler().fit_transform(rna_df[rna_features])
    atac_scaled = StandardScaler().fit_transform(atac_df[atac_features])

    joint = np.hstack([rna_scaled, atac_scaled])
    pca = PCA(n_components=n_factors, random_state=0)
    factors = pca.fit_transform(joint)

    factor_df = pd.DataFrame(factors, index=rna_df["sample_id"],
                              columns=[f"Factor_{i+1}" for i in range(n_factors)])
    return factor_df, pca
