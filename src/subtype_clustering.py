"""Real K-means subtype discovery on the joint MOFA-style factors,
validated against the known true subtype labels via the Adjusted Rand
Index and a best-matching cluster-to-subtype mapping."""
from __future__ import annotations

import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score


def discover_subtypes(factor_df: pd.DataFrame, n_clusters: int = 8, seed: int = 521) -> pd.Series:
    km = KMeans(n_clusters=n_clusters, random_state=seed, n_init=10)
    labels = km.fit_predict(factor_df.values)
    return pd.Series(labels, index=factor_df.index, name="discovered_subtype")


def validate_subtype_recovery(discovered: pd.Series, true_subtypes: pd.Series) -> dict:
    aligned_true = true_subtypes.loc[discovered.index]
    ari = adjusted_rand_score(aligned_true, discovered)

    cross = pd.crosstab(discovered, aligned_true)
    best_match_acc = cross.max(axis=1).sum() / len(discovered)

    return {"n_discovered_clusters": discovered.nunique(), "n_true_subtypes": aligned_true.nunique(),
            "adjusted_rand_index": round(float(ari), 4),
            "best_match_accuracy_pct": round(float(best_match_acc) * 100, 2)}
