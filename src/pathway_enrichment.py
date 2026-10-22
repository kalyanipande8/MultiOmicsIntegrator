"""Real hypergeometric pathway enrichment on the top-loading genes of
each MOFA-style factor, against a small simulated gene-set catalog."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

N_PATHWAYS = 5


def build_gene_sets(rna_features: list, seed: int = 521) -> dict:
    rng = np.random.default_rng(seed)
    gene_sets = {}
    for p in range(N_PATHWAYS):
        size = rng.integers(8, 16)
        gene_sets[f"Pathway_{p}"] = set(rng.choice(rna_features, size, replace=False))
    return gene_sets


def enrich_top_factor_genes(pca, rna_features: list, factor_idx: int, gene_sets: dict,
                             top_n: int = 15, alpha: float = 0.05) -> pd.DataFrame:
    loadings = pca.components_[factor_idx][: len(rna_features)]
    top_genes = set(pd.Series(np.abs(loadings), index=rna_features).sort_values(ascending=False).index[:top_n])

    results = []
    universe_size = len(rna_features)
    for pathway, genes in gene_sets.items():
        overlap = len(top_genes & genes)
        p_value = hypergeom.sf(overlap - 1, universe_size, len(genes), len(top_genes))
        results.append({"pathway": pathway, "overlap": overlap, "pathway_size": len(genes),
                         "p_value": p_value, "enriched": p_value < alpha})
    return pd.DataFrame(results)
