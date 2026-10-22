import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.simulate_multiomics import simulate
from src.mofa_integration import joint_factor_analysis
from src.subtype_clustering import discover_subtypes, validate_subtype_recovery
from src.pathway_enrichment import build_gene_sets, enrich_top_factor_genes


def test_cohort_has_200plus_samples():
    rna_df, atac_df, true_subtypes = simulate()
    assert len(rna_df) >= 200


def test_subtype_recovery_is_strong():
    rna_df, atac_df, true_subtypes = simulate()
    factor_df, pca = joint_factor_analysis(rna_df, atac_df)
    discovered = discover_subtypes(factor_df, n_clusters=true_subtypes.nunique())
    recovery = validate_subtype_recovery(discovered, true_subtypes)
    assert recovery["adjusted_rand_index"] >= 0.5


def test_discovers_8_subtypes():
    rna_df, atac_df, true_subtypes = simulate()
    factor_df, pca = joint_factor_analysis(rna_df, atac_df)
    discovered = discover_subtypes(factor_df, n_clusters=true_subtypes.nunique())
    assert discovered.nunique() == 8


def test_pathway_enrichment_produces_real_pvalues():
    rna_df, atac_df, true_subtypes = simulate()
    factor_df, pca = joint_factor_analysis(rna_df, atac_df)
    rna_features = [c for c in rna_df.columns if c != "sample_id"]
    gene_sets = build_gene_sets(rna_features)
    enrichment = enrich_top_factor_genes(pca, rna_features, factor_idx=0, gene_sets=gene_sets)
    assert (enrichment["p_value"] >= 0).all() and (enrichment["p_value"] <= 1).all()


def test_rna_and_atac_feature_counts():
    rna_df, atac_df, true_subtypes = simulate()
    assert len([c for c in rna_df.columns if c != "sample_id"]) == 60
    assert len([c for c in atac_df.columns if c != "sample_id"]) == 60
