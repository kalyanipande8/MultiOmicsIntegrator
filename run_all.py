import json
import os

from src.simulate_multiomics import simulate
from src.mofa_integration import joint_factor_analysis
from src.subtype_clustering import discover_subtypes, validate_subtype_recovery
from src.pathway_enrichment import build_gene_sets, enrich_top_factor_genes


def main():
    rna_df, atac_df, true_subtypes = simulate()
    factor_df, pca = joint_factor_analysis(rna_df, atac_df)

    discovered = discover_subtypes(factor_df, n_clusters=true_subtypes.nunique())
    subtype_recovery = validate_subtype_recovery(discovered, true_subtypes)

    rna_features = [c for c in rna_df.columns if c != "sample_id"]
    gene_sets = build_gene_sets(rna_features)
    enrichment = enrich_top_factor_genes(pca, rna_features, factor_idx=0, gene_sets=gene_sets)

    summary = {
        "n_samples": len(rna_df), "n_rna_features": len(rna_features),
        "n_atac_features": len([c for c in atac_df.columns if c != "sample_id"]),
        "n_factors": factor_df.shape[1], "subtype_recovery": subtype_recovery,
        "variance_explained_pct": [round(float(v) * 100, 2) for v in pca.explained_variance_ratio_],
        "factor1_pathway_enrichment": enrichment.to_dict(orient="records"),
    }

    os.makedirs("results", exist_ok=True)
    with open("results/summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
