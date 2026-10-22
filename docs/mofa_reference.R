# Authored-but-not-executed reference: real MOFA2 (R/Bioconductor)
# workflow for multi-omics factor analysis. Not run in this environment
# (no R/MOFA2 installation available here); src/mofa_integration.py
# implements the same joint-factor-analysis methodology in Python as
# the executed mirror, and all results in this repo come from that
# Python implementation, not from this script.

library(MOFA2)

# Build the MultiAssayExperiment-style input: one matrix per view
# (samples as columns, features as rows), with a shared sample order.
rna_matrix  <- t(as.matrix(read.csv("data/rna_features.csv", row.names = 1)))
atac_matrix <- t(as.matrix(read.csv("data/atac_features.csv", row.names = 1)))

mofa_data <- list(RNA = rna_matrix, ATAC = atac_matrix)

mofa_object <- create_mofa(mofa_data)

data_opts  <- get_default_data_options(mofa_object)
model_opts <- get_default_model_options(mofa_object)
model_opts$num_factors <- 10
train_opts <- get_default_training_options(mofa_object)

mofa_object <- prepare_mofa(mofa_object, data_options = data_opts,
                            model_options = model_opts,
                            training_options = train_opts)

mofa_trained <- run_mofa(mofa_object, outfile = "results/mofa_model.hdf5")

factors <- get_factors(mofa_trained, factors = "all")[[1]]
write.csv(factors, "results/mofa_factors.csv")

# Downstream: cluster samples on `factors` (e.g. kmeans) to discover
# subtypes, and inspect get_weights(mofa_trained) for top-loading
# genes per factor for pathway enrichment -- the same two steps
# src/subtype_clustering.py and src/pathway_enrichment.py perform on
# the Python-derived factors.
