# MultiOmicsIntegrator

Multi-omics integration pipeline combining simulated RNA-seq and
ATAC-seq data from 220 samples via a real MOFA-style joint
dimensionality reduction, discovering molecular subtypes validated
against ground truth.

> Built a multi-omics integration pipeline combining RNA-seq and
> ATAC-seq data across 220 simulated samples; MOFA-style joint factor
> analysis plus K-means clustering recovers all 8 true cancer subtypes
> with a perfect Adjusted Rand Index of 1.0 (100% sample-level
> best-match accuracy).

## What's real vs. synthetic — read this first

**Real:** `src/mofa_integration.py` implements genuine multi-view
joint dimensionality reduction — z-scoring each omics view
independently, concatenating them, and fitting real PCA across the
joint feature space — the same core idea as MOFA (Multi-Omics Factor
Analysis). The real R/Bioconductor `MOFA2` package is not installed in
this environment; `docs/mofa_reference.R` is an authored-but-not-
executed real MOFA2 script showing the equivalent R workflow, and
every result in this repo comes from the executed Python mirror, not
from that script. `src/subtype_clustering.py` runs genuine K-means and
validates against known true subtype labels via the Adjusted Rand
Index. `src/pathway_enrichment.py` runs a real hypergeometric
enrichment test.

**Synthetic:** there is no real GEO RNA-seq/ATAC-seq dataset in this
environment. `src/simulate_multiomics.py` generates 220 samples across
8 true subtypes, each subtype with its own latent factor loading onto
both the RNA and ATAC feature spaces.

## Results (real, from `run_all.py`)

```
220 samples, 60 RNA features, 60 ATAC features, 10 joint factors

Subtype recovery:
  8/8 true subtypes discovered by K-means on joint factors
  Adjusted Rand Index: 1.0
  Best-match sample accuracy: 100.0%

Factor 1 explains 15.82% of joint variance (top factor)
```

`tests/test_pipeline.py` (5 tests, all passing) checks cohort size
>=200, ARI >=0.5, exactly 8 discovered clusters, valid enrichment
p-values, and correct feature counts per omics layer.

## Honest limitation: pathway enrichment

`src/pathway_enrichment.py` tests Factor 1's top-loading genes against
a randomly generated 5-pathway gene-set catalog (there is no real
biological pathway database in this environment) — none of the 5
synthetic pathways reach significance (p<0.05) against this factor,
which is the real, honestly reported output, not a forced positive
result. The subtype-discovery result above does not depend on this
enrichment step.

## Usage

```bash
pip install -r requirements.txt
python3 run_all.py
python3 -m pytest tests/ -v
```

## Project structure

```
src/
  simulate_multiomics.py    # 220-sample RNA+ATAC simulation with 8 true subtypes
  mofa_integration.py       # MOFA-style joint PCA factor analysis
  subtype_clustering.py     # K-means subtype discovery + ARI validation
  pathway_enrichment.py     # hypergeometric enrichment on factor loadings
docs/mofa_reference.R       # authored-but-not-executed real MOFA2 R script
run_all.py
tests/test_pipeline.py
```

## Tech stack

Python 3.11, scikit-learn (PCA, KMeans), scipy (hypergeometric test),
pandas, numpy, pytest. (The resume names R; the real MOFA2 R package
workflow is authored in `docs/mofa_reference.R` but not executed since
R/Bioconductor is not available here — this Python implementation
mirrors the same joint-factor-analysis methodology.)

## Limitations

- No real GEO multi-omics dataset is used.
- The Python joint-PCA implementation is a simplified analogue of
  MOFA2's full probabilistic multi-view factor model (which also
  handles missing views and per-view noise models).
- The pathway gene-set catalog is randomly generated, not a real
  biological database (e.g. MSigDB/Reactome).
