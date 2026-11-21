# Workflow


----------


This section describes the methodology used by SigProfilerClusters to identify and classify clustered mutations. The tool follows a four-stage process that distinguishes genuinely clustered mutations from those distributed randomly across the genome.

!!! note "Minimum mutational burden"
    Reliable analysis requires a mutational burden greater than **0.005 mutations per megabase**. Samples below this threshold may not produce meaningful IMD-based clustering results.

----------

## Stage 1: Background Model Generation ##

A per-sample background model is created using [SigProfilerSimulator][1]. This tool randomizes each mutation across the genome while preserving:

- The desired trinucleotide sequence context
- Transcriptional strand bias ratios
- The mutational burden on each chromosome

A minimum of **100 simulations** is recommended to ensure a robust background distribution.

## Stage 2: IMD Threshold Calculation ##

SigProfilerClusters compares the distribution of real inter-mutational distances (IMD) with the simulated background to compute a **sample-specific IMD threshold**. The threshold is set at the point where at least **90% of mutations falling within it are unlikely to have occurred by chance**, i.e., their spacing cannot be explained by the background model.

This approach ensures that the clustering criterion adapts to the mutation burden and genomic context of each individual sample.

## Stage 3: Regional Correction ##

The global IMD threshold is further refined by applying a **localized correction across 1-megabase genomic windows**. This step accounts for the uneven distribution of mutations throughout the genome — regions with higher mutation density will have a proportionally adjusted threshold — preventing over- or under-calling of clustered events in heterogeneous mutational landscapes.

## Stage 4: Mutation Classification ##

After thresholding and correction, mutations are assigned to one of the following categories based on their inter-mutational spacing and, when available, variant allele frequencies (VAF):

| Class | Category | Description |
|-------|----------|-------------|
| Non-clustered | — | Mutations whose spacing is consistent with random expectation |
| Clustered | **DBS** (Doublet Base Substitution) | Exactly 2 adjacent SNVs on the same chromosomal position |
| Clustered | **MBS** (Multi-Base Substitution) | 3 or more adjacent SNVs on consecutive chromosomal positions |
| Clustered | **Omikli** | 2–3 clustered SNVs with IMD pattern consistent with a localized event |
| Clustered | **Kataegis** | 4 or more clustered SNVs with IMD pattern consistent with a regional hypermutation event |
| Clustered | **Other** | Clustered mutations with inconsistent allele frequencies; only identified when VAF/CCF data is available |

Indels are classified as clustered or non-clustered but are not further subclassified into the categories above.

  [1]: https://sigprofilersuite.github.io/SigProfilerSimulator/
