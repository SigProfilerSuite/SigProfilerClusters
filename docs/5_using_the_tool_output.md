# Using SigProfilerClusters - Output


----------


This section describes the output files and directories produced by SigProfilerClusters. All results are written under `[project_path]/output/`, organized into several subdirectories.

----------

## Output Overview ##

After a successful run, the output directory contains the following subdirectories:

```
[project_path]/output/
├── clustered/
├── nonClustered/
├── plots/
├── simulations/
└── vcf_files_corrected/
```

All errors and progress checkpoints are also saved as log files in the project directory:

- `SigProfilerClusters_[project]_[genome].err` — error log
- `SigProfilerClusters_[project]_[genome].out` — progress log

## `clustered/` Directory ##

Contains the partitioned clustered mutations, organized into subclass subdirectories. Within each subclass folder, one VCF file is saved per sample. The subclasses correspond to the categories described in the [Workflow][1] section:

| Subclass | Description |
|----------|-------------|
| `DBS/` | Doublet base substitutions (2 adjacent SNVs) |
| `MBS/` | Multi-base substitutions (3+ adjacent SNVs) |
| `omikli/` | Localized clustered events (2–3 SNVs) |
| `kataegis/` | Regional hypermutation events (4+ SNVs) |
| `other/` | Clustered events with inconsistent allele frequencies (only present when VAF/CCF data is available) |

When `probability=True`, an extra column is appended to each VCF file containing the probability of observing that clustered event within its local genomic region.

## `nonClustered/` Directory ##

Contains VCF files with the non-clustered mutations for each sample. Substitution and indel mutations are stored in separate files, one per sample.

## `plots/` Directory ##

Contains two types of visualizations generated per sample:

### IMD plots ###
Show the mutational spectra alongside the distribution of inter-mutational distances (IMD) for both the real data and the simulated background. These plots allow visual inspection of the clustering signal relative to the background model.

### Rainfall plots ###
Display the IMD distribution across genomic coordinates for each sample, with the sample-dependent IMD threshold overlaid. Clustered events are color-coded by subclass, providing a genome-wide view of mutational clustering patterns.

## `simulations/` Directory ##

Contains three subdirectories with intermediate simulation data used for threshold calculation:

| Subdirectory | Contents |
|--------------|----------|
| `imd_data/` | IMD values across all simulations |
| `original_imds/` | IMD values for the original (real) mutations |
| `sorted_sims/` | Simulation files sorted alphanumerically, used during threshold computation |

## `vcf_files_corrected/` Directory ##

Contains VCF files for both clustered and non-clustered mutations after the **localized IMD correction** step (Stage 3 of the [workflow][1]). These files reflect the final mutation assignments incorporating regional mutation density adjustments.

## Mutation Class Summary ##

Output files use the following class labels to describe the nature of each clustered event:

| Class | Categories included | Description |
|-------|---------------------|-------------|
| Class 1 | DBS, MBS, Omikli | Small-scale clustered events involving 2–3 mutations |
| Class 2 | Kataegis | Larger hypermutation events; further subdivided by processivity |
| Class 3 | Other clustered | Remaining clustered mutations not fitting Class 1 or Class 2 criteria |

Each output VCF row includes the mutation coordinates, alleles, IMD value, and allele frequency, enabling detailed downstream analysis of individual clustered events.

  [1]: https://sigprofilersuite.github.io/SigProfilerClusters/3_workflow.html
