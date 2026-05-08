# Quick Start Example


----------


This section provides a minimal example to get started with SigProfilerClusters. The following example uses a breast cancer sample (`BRCA`) or a melanoma sample (`MELA`) and demonstrates the complete workflow: generating simulations with [SigProfilerSimulator][2] and detecting clustered mutations with SigProfilerClusters.

----------

## Prerequisites ##

This tutorial requires that you have completed all steps in the [installation guide][1], specifically:

- Installed SigProfilerClusters
- Installed [SigProfilerSimulator][2]
- Downloaded the **GRCh37** reference genome using [SigProfilerMatrixGenerator][3]

## Downloading input example data ##

Download one of the example VCF files below and place it in a new project directory (e.g., `path/to/data/`):

- [BRCA_example_subs.vcf][5] — breast cancer substitutions
- [BRCA_example_indels.vcf][7] — breast cancer indels
- [MELA_example.vcf][6] — melanoma substitutions

These samples and their expected outputs are also available under the `examples/` directory within the [GitHub repository](https://github.com/SigProfilerSuite/SigProfilerClusters).

!!! note
    The `simulations/` folder is not included in the GitHub repository to reduce its memory footprint.

## Running SigProfilerClusters ##

First, start a Python interactive shell and import [SigProfilerMatrixGenerator][3], [SigProfilerSimulator][2], and SigProfilerClusters.

``` python
$ python
>>> from SigProfilerMatrixGenerator import install as genInstall
>>> from SigProfilerSimulator import SigProfilerSimulator as sigSim
>>> from SigProfilerClusters import SigProfilerClusters as hp
```

### Step 1: Run SigProfilerSimulator ###

Generate a background model by running at least 100 simulations of inter-mutational distances for your data. **Note**: Update `"path/to/data"` with the actual path to the directory containing the VCF file.

``` python
>>> sigSim.SigProfilerSimulator("BRCA", "path/to/data", "GRCh37", contexts=["288"], chrom_based=True, simulations=100)
```

### Step 2: Run SigProfilerClusters ###

Cluster mutations based on the simulated background distribution. **Note**: Update `"path/to/data"` with the actual path to your project directory.

``` python
>>> hp.analysis("BRCA", "GRCh37", "96", ["288"], "path/to/data",
                analysis="all", sortSims=True, subClassify=True,
                correction=True, calculateIMD=True, max_cpu=4,
                variant_caller="standard")
```

After SigProfilerClusters has finished running, the output will be organized under `path/to/data/output/`. Partitioned mutations are placed under `output/clustered/` and `output/nonClustered/`, and visualizations are found under `output/plots/`. To learn more about all output files, please refer to the [Using the Tool - Output][4] section.

## Additional Information ##

In the above example, unspecified parameters use their default values. All function arguments and their types are described in detail in the [Using the Tool - Input][8] section. To learn more about the generated output files, refer to [Using the Tool - Output][4].

  [1]: https://sigprofilersuite.github.io/SigProfilerClusters/1_installation.html
  [2]: https://osf.io/usxjz/wiki/home/
  [3]: https://sigprofilersuite.github.io/SigProfilerMatrixGenerator/
  [4]: https://sigprofilersuite.github.io/SigProfilerClusters/5_using_the_tool_output.html
  [5]: https://osf.io/gy2bh/
  [6]: https://osf.io/a8bq5/
  [7]: https://osf.io/tyw3u/
  [8]: https://sigprofilersuite.github.io/SigProfilerClusters/4_using_the_tool_input.html
