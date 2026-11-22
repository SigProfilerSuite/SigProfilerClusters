# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.2.5] - 2026-11-22

### Changed

- test: add unit tests for findClustersOfClusters_noVAF

- test: add unit tests for findClustersOfClusters with VAF

- test: add unit tests for pullCCF and pullVaf

- test: cover density correction in findClustersOfClusters

- test: cover re-segmentation of mixed-VAF events

- test: add unit tests for first_run IMD cutoff and clustered split

- test: add unit tests for densityCorrection and getResults

- test: add unit tests for the IMD histogram plotting helpers

- test: add unit tests for the IMD calculation of samples and simulations

- test: add unit tests for eventProbability

- test: cover the chrom_based guard of analysis

- test: add unit tests for the SBS96 and ID83 catalog plots

- test: add unit tests for calculateSampleIMDs

- test: add unit tests for the ID rainfall plots

- test: cover BED ranges and chromosome lengths in rainfall

- docs: align the documented defaults and VAF formats with the code

- build: make version.py the single version source and fix the changelog date

- test: make the variant-caller integration test fail on any mismatch

- test: add frozen golden tests for the IMD cutoff and clustered calls


### Fixed

- fix(convertToVCF): skip empty tables instead of raising IndexError

- fix(hotspot): start the sample IMD cutoff scan at the first observed bin

- fix(distance): compute the IMD between distinct sites so same-position mutations are kept

- fix(probability): make eventProbability robust, reproducible and non-destructive

## [1.2.4] - 2026-11-21

### Changed

- docs: add MkDocs documentation site

- docs: add logo and update docs workflow

- GHA workflow updated and deploy constrained to master

- Logo in readme updated

- Readme has been reduced

- contact updated

- Link to SPS wiki has been corrected

- build: bump version to 1.2.3

- test: add unit tests for convertToVCF path and line helpers

- test: extend seeds file to one seed per simulator process

- test: add unit tests for convertFiles ID and SBS conversion

- test: add unit tests for CLI argument handling, dispatch and entry point

- test: add unit tests for processivity subclassification and refineIMD

- refactor: remove unimportable SigProfilerHotSpots module


### Fixed

- fix: use fork context instead of spawn to avoid forkserver crash and slowdown

- fix: find input files in per-run subfolders created by newer MatrixGenerator

- fix(hotspot): import Iterable from collections.abc in flatten

## [1.2.2] - 2025-07-21

### Added
- Integrated `pytest` to ensure correct handling of the `variant_caller` parameter.

### Changed
- Update the standard names for the `variant_caller` parameter
  - Changed `sanger` to `caveman`
  - `TCGA` and `standardVC` have been merged into a more flexible `standard` option.
  - `standard` is now default and parses VAF from 8th and 10th columns of VCF files (`VAF=` or `AF=`).

### Fixed
- Plotting Stability: Added error handling to skip samples that are invalid and proceed with valid ones.
- Resolved a potential index out-of-bounds error in the `plot_hist` function.

## [1.2.1] - 2025-04-02

### Changed
- Add centromere coordinates for mm39 genome.

## [1.2.0] - 2025-02-24

### Changed
- Updated dependencies: Now requires **Pandas >= 2.0.0**, **NumPy >= 2.0.0**, and **Python >= 3.9**.
- Dropped support for **Python 3.8**

## [1.1.3] - 2024-10-31

### Added
- Added support for Variant Allele Frequency (VAF) values from the following variant callers:
  - Sanger
  - TCGA
  - standardVC
  - Mutect2
