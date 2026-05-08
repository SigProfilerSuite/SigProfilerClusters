# SigProfilerClusters

![Logo](assets/images/SigProfilerClusters.png)

----------

**SigProfilerClusters** is a Python framework for analyzing the inter-mutational distances (IMD) between SNV-SNV and INDEL-INDEL mutations. The tool separates mutations into clustered and non-clustered groups on a sample-dependent basis, and subclassifies all SNVs into one of four categories of clustered event: i) doublet base substitutions (DBS); ii) multi-base substitutions (MBS); iii) omikli; and iv) kataegis. Indels are identified as clustered or non-clustered but are not further subclassified.

SigProfilerClusters builds a sample-specific background model using [SigProfilerSimulator][1] and applies an IMD-based threshold — corrected for regional mutation density — to distinguish mutations that are unlikely to have co-occurred by chance. This approach allows reliable detection of clustered mutational events across a wide range of tumor types and mutation burdens.

**SigProfilerClusters** makes use of [SigProfilerSimulator][1] and [SigProfilerMatrixGenerator][2], enabling seamless integration with other tools in the SigProfiler suite.

The SigProfilerClusters library is available on [GitHub](https://github.com/SigProfilerSuite/SigProfilerClusters) and [PyPI](https://pypi.org/project/SigProfilerClusters).

!!! note "Previously known as SigProfilerHotSpots"
    This tool was previously distributed under the name **SigProfilerHotSpots**. For all usage instructions, `SigProfilerClusters` and `SigProfilerHotSpots` may be used interchangeably when working with the older version of the tool.

----------

### Citation

Bergstrom EN, Kundu M, Tbeileh N, Alexandrov LB. Examining clustered somatic mutations with SigProfilerClusters. *Bioinformatics*. 2022;38(13):3470–3473. [https://doi.org/10.1093/bioinformatics/btac335](https://doi.org/10.1093/bioinformatics/btac335)

Bergstrom EN, Luebeck J, Petljak M, et al. Mapping clustered mutations in cancer reveals APOBEC3 mutagenesis of ecDNA. *Nature*. 2022;602:510–517. [https://doi.org/10.1038/s41586-022-04398-6](https://doi.org/10.1038/s41586-022-04398-6)

### License

This software is copyrighted © 2022 by Erik Bergstrom, Alexandrov Lab. SigProfilerClusters is distributed under a BSD 2-Clause License. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the copyright notice and disclaimer are retained. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTIES OF ANY KIND, EXPRESS OR IMPLIED.

### Contact

For questions, support requests, or bug reports, please contact the development team via GitHub [issues](https://github.com/SigProfilerSuite/SigProfilerClusters/issues) or by email:

- Erik Bergstrom: [ebergstr@health.ucsd.edu](mailto:ebergstr@health.ucsd.edu)
- Mousumy Kundu: [mkundu@health.ucsd.edu](mailto:mkundu@health.ucsd.edu)

  [1]: https://osf.io/usxjz/wiki/home/
  [2]: https://osf.io/s93d5/wiki/home/
