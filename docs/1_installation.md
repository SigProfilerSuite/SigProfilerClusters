# Installation


----------


This section will help you set up the necessary software and packages required to run SigProfilerClusters.

----------


## Prerequisites ##

- [Python][1] version >= 3.9
- [SigProfilerMatrixGenerator][2] with a downloaded reference genome (required to generate mutational matrices and simulations)
- [SigProfilerSimulator][3] (required to generate the background model used for IMD threshold calculation)
- Other dependencies are installed automatically during package installation

## Installation ##

SigProfilerClusters can be executed on any Windows/macOS/Unix system. First follow the [SigProfilerMatrixGenerator][2] guide for installing `Python` and `pip`. Next, follow the instructions below for the latest stable release or the current GitHub version.

### Installation with `pip` ###

Install the latest `SigProfilerClusters` PyPI version using `pip`:
```
$ pip install SigProfilerClusters
```

To upgrade an existing installation to the most recent version:
```
$ pip install SigProfilerClusters --upgrade
```

### Install specific GitHub Release ###

First, download the [zip file][4] or clone the GitHub repository:
```
$ git clone https://github.com/SigProfilerSuite/SigProfilerClusters.git
```

Next, enter the downloaded directory and install the package:
```
$ cd SigProfilerClusters
$ pip install .
```

## Download Reference Genome ##

SigProfilerClusters requires a reference genome to generate simulations and mutational matrices. Current reference genomes supported include GRCh37, GRCh38, mm9, and mm10. To install the reference genome/s, you need to use [SigProfilerMatrixGenerator][2].

The last PyPI [SigProfilerMatrixGenerator][2] version is installed with SigProfilerClusters by default. You can also install a specific version following the instructions in the [SigProfilerMatrixGenerator Wiki][2].

Once [SigProfilerMatrixGenerator][2] is installed, install your desired reference genome from the command line/terminal as follows.

### Installation from command line ###

```
$ SigProfilerMatrixGenerator install GRCh37
```

### Installation from Python terminal ###

``` python
$ python
>>> from SigProfilerMatrixGenerator import install as genInstall
>>> genInstall.install('GRCh37', rsync=False, bash=True)
```

If you have a firewall on your server, you may need to install `rsync` and use the `rsync=True` parameter. Similarly, if you do not have bash available, use `bash=False`.

In case you prefer to install a reference genome from a local copy:
``` python
$ python
>>> from SigProfilerMatrixGenerator import install as genInstall
>>> genInstall.install('GRCh37', offline_files_path='path/to/directory/containing/GRCh37.tar.gz')
```

  [1]: https://www.python.org/downloads
  [2]: https://sigprofilersuite.github.io/SigProfilerMatrixGenerator/
  [3]: https://osf.io/usxjz/wiki/home/
  [4]: https://github.com/SigProfilerSuite/SigProfilerClusters/releases
