import pytest
import argparse
from SigProfilerClusters.controllers.cli_controller import (
    parse_arguments_clusters,
    str2bool,
    str2list
)


def test_argument_parsing():
    args = parse_arguments_clusters(
        [
            "dummy_project",
            "GRCh38",
            "96",
            "96",
            "path/to/files",
            "--output_type", "all",
            "--analysis", "all",
            "--interdistance", "96",
            "--exome", "True",
            "--clustering_vaf", "False",
            "--extraction", "False",
            "--correction","True",
            "--startProcess", "1",
            "--endProcess", "25",
            "--totalIterations","1000",
            "--calculateIMD","True",
            "--chrom_based","False",
            "--max_cpu", "1",
            "--subClassify","False",
            "--variant_caller","standard",
            "--includedVAFs","True",
            "--includedCCFs","False",
            "--windowSize","1000000",
            "--bedRanges","None",
            "--plotIMDfigure","True",
            "--plotRainfall","True",
            "--probability","False"
            ],
        "Test argument parsing",
    )

    assert args.project == "dummy_project"
    assert args.genome == "GRCh38"
    assert args.contexts == "96"
    assert args.simContext == ["96"]
    assert args.input_path == "path/to/files"
    assert args.output_type == "all"
    assert args.analysis == "all"
    assert args.interdistance == "96"
    assert args.exome == True
    assert args.clustering_vaf == False
    assert args.extraction == False
    assert args.correction == True
    assert args.startProcess == 1
    assert args.endProcess == 25
    assert args.totalIterations == 1000
    assert args.calculateIMD == True
    assert args.chrom_based == False
    assert args.max_cpu == 1
    assert args.subClassify == False
    assert args.variant_caller == "standard"
    assert args.includedVAFs == True
    assert args.includedCCFs == False
    assert args.windowSize == 1000000
    assert args.bedRanges == "None"
    assert args.plotIMDfigure == True
    assert args.plotRainfall == True
    assert args.probability == False


def test_boolean_conversion():
    assert str2bool("yes") == True
    assert str2bool("true") == True
    assert str2bool("t") == True
    assert str2bool("y") == True
    assert str2bool("1") == True
    assert str2bool("no") == False
    assert str2bool("false") == False
    assert str2bool("f") == False
    assert str2bool("n") == False
    assert str2bool("0") == False
    with pytest.raises(argparse.ArgumentTypeError):
        str2bool("maybe")

def test_str2list():
    assert str2list("arg1,arg2,arg3") == ["arg1", "arg2","arg3"]
    assert str2list("arg_unique") == ["arg_unique"]
    assert str2list("wrong.sepparator") == ["wrong.sepparator"]

if __name__ == "__main__":
    pytest.main()
