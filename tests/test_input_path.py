import os
import time

from SigProfilerClusters.classifyFunctions import get_input_vcf_path


def _touch(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("#CHROM\tPOS\n")


def test_legacy_layout_files_directly_in_input(tmp_path):
    _touch(str(tmp_path / "input" / "sample.vcf"))
    assert get_input_vcf_path(str(tmp_path) + "/") == str(tmp_path / "input") + "/"


def test_new_layout_files_in_per_run_subfolder(tmp_path):
    _touch(str(tmp_path / "input" / "proj_1234" / "sample.vcf"))
    assert (
        get_input_vcf_path(str(tmp_path) + "/")
        == str(tmp_path / "input" / "proj_1234") + "/"
    )


def test_new_layout_picks_most_recent_subfolder(tmp_path):
    old = tmp_path / "input" / "proj_old" / "a.vcf"
    new = tmp_path / "input" / "proj_new" / "a.vcf"
    _touch(str(old))
    time.sleep(0.05)
    _touch(str(new))
    os.utime(str(old.parent), (1, 1))
    assert (
        get_input_vcf_path(str(tmp_path) + "/")
        == str(tmp_path / "input" / "proj_new") + "/"
    )


def test_ignores_ds_store_and_empty_subfolders(tmp_path):
    _touch(str(tmp_path / "input" / ".DS_Store"))
    os.makedirs(str(tmp_path / "input" / "empty"))
    _touch(str(tmp_path / "input" / "proj_1" / "a.vcf"))
    assert (
        get_input_vcf_path(str(tmp_path) + "/")
        == str(tmp_path / "input" / "proj_1") + "/"
    )
