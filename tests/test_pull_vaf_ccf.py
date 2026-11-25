"""Unit tests for classifyFunctions.pullCCF and pullVaf: they read the VAF/CCF of every
mutation from the original input VCFs and append it to the clustered SBS table.
Inputs are small synthetic VCFs and tables built with tmp_path. No genome required.
Run with: pytest tests/test_pull_vaf_ccf.py -v
"""

import pytest

from clustered_table_helpers import HEADER, PROJECT, row, setup
from SigProfilerClusters.classifyFunctions import pullCCF, pullVaf

NO_VAF = "nan"  # value written when the input file carries no usable VAF/CCF (audit finding 6)


def _vcf(*records):
    """VCF text from tab-separated records; the ## header lines are comments."""
    return "##fileformat=VCFv4.2\n" + "".join("\t".join(r) + "\n" for r in records)


def _table(tmp_path, rows, correction=True):
    return setup(tmp_path, rows, HEADER, f"SNV/{PROJECT}_clustered.txt", None, correction)


def _write_vcf(tmp_path, name, text, subdir=None):
    folder = tmp_path / "input" if subdir is None else tmp_path / "input" / subdir
    folder.mkdir(parents=True, exist_ok=True)
    (folder / name).write_text(text)


def _vaf_table(clustered):
    lines = (clustered / f"{PROJECT}_clustered_vaf.txt").read_text().splitlines()
    return [line.split("\t") for line in lines]


def _values(clustered):
    """The appended VAF/CCF column of every data row."""
    return [cols[-1] for cols in _vaf_table(clustered)[1:]]


def _snv(chrom, pos, ref, alt, info=".", fmt=None, tumor=None, last=None):
    """A VCF record with the 8 fixed columns plus optional FORMAT/sample columns."""
    rec = [chrom, str(pos), ".", ref, alt, ".", "PASS", info]
    if fmt is not None:
        rec += [fmt, tumor]
    if last is not None:
        rec += [last]
    return rec


# ─── pullCCF ─────────────────────────────────────────────────


class TestPullCCF:
    def test_ccf_is_appended_to_each_clustered_mutation(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T"), row(200, "G", "A")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", last="0.4"), _snv("1", 200, "G", "A", last="0.7")))
        pullCCF(PROJECT, str(tmp_path) + "/")
        assert _values(clustered) == ["0.4", "0.7"]

    def test_output_has_the_vaf_header(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", last="0.4")))
        pullCCF(PROJECT, str(tmp_path) + "/")
        header = _vaf_table(clustered)[0]
        assert header[0] == "project" and header[-1] == "VAF/CCF"

    def test_original_row_columns_are_kept(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", last="0.4")))
        pullCCF(PROJECT, str(tmp_path) + "/")
        cols = _vaf_table(clustered)[1]
        assert cols[1] == "S1" and cols[5:10] == ["1", "100", "100", "C", "T"]

    def test_chr_prefix_in_the_vcf_is_ignored(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("chr1", 100, "C", "T", last="0.4")))
        pullCCF(PROJECT, str(tmp_path) + "/")
        assert _values(clustered) == ["0.4"]

    def test_record_without_a_numeric_ccf_gets_the_placeholder(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", info="DP=30")))
        pullCCF(PROJECT, str(tmp_path) + "/")
        assert _values(clustered) == [NO_VAF]

    def test_doublet_ccf_is_assigned_to_both_bases(self, tmp_path):
        rows = [row(100, "A", "G"), row(101, "C", "T")]
        clustered = _table(tmp_path, rows)
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "AC", "GT", last="0.55")))
        pullCCF(PROJECT, str(tmp_path) + "/")
        assert _values(clustered) == ["0.55", "0.55"]

    def test_each_sample_uses_its_own_file(self, tmp_path):
        rows = [row(100, "C", "T", sample="S1"), row(100, "C", "T", sample="S2")]
        clustered = _table(tmp_path, rows)
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", last="0.2")))
        _write_vcf(tmp_path, "S2.vcf", _vcf(_snv("1", 100, "C", "T", last="0.9")))
        pullCCF(PROJECT, str(tmp_path) + "/")
        assert _values(clustered) == ["0.2", "0.9"]

    def test_sample_name_is_the_file_name_before_the_first_dot(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T", sample="S1")])
        _write_vcf(tmp_path, "S1.somatic.vcf", _vcf(_snv("1", 100, "C", "T", last="0.4")))
        pullCCF(PROJECT, str(tmp_path) + "/")
        assert _values(clustered) == ["0.4"]

    def test_uncorrected_table_is_read_from_vcf_files(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")], correction=False)
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", last="0.4")))
        pullCCF(PROJECT, str(tmp_path) + "/", correction=False)
        assert "vcf_files/" in str(clustered) and _values(clustered) == ["0.4"]

    def test_inputs_in_a_per_run_subfolder_are_found(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", last="0.4")), subdir="proj_1234")
        pullCCF(PROJECT, str(tmp_path) + "/")
        assert _values(clustered) == ["0.4"]

    def test_mutation_missing_from_the_input_raises(self, tmp_path):
        _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 999, "C", "T", last="0.4")))
        with pytest.raises(KeyError):
            pullCCF(PROJECT, str(tmp_path) + "/")


# ─── pullVaf: standard ───────────────────────────────────────


class TestPullVafStandard:
    def test_vaf_in_the_info_column(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", info="DP=30;VAF=0.35;X=1")))
        pullVaf(PROJECT, str(tmp_path) + "/", "standard")
        assert _values(clustered) == ["0.35"]

    def test_af_in_the_sample_column(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", fmt="GT", tumor="AF=0.6;DP=20")))
        pullVaf(PROJECT, str(tmp_path) + "/", "standard")
        assert _values(clustered) == ["0.6"]

    def test_info_vaf_takes_precedence_over_the_sample_column(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        rec = _snv("1", 100, "C", "T", info="VAF=0.1", fmt="GT", tumor="AF=0.9")
        _write_vcf(tmp_path, "S1.vcf", _vcf(rec))
        pullVaf(PROJECT, str(tmp_path) + "/", "standard")
        assert _values(clustered) == ["0.1"]

    def test_record_without_vaf_gets_the_placeholder(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", fmt="GT", tumor="0/1")))
        pullVaf(PROJECT, str(tmp_path) + "/", "standard")
        assert _values(clustered) == [NO_VAF]

    def test_unparsable_vaf_gets_the_placeholder(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", info="VAF=high")))
        pullVaf(PROJECT, str(tmp_path) + "/", "standard")
        assert _values(clustered) == [NO_VAF]

    def test_af_in_the_info_column_is_read(self, tmp_path):
        # audit finding 6: it used to be ignored and recorded as missing
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", info="DP=20;AF=0.31")))
        pullVaf(PROJECT, str(tmp_path) + "/", "standard")
        assert _values(clustered) == ["0.31"]

    def test_tvaf_is_not_mistaken_for_vaf(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", info="TVAF=0.9;VAF=0.31")))
        pullVaf(PROJECT, str(tmp_path) + "/", "standard")
        assert _values(clustered) == ["0.31"]

    def test_only_tvaf_present_means_no_vaf(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", info="TVAF=0.9")))
        pullVaf(PROJECT, str(tmp_path) + "/", "standard")
        assert _values(clustered) == [NO_VAF]

    def test_missing_vafs_are_reported_per_sample(self, tmp_path):
        _table(tmp_path, [row(100, "C", "T"), row(200, "G", "A")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", info="VAF=0.3"), _snv("1", 200, "G", "A", info="DP=9")))
        with pytest.warns(UserWarning, match="1 of 2 clustered mutations of sample S1 have no VAF"):
            pullVaf(PROJECT, str(tmp_path) + "/", "standard")

    def test_no_warning_when_every_vaf_is_known(self, tmp_path):
        import warnings

        _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", info="VAF=0.3")))
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            pullVaf(PROJECT, str(tmp_path) + "/", "standard")

    def test_chr_prefixed_vcf_still_matches_the_table(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("chr1", 100, "C", "T", info="VAF=0.35")))
        pullVaf(PROJECT, str(tmp_path) + "/", "standard")
        assert _values(clustered) == ["0.35"]

    def test_doublet_vaf_is_assigned_to_both_bases(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "A", "G"), row(101, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "AC", "GT", info="VAF=0.42")))
        pullVaf(PROJECT, str(tmp_path) + "/", "standard")
        assert _values(clustered) == ["0.42", "0.42"]

    def test_mutation_missing_from_the_input_raises(self, tmp_path):
        _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 999, "C", "T", info="VAF=0.3")))
        with pytest.raises(KeyError):
            pullVaf(PROJECT, str(tmp_path) + "/", "standard")


# ─── pullVaf: caveman ────────────────────────────────────────


class TestPullVafCaveman:
    def test_vaf_is_the_last_colon_separated_value(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        rec = _snv("1", 100, "C", "T", last="0:12:3:0.25")
        _write_vcf(tmp_path, "S1.vcf", _vcf(rec))
        pullVaf(PROJECT, str(tmp_path) + "/", "caveman")
        assert _values(clustered) == ["0.25"]

    def test_chr_prefix_is_removed(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("chr1", 100, "C", "T", last="a:b:0.5")))
        pullVaf(PROJECT, str(tmp_path) + "/", "caveman")
        assert _values(clustered) == ["0.5"]

    def test_a_line_without_vaf_does_not_stop_the_file(self, tmp_path):
        # audit finding 6: the first line without VAF used to stop the whole file
        clustered = _table(tmp_path, [row(100, "C", "T"), row(200, "G", "A")])
        records = _vcf(_snv("1", 100, "C", "T", last="not_a_number"), _snv("1", 200, "G", "A", last="a:b:0.5"))
        _write_vcf(tmp_path, "S1.vcf", records)
        pullVaf(PROJECT, str(tmp_path) + "/", "caveman")
        assert _values(clustered) == [NO_VAF, "0.5"]


# ─── pullVaf: mutect2 and custom ─────────────────────────────

MUTECT_HEADER = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR"]


def _mutect_vcf(*records):
    return "##fileformat=VCFv4.2\n" + "\t".join(MUTECT_HEADER) + "\n" + "".join("\t".join(r) + "\n" for r in records)


def _mutect_rec(chrom, pos, ref, alt, fmt, tumor):
    return [chrom, str(pos), ".", ref, alt, ".", "PASS", ".", fmt, "0/0:0", tumor]


class TestPullVafMutect2:
    def test_af_is_read_from_the_tumor_column(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _mutect_vcf(_mutect_rec("1", 100, "C", "T", "GT:AF", "0/1:0.31")))
        pullVaf(PROJECT, str(tmp_path) + "/", "mutect2")
        assert _values(clustered) == ["0.31"]

    def test_chr_prefix_is_removed_case_insensitively(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _mutect_vcf(_mutect_rec("CHR1", 100, "C", "T", "GT:AF", "0/1:0.31")))
        pullVaf(PROJECT, str(tmp_path) + "/", "mutect2")
        assert _values(clustered) == ["0.31"]

    def test_af_position_follows_the_format_column(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _mutect_vcf(_mutect_rec("1", 100, "C", "T", "GT:DP:AF", "0/1:40:0.12")))
        pullVaf(PROJECT, str(tmp_path) + "/", "mutect2")
        assert _values(clustered) == ["0.12"]

    def test_doublet_af_is_assigned_to_both_bases(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "A", "G"), row(101, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _mutect_vcf(_mutect_rec("1", 100, "AC", "GT", "GT:AF", "0/1:0.2")))
        pullVaf(PROJECT, str(tmp_path) + "/", "mutect2")
        assert _values(clustered) == ["0.2", "0.2"]

    def test_record_without_af_in_format_has_an_unknown_vaf(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _mutect_vcf(_mutect_rec("1", 100, "C", "T", "GT:DP", "0/1:40")))
        pullVaf(PROJECT, str(tmp_path) + "/", "mutect2")
        assert _values(clustered) == [NO_VAF]

    def test_file_without_a_tumor_column_raises_a_clear_error(self, tmp_path):
        _table(tmp_path, [row(100, "C", "T")])
        text = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n1\t100\t.\tC\tT\t.\tPASS\t.\tGT:AF\t0/1:0.3\n"
        _write_vcf(tmp_path, "S1.vcf", text)
        with pytest.raises(ValueError, match="tumor sample column of S1.vcf was not found"):
            pullVaf(PROJECT, str(tmp_path) + "/", "mutect2")

    def test_tumor_sample_named_in_the_header_is_used(self, tmp_path):
        # audit finding 6: a tumor sample that is not called TUMOR
        clustered = _table(tmp_path, [row(100, "C", "T")])
        text = (
            "##tumor_sample=SAMPLE_01\n"
            + "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL_01", "SAMPLE_01"])
            + "\n1\t100\t.\tC\tT\t.\tPASS\t.\tGT:AF\t0/0:0.9\t0/1:0.3\n"
        )
        _write_vcf(tmp_path, "S1.vcf", text)
        pullVaf(PROJECT, str(tmp_path) + "/", "mutect2")
        assert _values(clustered) == ["0.3"]


class TestPullVafCustom:
    def _vcf(self, fmt, tumor):
        return _mutect_vcf(_mutect_rec("1", 100, "C", "T", fmt, tumor))

    def test_default_column_and_field(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", self._vcf("GT:AF", "0/1:0.44"))
        pullVaf(PROJECT, str(tmp_path) + "/", "custom")
        assert _values(clustered) == ["0.44"]

    def test_custom_field_name(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", self._vcf("GT:VF", "0/1:0.27"))
        pullVaf(PROJECT, str(tmp_path) + "/", "custom", vaf_field="VF")
        assert _values(clustered) == ["0.27"]

    def test_custom_column_index(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        # column 9 is the NORMAL sample in the synthetic VCF, column 10 the TUMOR
        _write_vcf(tmp_path, "S1.vcf", self._vcf("GT:AF", "0/1:0.44"))
        pullVaf(PROJECT, str(tmp_path) + "/", "custom", tumor_vaf_column=10)
        assert _values(clustered) == ["0.44"]

    def test_column_beyond_the_header_raises_a_clear_error(self, tmp_path):
        _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", self._vcf("GT:AF", "0/1:0.44"))
        with pytest.raises(ValueError, match="Tumor column 20 not found in the header of S1.vcf"):
            pullVaf(PROJECT, str(tmp_path) + "/", "custom", tumor_vaf_column=20)


# ─── common behavior ─────────────────────────────────────────


class TestPullVafCommon:
    def test_output_has_header_and_original_columns(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", info="VAF=0.35")))
        pullVaf(PROJECT, str(tmp_path) + "/", "standard")
        table = _vaf_table(clustered)
        assert table[0][-1] == "VAF/CCF"
        assert table[1][1] == "S1" and table[1][-1] == "0.35"

    def test_inputs_in_a_per_run_subfolder_are_found(self, tmp_path):
        clustered = _table(tmp_path, [row(100, "C", "T")])
        text = _vcf(_snv("1", 100, "C", "T", info="VAF=0.35"))
        _write_vcf(tmp_path, "S1.vcf", text, subdir="proj_1234")
        pullVaf(PROJECT, str(tmp_path) + "/", "standard")
        assert _values(clustered) == ["0.35"]

    def test_unknown_variant_caller_fails_with_a_clear_error(self, tmp_path):
        _table(tmp_path, [row(100, "C", "T")])
        _write_vcf(tmp_path, "S1.vcf", _vcf(_snv("1", 100, "C", "T", info="VAF=0.35")))
        with pytest.raises(ValueError, match="variant_caller must be"):
            pullVaf(PROJECT, str(tmp_path) + "/", "unknown_caller")
