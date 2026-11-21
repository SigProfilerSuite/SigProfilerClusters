"""Shared helpers for the tests of classifyFunctions.findClustersOfClusters(_noVAF).

They build the same small synthetic clustered SBS table (with or without a VAF column)
inside tmp_path and read back the subclass files written by the functions.
No genome reference is required.
"""

import pickle

PROJECT = "proj"
GENOME = "GRCh37"
IMD_CUTOFF = 100  # sample-wide IMD: a distance above it is not clustered
HEADER = (
    "project\tsamples\tID\tgenome\tmutType\tchr\tstart\tend\tref\talt\t"
    "mutClass\tIMDplot\tIMD\tprobability\n"
)
HEADER_VAF = (
    "project\tsamples\tID\tgenome\tmutType\tchr\tstart\tend\tref\talt\t"
    "mutClass\tIMDplot\tIMD\tVAF/CCF\n"
)
SUBCLASSES = [
    "class1",
    "class2",
    "class3",
    "class1a",
    "class1b",
    "class1c",
    "class2Y",
    "class2K",
    "class2S",
    "class2N",
]


def no_matrices(*args, **kwargs):
    """Stand-in for generateMatrices (needs a genome); must be importable by the pool."""


def row(pos, ref="C", alt="T", sample="S1", chrom="1", vaf=None):
    """One clustered SBS: project, sample, ID, genome, type, chr, start, end, ref, alt,
    class, IMDplot, IMD and a last column: the VAF when given, else a probability of 1."""
    cols = [PROJECT, sample, ".", GENOME, "SNP", chrom, str(pos), str(pos), ref, alt]
    cols += ["SOMATIC", "1", "1", "1" if vaf is None else str(vaf)]
    return "\t".join(cols) + "\n"


def setup(tmp_path, rows, header, table_name, imds=None, correction=False, chrom_based=False):
    """Write the clustered table, named `table_name` inside the project folder, and the
    IMD pickles read by the functions. Returns the clustered project folder."""
    imds = {"S1": IMD_CUTOFF} if imds is None else imds
    vcf_dir = "vcf_files_corrected" if correction else "vcf_files"
    clustered = tmp_path / "output" / vcf_dir / f"{PROJECT}_clustered"
    (clustered / "SNV").mkdir(parents=True)
    table = clustered / table_name
    table.parent.mkdir(parents=True, exist_ok=True)
    table.write_text(header + "".join(rows))

    data = tmp_path / "output" / "simulations" / "data"
    data.mkdir(parents=True)
    suffix = "_chrom" if chrom_based else ""
    with open(str(data / f"imds{suffix}.pickle"), "wb") as f:
        pickle.dump(imds, f)
    if correction:
        with open(str(data / "imds_corrected.pickle"), "wb") as f:
            pickle.dump({sample: {} for sample in imds}, f)
    return clustered


def subclass_path(clustered, name, chrom_based=False):
    suffix = "_chrom" if chrom_based else ""
    return clustered / f"subclasses{suffix}" / name / f"{PROJECT}_clustered_{name}.txt"


def entries(clustered, name, chrom_based=False):
    """Data rows (header skipped) of one subclass file, as lists of columns."""
    lines = subclass_path(clustered, name, chrom_based).read_text().splitlines()
    return [line.split("\t") for line in lines[1:]]


def positions(clustered, name, chrom_based=False):
    return [int(cols[6]) for cols in entries(clustered, name, chrom_based)]


def populated(clustered):
    return sorted(name for name in SUBCLASSES if entries(clustered, name))
