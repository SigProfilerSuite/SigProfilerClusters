import os
import re
import shutil

import pytest
from SigProfilerSimulator import SigProfilerSimulator as sigSim
from SigProfilerClusters import SigProfilerClusters as hp


def chrom_sort_key(line):
    """
    Generate a sort key that handles chromosomes in a natural order.
    Example: chr1 < chr2 < ... < chr10 < chrX < chrY
    """
    parts = line.strip().split("\t")
    chrom = parts[0]

    # Remove 'chr' if present
    chrom = chrom.replace("chr", "")

    # Map special chromosomes to a numeric order
    if chrom == "X":
        return (23, parts)
    elif chrom == "Y":
        return (24, parts)
    elif chrom == "M" or chrom == "MT":
        return (25, parts)
    else:
        try:
            return (int(chrom), parts)
        except ValueError:
            return (99, parts)  # For unknown chroms, put them at the end


def compare_files_strict(file1, file2):
    """
    Compare two files line-by-line and report mismatches.
    """
    with open(file1) as f1, open(file2) as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

    if len(lines1) != len(lines2):
        print(f"Line count mismatch: {len(lines1)} vs {len(lines2)}")
        return False

    all_match = True
    for i, (line1, line2) in enumerate(zip(lines1, lines2), start=1):
        if line1.strip() != line2.strip():
            print(
                f"Line {i} mismatch:\nExpected: {line1.strip()}\nObserved: {line2.strip()}"
            )
            all_match = False

    if all_match:
        print("Files match exactly line-by-line.")
    return all_match


def run_analysis_for_variant_callers_separate_steps(
    variant_callers, base_path, genome, simulations
):
    # Step 1: Run simulation + clustering for all callers
    for caller in variant_callers:
        print(f"Processing variant caller: {caller}")

        input_path = os.path.join(base_path, "input_vcf_sim", caller)
        if not input_path.endswith(os.sep):
            input_path += os.sep
        assert os.path.exists(input_path), f"Input path does not exist for {caller}: {input_path}"
        seed_file_path = os.path.join(base_path, "seeds_file", "Simulator_seeds.txt")

        # From SigProfilerSimulator 1.3.2 the seed file alone fixes the simulations (the
        # distribution of the mutations over the chromosomes no longer uses global generators).
        # Before, this test had to seed random and numpy as well (audit finding 114).

        # Run SigProfilerSimulator
        print(f"Running SigProfilerSimulator for {caller}...")
        sigSim.SigProfilerSimulator(
            project=caller,
            project_path=input_path,
            genome=genome,
            contexts=["96"],
            seed_file=seed_file_path,
            simulations=simulations,
        )
        print(f"SigProfilerSimulator completed for {caller}.")

        # Run SigProfilerClusters
        print(f"Running SigProfilerClusters analysis for {caller}...")
        hp.analysis(
            project=caller,
            genome=genome,
            contexts="SBS96",
            simContext=["96"],
            input_path=input_path,
            analysis="all",
            sortSims=True,
            subClassify=True,
            correction=True,
            calculateIMD=True,
            variant_caller=caller,
        )
        print(f"SigProfilerClusters analysis completed for {caller}.")

    # Step 2: After all simulations and clustering are done, compare the outputs
    for caller in variant_callers:
        print(f"Comparing output for {caller}...")
        output_file = os.path.join(
            base_path,
            "input_vcf_sim",
            caller,
            "output",
            "vcf_files_corrected",
            f"{caller}_clustered",
            f"{caller}_clustered_vaf.txt",
        )
        reference_file = os.path.join(
            base_path, "output_ref_sim", caller, f"{caller}_clustered_vaf.txt"
        )
        assert os.path.exists(output_file), f"No output was generated for {caller}: {output_file}"
        assert compare_files_strict(
            output_file, reference_file
        ), f"Files for {caller} differ: {output_file} vs {reference_file}"
        print(f"Output matches reference for {caller}.")


# The reload of the saved results does not depend on the variant caller (each caller's VAF
# reading is covered by the main run), so one of them is enough.
RERUN_SUBCLASSIFY_CALLER = "standard"


def rerun_subclassify(caller, base_path, genome):
    """Audit finding 9: the subclassification can be repeated on the saved results, without
    the hotspot analysis (it used to fail with an undefined variable). The table that it
    writes must be the same as in the first run. The path has no trailing slash. The context
    is given with its alias "96" (the first run used "SBS96"): both mean the same."""
    # Without the trailing slash (audit finding 10)
    input_path = os.path.join(base_path, "input_vcf_sim", caller)
    # The distance files of the first run are named with the standard context (audit finding 119)
    assert os.path.isdir(
        os.path.join(
            input_path, "output", "simulations", f"{caller}_intradistance_original_{genome}_SBS96"
        )
    )
    clustered_vaf = os.path.join(
        input_path, "output", "vcf_files_corrected", f"{caller}_clustered",
        f"{caller}_clustered_vaf.txt",
    )
    with open(clustered_vaf) as f:
        first_run = f.read()
    hp.analysis(
        project=caller,
        genome=genome,
        contexts="96",
        simContext=["96"],
        input_path=input_path,
        analysis="subClassify",
        subClassify=True,
        correction=True,
        variant_caller=caller,
    )
    with open(clustered_vaf) as f:
        assert f.read() == first_run, f"subClassify gave another table for {caller}"


def run_one_genome_test(caller):
    # The data live next to this file, whatever the checkout directory is called
    base_path = os.path.dirname(os.path.abspath(__file__))

    genome = "GRCh37"
    simulations = 100

    # Any exception, including an AssertionError for an output that differs from its
    # reference, propagates so that pytest reports the failure.
    run_analysis_for_variant_callers_separate_steps(
        [caller], base_path, genome, simulations
    )
    if caller == RERUN_SUBCLASSIFY_CALLER:
        rerun_subclassify(caller, base_path, genome)
    return True


@pytest.mark.parametrize("caller", ["caveman", "mutect2", "standard"])
def test_caller_variants(caller):
    assert run_one_genome_test(caller) is True


def test_simulations_are_reproducible_with_the_seed_file(tmp_path):
    """Two simulations with the same seed file and the same input give the same files, without
    seeding random or numpy. SigProfilerSimulator 1.3.2 is the first version that does it
    (audit findings 4, 114 and 116); the integration tests above rely on it."""
    base_path = os.path.dirname(os.path.abspath(__file__))
    seed_file = os.path.join(base_path, "seeds_file", "Simulator_seeds.txt")
    vcf = os.path.join(base_path, "input_vcf_sim", "standard", "standard.vcf")

    runs = []
    for number in range(2):
        project_path = tmp_path / f"run{number}"
        project_path.mkdir()
        shutil.copy(vcf, project_path)
        sigSim.SigProfilerSimulator(
            project="standard",
            project_path=str(project_path) + os.sep,
            genome="GRCh37",
            contexts=["96"],
            seed_file=seed_file,
            simulations=5,
        )
        runs.append(project_path / "output" / "simulations")

    files = sorted(p.relative_to(runs[0]) for p in runs[0].rglob("*") if p.is_file())
    assert len(files) == 5, files
    assert files == sorted(p.relative_to(runs[1]) for p in runs[1].rglob("*") if p.is_file())
    different = [str(f) for f in files if (runs[0] / f).read_bytes() != (runs[1] / f).read_bytes()]
    assert not different, f"simulations differ between two runs with the same seed file: {different}"


UUID = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")


def _data_files(output):
    """{path relative to output: file} of the data files of a run. The logs (time stamps, paths)
    and the PDFs (creation date) change with every run and are not results. Newer versions of
    SigProfilerMatrixGenerator copy their input into folders named <name>_<uuid>, with a new
    uuid in every run: it is replaced in the key so that the files of two runs can be paired."""
    found = {}
    for path in output.rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(output)
        if "logs" in relative.parts or path.suffix == ".pdf":
            continue
        key = UUID.sub("UUID", relative.as_posix())
        assert key not in found, f"two files with the same name once the uuid is removed: {key}"
        found[key] = path
    return found


def test_data_files_pairs_the_folders_that_only_differ_in_their_uuid(tmp_path):
    for run, uuid in (("a", "55b841a4-cf24-4571-a1ba-7fe55cfb274a"), ("b", "4b44f91d-cfda-4b3e-9231-a4ed02bd1a80")):
        folder = tmp_path / run / "vcf_files" / "SNV" / "input" / f"standard_clustered_{uuid}"
        folder.mkdir(parents=True)
        (folder / "standard_clustered.txt").write_text("same")
        (tmp_path / run / "logs").mkdir()
        (tmp_path / run / "logs" / f"{run}.out").write_text("different")
        (tmp_path / run / "plot.pdf").write_bytes(b"date")
    first, second = _data_files(tmp_path / "a"), _data_files(tmp_path / "b")
    assert sorted(first) == sorted(second)
    assert list(first) == ["vcf_files/SNV/input/standard_clustered_UUID/standard_clustered.txt"]


def test_full_pipeline_is_repeatable_with_the_seed_file(tmp_path):
    """Audit finding 4: two complete runs (SigProfilerSimulator and SigProfilerClusters) on the
    same input and the same seed file give identical results. It needs SigProfilerSimulator
    1.3.2 (before, the simulations differed between runs)."""
    base_path = os.path.dirname(os.path.abspath(__file__))
    seed_file = os.path.join(base_path, "seeds_file", "Simulator_seeds.txt")
    vcf = os.path.join(base_path, "input_vcf_sim", "standard", "standard.vcf")

    outputs = []
    for number in range(2):
        project_path = tmp_path / f"run{number}"
        project_path.mkdir()
        shutil.copy(vcf, project_path)
        path = str(project_path) + os.sep
        sigSim.SigProfilerSimulator(
            project="standard", project_path=path, genome="GRCh37", contexts=["96"],
            seed_file=seed_file, simulations=100,
        )
        hp.analysis(
            project="standard", genome="GRCh37", contexts="SBS96", simContext=["96"],
            input_path=path, analysis="all", sortSims=True, subClassify=True, correction=True,
            calculateIMD=True, variant_caller="standard", max_cpu=2,
        )
        outputs.append(project_path / "output")

    files, other = _data_files(outputs[0]), _data_files(outputs[1])
    assert len(files) > 100, "too few result files to be a meaningful comparison"
    assert sorted(files) == sorted(other)
    different = [key for key in files if files[key].read_bytes() != other[key].read_bytes()]
    assert not different, f"results differ between two runs with the same seed file: {different[:10]}"
