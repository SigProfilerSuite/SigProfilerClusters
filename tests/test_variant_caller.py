import os
import random

import numpy as np
import pytest
from SigProfilerSimulator import SigProfilerSimulator as sigSim
from SigProfilerClusters import SigProfilerClusters as hp


# With this seed the three example inputs get mutations repeated on three chromosomes
SIMULATION_SEED = 11


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

        # SigProfilerSimulator distributes the mutations over the chromosomes with the global
        # random and numpy generators, not with the seed file. Unseeded, some draws put every
        # mutation on a different chromosome: no simulation has an adjacent pair, there are no
        # IMDs and the analysis fails (audit register, finding 114). Seeding makes it reproducible.
        random.seed(SIMULATION_SEED)
        np.random.seed(SIMULATION_SEED)

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
            contexts="96",
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


def rerun_subclassify(caller, base_path, genome):
    """Audit finding 9: the subclassification can be repeated on the saved results, without
    the hotspot analysis (it used to fail with an undefined variable). The table that it
    writes must be the same as in the first run."""
    input_path = os.path.join(base_path, "input_vcf_sim", caller) + os.sep
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
    rerun_subclassify(caller, base_path, genome)
    return True


@pytest.mark.parametrize("caller", ["caveman", "mutect2", "standard"])
def test_caller_variants(caller):
    assert run_one_genome_test(caller) is True
