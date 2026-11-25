#!/usr/bin/env python3
from . import constants
from .exceptions import (
    InvalidParameterError,
    MissingInputError,
    MissingSimulationsError,
)
import os
import re
import argparse
from SigProfilerMatrixGenerator.scripts import (
    convert_input_to_simple_files as convertIn,
)
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from . import hotspot
from . import classifyFunctions
import SigProfilerMatrixGenerator as sig
import matplotlib as plt
import sigProfilerPlotting as sigPlt
import datetime
import shutil
import platform
import numpy as np
import statistics
import scipy
import time
import functools
import numbers
import sys

# from SigProfilerExtractor import sigpro as sigs
import multiprocessing as mp
from . import plottingFunctions
from SigProfilerMatrixGenerator.scripts import MutationMatrixGenerator as matRef
from . import convertToVCF
import random
import pickle
import warnings
import pandas as pd
import bisect
import seaborn as sns


# Values of output_type that analysis() accepts. None is what analysis="hotspot" and
# analysis="subClassify" use internally to skip the calculation of the distances.
SUPPORTED_OUTPUT_TYPES = ("all", None)

def log_stamp():
    """Date and time for the names of the log files, so that two runs of the same project on
    the same day do not overwrite each other's logs (audit finding 31)"""
    return datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")


_CONTEXT_NAMES = (
    "SBS96 (also accepted: '96', 'SBS') or ID83 (also accepted: 'ID', '83', 'INDEL')"
)
_DBS_NOTE = " (DBS is one of the outputs of SigProfilerClusters, not an input context)"


def validate_arguments(
    analysis, contexts, simContext, windowSize, variant_caller, interdistance="96"
):
    """Checks the arguments of analysis() that have a closed set of valid values, before any
    file is read, created or deleted. Raises InvalidParameterError naming the value and the valid
    ones. (genome, output_type and chrom_based have their own checks.)"""
    if analysis not in constants.SUPPORTED_ANALYSES:
        raise InvalidParameterError(
            "analysis must be one of %s, not %r"
            % (", ".join(map(repr, constants.SUPPORTED_ANALYSES)), analysis)
        )
    if constants.canonical_context(contexts) is None:
        raise InvalidParameterError(
            "contexts must be %s, not %r%s"
            % (_CONTEXT_NAMES, contexts, _DBS_NOTE if contexts == "DBS" else "")
        )
    if interdistance != "SI" and constants.canonical_context(interdistance) is None:
        raise InvalidParameterError(
            "interdistance must be %s or 'SI', not %r" % (_CONTEXT_NAMES, interdistance)
        )
    if not isinstance(simContext, (list, tuple)) or not all(
        isinstance(name, str) and name for name in simContext
    ):
        raise InvalidParameterError(
            "simContext must be a list with the name of the simulated context, "
            "for example ['96'], not %r" % (simContext,)
        )
    if len(simContext) != 1:
        raise InvalidParameterError(
            "Please only provide one type of context for analysis (simContext has "
            + str(len(simContext))
            + ")."
        )
    if (
        not isinstance(windowSize, numbers.Integral)
        or isinstance(windowSize, bool)
        or windowSize <= 0
    ):
        raise InvalidParameterError(
            "windowSize must be a positive integer, not %r" % (windowSize,)
        )
    if variant_caller not in constants.VARIANT_CALLERS:
        raise InvalidParameterError(
            "variant_caller must be one of %s, not %r"
            % (", ".join(map(repr, constants.VARIANT_CALLERS)), variant_caller)
        )


def check_input_folder(input_path):
    """The project folder must exist and have an input folder with the mutation files (the
    conversion reads them in every kind of analysis). Raises MissingInputError."""
    if not isinstance(input_path, str) or not os.path.isdir(input_path):
        raise MissingInputError(
            "The project folder input_path does not exist: %r" % (input_path,)
        )
    input_folder = os.path.join(input_path, "input")
    if not os.path.isdir(input_folder):
        raise MissingInputError(
            "The project folder %r has no 'input' folder with the mutation files "
            "(expected %r)" % (input_path, input_folder)
        )
    files_folder = classifyFunctions.get_input_vcf_path(input_path)
    if not any(
        x != ".DS_Store" and os.path.isfile(os.path.join(files_folder, x))
        for x in os.listdir(files_folder)
    ):
        raise MissingInputError(
            "The input folder %r has no input files" % (input_folder,)
        )


# Seed used to choose the simulations of eventProbability (audit finding 11)
EVENT_PROBABILITY_SEED = 0


def distance_multiple_files_sims(
    output_path,
    simulations,
    simulation_path,
    simulation_path_sorted,
    file_context2,
    genome,
    centromeres,
    sortSims,
):
    """
    Calculates the IMDs across all mutations in the simulated samples.

    Parameters:
                               output_path	->	path where the IMD distance files will be written (string)
                               simulations	->	list of simulations for a given processor (list)
                       simulation_path	->	path to simulation files (string)
            simulation_path_sorted	->	path to sorted simulation files (string)
                             file_context2	->	combined context file suffix for file naming (string)
                                            genome	->	the reference genome used for the given analysis (string)
                               centromeres	->	a dictionary of the locations of all centromeres for a given reference genome (dictionary)
                                      sortSims	->  optional parameter that requires the function to sort the simulations or not (boolean; default=True)

    Returns:
            None

    Outputs:
            Distance files for all mutations in each sample across each simulation.
    """

    # Start iterating through the simulated files (parallelized)
    for file in simulations:

        # Reads in the all of the mutations and sorts them by sample, chromosome, and position if sortSims=True
        with open(simulation_path + file) as f:
            lines = [line.strip().split() for line in f]
        if sortSims:
            original_lines = sorted(
                lines[1:],
                key=lambda x: (
                    x[15],
                    [
                        "X",
                        "Y",
                        "1",
                        "2",
                        "3",
                        "4",
                        "5",
                        "6",
                        "7",
                        "8",
                        "9",
                        "10",
                        "11",
                        "12",
                        "13",
                        "14",
                        "15",
                        "16",
                        "17",
                        "18",
                        "19",
                        "20",
                        "21",
                        "22",
                    ].index(x[4]),
                    int(x[5]),
                ),
            )
            with open(simulation_path_sorted + file, "w") as f2:
                for line in original_lines:
                    print("\t".join([x for x in line]), file=f2)
        else:
            original_lines = lines[
                1:
            ]  # Skips header that should be in each simulated files

        # Clustered nalysis requires >1 mutation
        if len(original_lines) > 1:
            distances = [None] * (
                len(original_lines) - 1
            )  # Instantiates an empty dataframe to place the IMD's
            i = 0

            # Iterate through each mutation and determine the minimum distance between the upstream and downstream adjacent mutations.
            # Distances across centromeres are not considered (ie the final mutation on the p-arm only has a single IMD to the upstream mutations,
            # while the first mutation on the q-arm only has a single IMD to the downstream mutation).
            for x, y in zip(original_lines, original_lines[1:]):

                # Collect the relecant information the two current mutations
                prevSamp = x[15]
                prevChrom = x[4]
                prevPos = int(x[5])
                prevRef = x[10]
                prevMut = x[12]
                currentSamp = y[15]
                currentChrom = y[4]
                currentPos = int(y[5])
                centroStart = centromeres[genome][prevChrom][0]
                centroEnd = centromeres[genome][prevChrom][1]

                # Determine the IMD for the pair of mutations
                if (
                    prevSamp == currentSamp
                ):  # Ensure that the same sample is being considered
                    if (
                        prevChrom == currentChrom
                    ):  # Ensure that the same chromosome is being considered
                        if (
                            currentPos < centroStart
                        ):  # Ensure both mutations are on the p-arm
                            distances[i] = [
                                currentPos
                                - (prevPos - 2 + len(prevRef) + len(prevMut)),
                                prevSamp,
                                prevChrom,
                                prevPos,
                                prevRef,
                                prevMut,
                                "c",
                            ]
                        elif (
                            prevPos > centroEnd
                        ):  # Or ensure both mutations are on the q-arm
                            distances[i] = [
                                currentPos
                                - (prevPos - 2 + len(prevRef) + len(prevMut)),
                                prevSamp,
                                prevChrom,
                                prevPos,
                                prevRef,
                                prevMut,
                                "c",
                            ]
                        else:  # Else current mutations are on different arms, so calculate IMD's accordingly
                            distances[i] = [
                                prevPos
                                - (
                                    int(original_lines[i - 1][5])
                                    - 2
                                    + len(original_lines[i - 1][10])
                                    + len(original_lines[i - 1][12])
                                ),
                                prevSamp,
                                prevChrom,
                                prevPos,
                                prevRef,
                                prevMut,
                                "d",
                            ]  # Place a 'd' marker to represent the start of the q-arm
                    else:
                        distances[i] = [
                            "a",
                            prevSamp,
                            prevChrom,
                            prevPos,
                            prevRef,
                            prevMut,
                        ]  # Place an 'a' marker to reflect a new chromosome
                else:
                    distances[i] = [
                        "a",
                        prevSamp,
                        prevChrom,
                        prevPos,
                        prevRef,
                        prevMut,
                    ]  # Place an 'a' marker to reflect a new sample
                i += 1

            # Save the final distances that reflect the minimum IMD for each mutations
            # ensuring that the mutations are not on different samples, chromosomes, or arms of the chromosomes
            # using the simple character code in the above code
            distances_new = distances + [distances[-1]]
            final_distances = [
                (
                    "bb"
                    if x[1] != y[1]
                    else (
                        [min(x[0], y[0])] + y[1:-1]
                        if x[0] != "a" and y[0] != "a" and x[-1] != "d"
                        else (
                            [y[0]] + y[1:-1]
                            if x[-1] == "d" and x[0] != "a" and y[0] != "a"
                            else (
                                x
                                if y[0] == "a" and x[0] != "a"
                                else y if x[0] == "a" and y[0] != "a" else "bb"
                            )
                        )
                    )
                )
                for x, y in zip(distances_new, distances_new[1:])
            ]
            final_distances = [distances_new[0]] + final_distances
            final_distances_filtered = [
                x for x in final_distances if x[0] != "b" and x[0] != "a"
            ]

            if (
                len(final_distances_filtered) == 0
            ):  # Continue if there are no IMDs (ex. only two mutation are a chromosome and occured on different arms)
                continue

            # Format the sample name to use for the final file name
            # Also create the respective output folders
            sample = final_distances_filtered[0][1]
            file_sample = sample.split("_")[:-1]
            file_sample = ("_").join([x for x in file_sample])
            # exist_ok: several processes may create the same folder
            os.makedirs(output_path + file_sample + "/", exist_ok=True)

            # Open the current sample file and start iterating through the final distances
            # If a new sample is reached, close the current file and open the new one
            # The files are written under the simulations/ folder with the following fields:
            # <IMD> <Sample_sim#> <chrom> <pos> <ref> <alt>
            out_file = open(
                output_path
                + file_sample
                + "/"
                + sample
                + "_"
                + file_context2
                + "_intradistance.txt",
                "a",
                10000000,
            )
            for x in final_distances_filtered:
                if x[1] == sample:
                    print("\t".join([str(y) for y in x]), file=out_file)
                else:
                    out_file.close()
                    sample = x[1]
                    file_sample = sample.split("_")[:-1]
                    file_sample = ("_").join([x for x in file_sample])
                    if sample.split("_")[0] != file_sample:
                        file_sample = sample.split("_")[:-1]
                        file_sample = ("_").join([x for x in file_sample])
                        os.makedirs(output_path + file_sample + "/", exist_ok=True)
                        out_file = open(
                            output_path
                            + file_sample
                            + "/"
                            + sample
                            + "_"
                            + file_context2
                            + "_intradistance.txt",
                            "a",
                            10000000,
                        )
                        print("\t".join([str(y) for y in x]), file=out_file)
                    else:
                        os.makedirs(output_path + file_sample + "/", exist_ok=True)
                        out_file = open(
                            output_path
                            + file_sample
                            + "/"
                            + sample
                            + "_"
                            + file_context2
                            + "_intradistance.txt",
                            "a",
                        )
                        print("\t".join([str(y) for y in x]), file=out_file)
            out_file.close()


def chromosome_rank(chrom):
    """Sort key of a chromosome name in the usual order: 1, 2, ..., 22, X, Y, MT, then the rest"""
    name = str(chrom)
    if name.lower().startswith("chr"):
        name = name[3:]
    if name.isdigit():
        return (0, int(name), "")
    if name in ("X", "Y", "M", "MT"):
        return (1, ("X", "Y", "M", "MT").index(name), "")
    return (2, 0, name)


def sort_distance_files(directory):
    """Sorts the rows of every <sample>_<context>_intradistance.txt of a folder by chromosome
    (natural order) and position.

    The distances of an original sample are written by several workers, each one appending the
    chromosomes that it was given, so the order of the rows depended on which worker finished
    first and changed from one run to another, and with it the order of the tables and the
    numbers of the groups built from them (audit findings 121 and 4). The rows of one
    chromosome are already in position order; the sort is stable, so alleles of the same site
    keep their order."""
    for name in sorted(os.listdir(directory)):
        if not name.endswith("_intradistance.txt"):
            continue
        path = os.path.join(directory, name)
        with open(path) as f:
            rows = [line for line in f.read().splitlines() if line]
        rows.sort(key=lambda row: (chromosome_rank(row.split("\t")[2]), int(row.split("\t")[3])))
        temporary = path + ".sorting"
        with open(temporary, "w") as f:
            f.write("".join(row + "\n" for row in rows))
        os.replace(temporary, path)


def distance_one_file(
    sample_path,
    original_samples,
    output_path_original,
    file_context2,
    genome,
    centromeres,
):
    """
    #################################################################
    ############### Needs to be parallelized ########################
    #################################################################

    Calculates the IMDs across all mutations in the original samples.

    Parameters:
                             sample_path	->	path to the newly generated mutation files for the original samples (string)
            output_path_original	->	path where the IMD distance files will be written for the original samples (string)
                       file_context2	->	combined context file suffix for file naming (string)
                                      genome	->	the reference genome used for the given analysis (string)
                             centromeres	->	a dictionary of the locations of all centromeres for a given reference genome (dictionary)

    Returns:
            None

    Outputs:
            Distance files for all mutations in each original sample.
    """

    # # Organize the output paths (ie remove if they exist, and then create and emtpy folder)
    # if os.path.exists(output_path_original):
    # 	shutil.rmtree(output_path_original)
    # 	os.makedirs(output_path_original)
    # else:
    # 	os.makedirs(output_path_original)

    # Iterate through each chromomse file that contains all of the mutations
    # from the original samples. Sample_path is found under vcf_files_[suffix]
    # for file in os.listdir(sample_path):
    for file in sample_path:

        # Read in the files and sort them
        # with open(sample_path + file) as f:
        with open(original_samples + file) as f:

            lines = [line.strip().split() for line in f]
        sorted_lines = sorted(lines, key=lambda x: (x[0], int(x[2])))

        # Audit finding 5: the file of a chromosome holds the mutations of every sample, and the
        # IMD of the first and last mutation of a sample, and of a sample with a single
        # mutation, used to be computed with a neighbour of another sample. Each sample is now
        # processed on its own, so no distance crosses a sample boundary.
        lines_by_sample = {}
        for line in sorted_lines:
            lines_by_sample.setdefault(line[0], []).append(line)
        for original_lines in lines_by_sample.values():

            # Audit finding 7: several alleles at one position of a sample (multi-allelic sites or
            # duplicated records) are one site for the IMD. The IMD is computed between distinct
            # sites and every allele of a site keeps the distance of its site; before, each allele
            # got an IMD of 0, which hid its distance to the neighbors and made the rows disappear.
            site_rows = {}
            for line in original_lines:
                site_rows.setdefault((line[0], line[2]), []).append(line)
            original_lines = [rows[0] for rows in site_rows.values()]

            def write_site(x, out_file):
                """Write the row of a site once per allele found at that position."""
                for allele in site_rows[(x[1], x[3])]:
                    row = list(x)
                    row[4], row[5] = allele[3], allele[4]
                    print("\t".join([str(y) for y in row]), file=out_file)

            # Ensure there are at least two mutations present for a given chromosome
            if len(original_lines) > 1:
                centro_start = centromeres[genome][original_lines[0][1]][
                    0
                ]  # Position of the start of the centromere
                centro_end = centromeres[genome][original_lines[0][1]][
                    1
                ]  # Position of the end of the centromere

                # Calculate the minimum IMDs ensuring that the samples, chromosome, and chromosome arms match
                # distances = [[int(y[2])-(int(x[2])-2+len(x[3])+len(x[4]))] + x + ['c'] if x[0] == y[0] and int(y[2])<centro_start else [int(y[2])-(int(x[2])-2+len(x[3])+len(x[4]))] + x +['c'] if x[0] == y[0] and centro_end<int(x[2]) else 'aa' if x[0] != y[0] else [int(x[2]) - (int(original_lines[original_lines.index(x)-1][2])-2+len(original_lines[original_lines.index(x)-1][3])+len(original_lines[original_lines.index(x)-1][4]))] + x +['d'] for x,y in zip(original_lines, original_lines[1:])]
                distances = [
                    (
                        [int(y[2]) - (int(x[2]) - 2 + len(x[3]) + len(x[4]))] + x + ["c"]
                        if x[0] == y[0] and int(y[2]) < centro_start
                        else (
                            [int(y[2]) - (int(x[2]) - 2 + len(x[3]) + len(x[4]))]
                            + x
                            + ["c"]
                            if x[0] == y[0] and centro_end < int(x[2])
                            else [
                                abs(
                                    int(x[2])
                                    - (
                                        int(original_lines[original_lines.index(x) - 1][2])
                                        - 2
                                        + len(
                                            original_lines[original_lines.index(x) - 1][3]
                                        )
                                        + len(
                                            original_lines[original_lines.index(x) - 1][4]
                                        )
                                    )
                                )
                            ]
                            + x
                            + ["d"]
                        )
                    )
                    for x, y in zip(original_lines, original_lines[1:])
                ]
                distances_new = distances[:] + [
                    [distances[-1][0]] + original_lines[-1] + [distances[-1][-1]]
                ]
                # final_distances = ['bb' if x[1] != y[1] else [min(x[0],y[0])]+y[1:-1]+[x[0]] if x[0] != 'a' and y[0] != 'a' and x[-1] != 'd' else [y[0]] + y[1:-1]+[y[0]] if x[-1] == 'd' and x[0] != 'a' and y[0] != 'a' else x if y[0] == 'a' and x[0] != 'a' else y if x[0] == 'a' and y[0] != 'a' else 'bb' for x,y in zip(distances_new[:], distances_new[1:])]
                final_distances = [
                    (
                        [x[0]] + x[1:-1] + [x[0]]
                        if x[1] != y[1] and y[1] == "a"
                        else (
                            [min(x[0], y[0])] + y[1:-1] + [x[0]]
                            if x[0] != "a" and y[0] != "a" and x[-1] != "d"
                            else (
                                [y[0]] + y[1:-1] + [y[0]]
                                if x[-1] == "d" and x[0] != "a" and y[0] != "a"
                                else (
                                    x
                                    if y[0] == "a" and x[0] != "a"
                                    else y if x[0] == "a" and y[0] != "a" else "bb"
                                )
                            )
                        )
                    )
                    for x, y in zip(distances_new[:], distances_new[1:])
                ]
                final_distances = [distances_new[0]] + final_distances
                final_distances_filtered = [
                    x for x in final_distances if x[0] != "b" and x[0] != "a"
                ]

                if (
                    len(final_distances_filtered) == 0
                ):  # Continue if there are no IMDs (ex. only two mutation are a chromosome and occured on different arms)
                    continue

                # Open the current sample file and start iterating through the final distances
                # If a new sample is reached, close the current file and open the new one
                # The files are written under the simulations/ folder with the following fiels:
                # <IMD> <Sample_sim#> <chrom> <pos> <ref> <alt> <plotIMD_for_rainfall plots>
                sample = final_distances_filtered[0][1]
                out_file = open(
                    output_path_original
                    + sample
                    + "_"
                    + file_context2
                    + "_intradistance.txt",
                    "a",
                    10000000,
                )
                for x in final_distances_filtered:
                    if x[1] == sample:
                        write_site(x, out_file)
                    else:
                        out_file.close()
                        sample = x[1]
                        out_file = open(
                            output_path_original
                            + sample
                            + "_"
                            + file_context2
                            + "_intradistance.txt",
                            "a",
                            10000000,
                        )
                        write_site(x, out_file)
                out_file.close()

            elif len(original_lines) == 1:
                # Handle the case when there is only one mutation
                sample = original_lines[0][0]
                out_file = (
                    output_path_original
                    + sample
                    + "_"
                    + file_context2
                    + "_intradistance.txt"
                )

                with open(out_file, "a", 10000000) as out_file:
                    # Add a new column with the value 1 million
                    # no neighbor: every allele of the site gets the singleton distance
                    for mutation in site_rows[(sample, original_lines[0][2])]:
                        print("\t".join([str(y) for y in [1000000] + mutation]), file=out_file)


def eventProbability(
    project, input_path, simulation_path, output_path_original, windowSize, k=10
):
    # Audit finding 11: EVENT_PROBABILITY is a product of KDE density heights, not a probability.
    warnings.warn(
        "EVENT_PROBABILITY is an experimental score (a product of kernel density estimates "
        "of the simulated mutations), not a probability; the model is under review.",
        UserWarning,
        stacklevel=2,
    )
    windowSize *= 10
    subclasses = ["DBS", "MBS", "omikli", "kataegis", "other"]

    samples = sorted(x for x in os.listdir(simulation_path) if x[0] != ".")
    for sample in samples:
        samplePath = os.path.join(simulation_path, sample)
        # A fixed seed over the sorted file names makes the choice of simulations reproducible
        simFiles = random.Random(EVENT_PROBABILITY_SEED).sample(
            sorted(x for x in os.listdir(samplePath) if x[0] != "."), k
        )
        simIMDs = {}

        for sim in simFiles:
            muts = pd.read_csv(
                os.path.join(samplePath, sim),
                sep="\t",
                header=None,
                names=["IMD", "sample", "chrom", "pos", "ref", "alt", "plotIMD"],
                dtype={"chrom": str},  # all-numeric names would be read as integers
            )
            chroms = list(set(muts["chrom"]))
            for chrom in chroms:
                chromMuts = muts[muts["chrom"] == chrom]
                if chrom not in simIMDs:
                    simIMDs[chrom] = {}

                maxDistance = max(chromMuts["pos"]) + windowSize
                slideSize = windowSize
                imdsProbabilities = []
                currentWindow = 0
                densities = []
                while currentWindow < maxDistance:
                    currentMuts = chromMuts[
                        (chromMuts["pos"] < currentWindow + slideSize)
                        & (chromMuts["pos"] >= currentWindow)
                    ]
                    if (
                        str(currentWindow) + ":" + str(currentWindow + slideSize)
                        not in simIMDs[chrom]
                    ):
                        simIMDs[chrom][
                            str(currentWindow) + ":" + str(currentWindow + slideSize)
                        ] = []
                    simIMDs[chrom][
                        str(currentWindow) + ":" + str(currentWindow + slideSize)
                    ] += list(currentMuts["pos"])
                    currentWindow += slideSize

        distributions = {}
        for chrom in simIMDs:
            distributions[chrom] = {}
            for subwindow in simIMDs[chrom]:
                if (
                    len(simIMDs[chrom][subwindow]) == 0
                    or len(set(simIMDs[chrom][subwindow])) <= 1
                ):
                    continue
                x, y = (
                    sns.kdeplot(
                        np.log10(simIMDs[chrom][subwindow]), bw_adjust=20, gridsize=200
                    )
                    .get_lines()[0]
                    .get_data()
                )
                distributions[chrom][subwindow] = [x, y]
                plt.pyplot.close()
        for subclass in subclasses:
            subclassPath = os.path.join(
                input_path, "output", "clustered", subclass, sample + ".vcf"
            )

            if os.path.exists(subclassPath):
                out = open(
                    os.path.join(
                        input_path,
                        "output",
                        "clustered",
                        subclass,
                        sample + "_prob.vcf",
                    ),
                    "w",
                )
                print(
                    "\t".join(
                        [
                            "#CHROM",
                            "POS",
                            "ID",
                            "REF",
                            "ALT",
                            "QUAL",
                            "FILTER",
                            "INFO",
                            "EVENT_PROBABILITY",
                        ]
                    ),
                    file=out,
                )
                clusteredMuts = pd.read_csv(subclassPath, header=0, sep="\t")
                clusteredMuts["EVENTS"] = [
                    x.split("groupNumber")[1] for x in clusteredMuts["INFO"]
                ]

                for chrom in sorted(set(clusteredMuts["#CHROM"]), key=str):
                    currentClusteredMuts = clusteredMuts[
                        clusteredMuts["#CHROM"] == chrom
                    ]
                    for event in sorted(
                        set(currentClusteredMuts["EVENTS"]),
                        key=lambda e: (int(e) if str(e).strip().isdigit() else float("inf"), str(e)),
                    ):
                        eventMuts = currentClusteredMuts[
                            currentClusteredMuts["EVENTS"] == event
                        ]
                        locations = list(eventMuts["POS"])
                        imds = np.log10(
                            [y - x for x, y in zip(locations, locations[1:])]
                        )
                        window = int(locations[0] / windowSize) * windowSize
                        window = str(window) + ":" + str(window + windowSize)
                        currentDist = distributions[str(chrom)][window]
                        probs = [
                            currentDist[1][bisect.bisect_left(currentDist[0], x)]
                            for x in imds
                        ]
                        probs = np.prod(np.array(probs))
                        if probs < constants.EVENT_PROBABILITY_FLOOR:
                            probs = constants.EVENT_PROBABILITY_FLOOR

                        for i, x in eventMuts.iterrows():
                            x = list(x)[:-1] + [probs]
                            print("\t".join([str(y) for y in x]), file=out)
                out.close()


def load_previous_analysis(project, ref_dir, contexts, correction, chrom_based):
    """
    Loads what a previous run of the hotspot analysis saved and the subclassification
    needs (analysis="subClassify"): the dense regions and the corrected IMD cutoffs of every
    sample. Without correction both are empty, as in a run that does not correct. Raises a
    FileNotFoundError that lists every missing file before anything is loaded.
    """
    data = ref_dir + "output/simulations/data/"
    suffix = "_chrom" if chrom_based else ""
    required = [data + "imds" + suffix + ".pickle"]
    if correction:
        required += [
            data + "imds_corrected" + suffix + ".pickle",
            data + "regions_corrected" + suffix + ".pickle",
        ]
    if not constants.is_indel_context(contexts):
        vcf_dir = "output/vcf_files" + ("_corrected" if correction else "")
        required.append(
            ref_dir + vcf_dir + "/" + project + "_clustered/SNV/" + project + "_clustered.txt"
        )
    missing = [path for path in required if not os.path.exists(path)]
    if missing:
        raise FileNotFoundError(
            'analysis="subClassify" needs the results of a previous run of the hotspot '
            'analysis (run analysis="hotspot" or analysis="all" first). Missing files:\n  '
            + "\n  ".join(missing)
        )
    if not correction:
        return {}, {}
    with open(data + "regions_corrected" + suffix + ".pickle", "rb") as handle:
        regions = pickle.load(handle)
    with open(data + "imds_corrected" + suffix + ".pickle", "rb") as handle:
        imds = pickle.load(handle)
    return regions, imds


def _run_analysis(
    project,
    genome,
    contexts,
    simContext,
    input_path,
    output_type="all",
    analysis="all",
    interdistance="96",
    exome=False,
    clustering_vaf=False,
    sortSims=True,
    extraction=False,
    correction=True,
    startProcess=1,
    endProcess=25,
    totalIterations=1000,
    calculateIMD=True,
    chrom_based=False,
    max_cpu=None,
    subClassify=False,
    variant_caller="standard",
    tumor_vaf_column = 10, 
    vaf_field = "AF",
    includedVAFs=True,
    includedCCFs=False,
    windowSize=1000000,
    bedRanges=None,
    plotIMDfigure=True,
    plotRainfall=True,
    probability=False,
):
    """
    Organizes all of the data structures and calls all of the sub-functions. This is the main function called when running SigProfilerClusters.

    Parameters:
                            project	->	user provided project name (string)
                             genome	->	the reference genome used for the given analysis (string)
                       contexts	->	the context used to calculate IMDs: "SBS96" or "ID83" (aliases: "96", "SBS", "ID", "83", "INDEL"; string; DBS is not an input)
                     simContext	->	the simulated context used for the background model (list of strings; ie ["6144"])
                     input_path	->	the directory for the given project. This should contain the output from SigProfilerMatrixGenerator/Simulator (string)
                    output_type	->	only 'all' is supported (the distances of the original samples and of the simulations are calculated together; default='all')
                       analysis	->	optional parameter for specifying which type of analysis to perform. Helpful for when the IMDs have already been calcuated. (i.e. "all", "hotspot", "all"; default="all")
              interdistance	->	the mutation types to calculate IMDs between: the same names as contexts, or "SI" (NOT STABLE; RECOMMENDED NOT TO USE; default='96')
             clustering_vaf	->	optional paramater to save the VAF from the original mutation files. Required if performing subclassification (default=False; switched to True if subclassify=True)
                       sortSims	->	optional parameter that requires the function to sort the simulations or not (boolean; default=True)
                     extraction	[DEPRACATED]->	optional parameter to perform a signature extraction directly after the hotspot analysis (boolean; default=True)
                     correction	->	optional parameter to perform a genome-wide mutational density correction (boolean; default=False)
               startProcess	[DEPRACATED] ->	the number of starting processes for performing signature extraction (int; default=1)
                     endProcess	[DEPRACATED] ->	the number of final processes for performing signature extraction (int; default=25)
            totalIterations	[DEPRACATED] ->	the number of iterations for performing signature extraction (int; default=1000)
               calculateIMD	->	optional parameter to calculate the IMDs. This will save time if you need to rerun the subclassification step only (boolean; default=True)
                    chrom_based	->	optional parameter to perform the analysis on a chromosome basis (boolean; default=False)
                            max_cpu	->	optional parameter to specify the number of maximum cpu's to use for parallelizing the code (integer; default=None: uses all available cpu's)
                    subClassify	->	optional parameter to subclassify the clustered mutations into refinded classes including DBSs, extended MBSs, kataegis, etc. (boolean; default=False)
                variant_caller	->	optional parameter that informs the tool of what format the VAF scores are provided (boolean; default=None). This is required when subClassify=True. Currently, there are four supported formats:
                                -> caveman: If your VAF is recorded in the 11th column of your VCF as the last number of the colon delimited values, set variant_caller="caveman".
                                -> standard: If your VAF is recorded in the 8th or 10th column of your VCF as VAF=xx or AF=xx, set variant_caller="standard".
                                -> mutect2: If your VAF is recorded in the 10th or 11th column of your VCF as AF=xx, set variant_caller="mutect2".
                                -> custom: Configure tumor_vaf_column as 0 based index (Default = 10) and the vaf_field tag (Default = AF)
            tumor_vaf_column -> Activated when variant_caller = "custom". 0 based index of the column in VCF containing VAF (Default 10)
            vaf_field       -> Activated when variant_caller = "custom". Set the VAF tag (Default AF)

               includedVAFs ->  optional parameter that informs the tool of the inclusion of VAFs in the dataset (boolean; default=True)
               includedCCFs ->  optional parameter that informs the tool of the inclusion of cancer cell fractions in the dataset (boolean; default=True)
                     windowSize	->	the size of the window used for correcting the IMDs based upon mutational density within a given genomic range (integer; default=10000000)
              plotIMDfigure	->	optional parameter that generates IMD and mutational spectra plots for each sample (boolean; default=True).
               plotRainfall	->	optional parameter that generates rainfall plots for each sample using the subclassification of clustered events (boolean; default=True).

    Returns:
            None

    Outputs:
            Filtered clustered and non-clustered mutations along with optional subclassification of these given mutations. The IMD plots are provided and rainfall plots are provided when
            subClassify=True. Note that the IMD plots do not have the resolution to show the corrected IMDs if correction=True. The corrected mutations are shown in the SBS96, but the
            distances are not used in the histogram.
    """
    if output_type not in SUPPORTED_OUTPUT_TYPES:
        raise InvalidParameterError(
            f"output_type {output_type!r} is not supported; valid values: "
            + ", ".join(repr(x) for x in SUPPORTED_OUTPUT_TYPES)
            + " (the distances of the original samples and of the simulations are "
            "always calculated together)"
        )

    if includedCCFs:
        includedVAFs = False

    if chrom_based:
        raise InvalidParameterError(
            "The parameter chrom_based has been deprecated within the SigProfilerClusters "
            "pipeline. Please rerun with chrom_based=False. Note, you may still run the "
            "simulations using chrom_based=True, which will be used as the background model "
            "for the clustered analysis."
        )
    validate_arguments(
        analysis, contexts, simContext, windowSize, variant_caller, interdistance
    )
    # From here on the contexts have their standard names (SBS96, ID83): they are part of the
    # names of the files that are written and read
    contexts = constants.canonical_context(contexts)
    if interdistance != "SI":
        interdistance = constants.canonical_context(interdistance)
    # Conversions for NCBI chromosome names to standard notation
    ncbi_chrom = {
        "NC_000067.6": "1",
        "NC_000068.7": "2",
        "NC_000069.6": "3",
        "NC_000070.6": "4",
        "NC_000071.6": "5",
        "NC_000072.6": "6",
        "NC_000073.6": "7",
        "NC_000074.6": "8",
        "NC_000075.6": "9",
        "NC_000076.6": "10",
        "NC_000077.6": "11",
        "NC_000078.6": "12",
        "NC_000079.6": "13",
        "NC_000080.6": "14",
        "NC_000081.6": "15",
        "NC_000082.6": "16",
        "NC_000083.6": "17",
        "NC_000084.6": "18",
        "NC_000085.6": "19",
        "NC_000086.7": "X",
        "NC_000087.7": "Y",
    }

    # Centromere locations for a given set of genomes for each chromosome
    # This needs to be updated for a new genome to work properly
    # Should incorporate this as an optional parameter to the tool, that allows the user to specify a centromere file
    centromeres = {
        "GRCh38": {
            "1": [122026460, 125184587],
            "10": [39686683, 41593521],
            "11": [51078349, 54425074],
            "12": [34769408, 37185252],
            "13": [16000001, 18051248],
            "14": [16000001, 18173523],
            "15": [17000001, 19725254],
            "16": [36311159, 38280682],
            "17": [22813680, 26885980],
            "18": [15460900, 20861206],
            "19": [24498981, 27190874],
            "2": [92188146, 94090557],
            "20": [26436233, 30038348],
            "21": [10864561, 12915808],
            "22": [12954789, 15054318],
            "3": [90772459, 93655574],
            "4": [49708101, 51743951],
            "5": [46485901, 50059807],
            "6": [58553889, 59829934],
            "7": [58169654, 60828234],
            "8": [44033745, 45877265],
            "9": [43236168, 45518558],
            "X": [58605580, 62412542],
            "Y": [10316945, 10544039],
        },
        "GRCh37": {
            "1": [121535434, 124535434],
            "10": [39254935, 42254935],
            "11": [51644205, 54644205],
            "12": [34856694, 37856694],
            "13": [16000000, 19000000],
            "14": [16000000, 19000000],
            "15": [17000000, 20000000],
            "16": [35335801, 38335801],
            "17": [22263006, 25263006],
            "18": [15460898, 18460898],
            "19": [24681782, 27681782],
            "2": [92326171, 95326171],
            "20": [26369569, 29369569],
            "21": [11288129, 14288129],
            "22": [13000000, 16000000],
            "3": [90504854, 93504854],
            "4": [49660117, 52660117],
            "5": [46405641, 49405641],
            "6": [58830166, 61830166],
            "7": [58054331, 61054331],
            "8": [43838887, 46838887],
            "9": [47367679, 50367679],
            "X": [58632012, 61632012],
            "Y": [10316945, 10544039],
        },
        "mm39": {
            "1": [110000, 3000000],
            "10": [110000, 3000000],
            "11": [110000, 3000000],
            "12": [110000, 3000000],
            "13": [110000, 3000000],
            "14": [110000, 3000000],
            "15": [110000, 3000000],
            "16": [110000, 3000000],
            "17": [110000, 3000000],
            "18": [110000, 3000000],
            "19": [110000, 3000000],
            "2": [110000, 3000000],
            "3": [110000, 3000000],
            "4": [110000, 3000000],
            "5": [110000, 3000000],
            "6": [110000, 3000000],
            "7": [110000, 3000000],
            "8": [110000, 3000000],
            "9": [110000, 3000000],
            "X": [110000, 3000000],
            "Y": [110000, 3000000],
        },
        "mm10": {
            "1": [110000, 3000000],
            "10": [110000, 3000000],
            "11": [110000, 3000000],
            "12": [110000, 3000000],
            "13": [110000, 3000000],
            "14": [110000, 3000000],
            "15": [110000, 3000000],
            "16": [110000, 3000000],
            "17": [110000, 3000000],
            "18": [110000, 3000000],
            "19": [110000, 3000000],
            "2": [110000, 3000000],
            "3": [110000, 3000000],
            "4": [110000, 3000000],
            "5": [110000, 3000000],
            "6": [110000, 3000000],
            "7": [110000, 3000000],
            "8": [110000, 3000000],
            "9": [110000, 3000000],
            "X": [110000, 3000000],
            "Y": [110000, 3000000],
        },
        "mm9": {
            "1": [0, 3000000],
            "10": [0, 3000000],
            "11": [0, 3000000],
            "12": [0, 3000000],
            "13": [0, 3000000],
            "14": [0, 3000000],
            "15": [0, 3000000],
            "16": [0, 3000000],
            "17": [0, 3000000],
            "18": [0, 3000000],
            "19": [0, 3000000],
            "2": [0, 3000000],
            "3": [0, 3000000],
            "4": [0, 3000000],
            "5": [0, 3000000],
            "6": [0, 3000000],
            "7": [0, 3000000],
            "8": [0, 3000000],
            "9": [0, 3000000],
            "X": [0, 3000000],
            "Y": [0, 3000000],
        },
        "rn6": {
            "1": [1000000000, 1000000000],
            "10": [1000000000, 1000000000],
            "11": [1000000000, 1000000000],
            "12": [1000000000, 1000000000],
            "13": [1000000000, 1000000000],
            "14": [1000000000, 1000000000],
            "15": [1000000000, 1000000000],
            "16": [1000000000, 1000000000],
            "17": [1000000000, 1000000000],
            "18": [1000000000, 1000000000],
            "19": [1000000000, 1000000000],
            "2": [1000000000, 1000000000],
            "20": [1000000000, 1000000000],
            "21": [1000000000, 1000000000],
            "22": [1000000000, 1000000000],
            "3": [1000000000, 1000000000],
            "4": [1000000000, 1000000000],
            "5": [1000000000, 1000000000],
            "6": [1000000000, 1000000000],
            "7": [1000000000, 1000000000],
            "8": [1000000000, 1000000000],
            "9": [1000000000, 1000000000],
            "X": [1000000000, 1000000000],
            "Y": [1000000000, 1000000000],
        },
    }

    # The centromere positions are needed by the distance calculation of every genome: fail now,
    # with the genomes that can be used, instead of with a KeyError after the work was done
    if genome not in centromeres:
        raise InvalidParameterError(
            "The genome "
            + repr(genome)
            + " is not supported. Supported genomes: "
            + ", ".join(sorted(centromeres))
            + "."
        )
    check_input_folder(input_path)

    # Collect the paths for the respective reference genome and quantiy the genome length and respective bins
    # using the default/provided window size
    chrom_path, ref_dir = matRef.reference_paths(genome)
    chromLengths = {}
    binsDensity = []
    if correction or subClassify:
        chroms = [
            x.split(".")[0]
            for x in os.listdir(chrom_path)
            if x != ".DS_Store"
            and x[0] != "G"
            and x[0] != "M"
            and x != "BED_" + genome + "_proportions.txt"
            and x[0] != "m"
            and "_proportions" not in x
        ]
        chroms = sorted(
            chroms,
            key=lambda x: [
                "1",
                "2",
                "3",
                "4",
                "5",
                "6",
                "7",
                "8",
                "9",
                "10",
                "11",
                "12",
                "13",
                "14",
                "15",
                "16",
                "17",
                "18",
                "19",
                "20",
                "21",
                "22",
                "X",
                "Y",
            ].index(x[0:]),
        )
        chromLengths[genome] = {}
        totalLength = 0
        for i in range(0, len(chroms), 1):
            with open(chrom_path + chroms[i] + ".txt", "rb") as f:
                chromLengths[genome][chroms[i]] = totalLength
                chrom_length = len(f.read())
                totalLength += chrom_length
        for i in range(0, totalLength, windowSize):
            binsDensity.append(i)
        binsDensity.append(totalLength)

    # Begin formatting the various parameters; especially for file naming conventions
    inter = ""
    inter_label = interdistance
    interdistance = False
    ref_dir = input_path
    # Determine the appropriate analysis/file extensions given the provided parameters
    if analysis == "hotspot" or analysis == "subClassify":
        output_type = None
    if inter_label == "SI":
        interdistance = True
    else:
        if constants.is_indel_context(inter_label):
            inter = "ID"
        else:
            inter = "SNP"
    simContext = sorted(simContext, reverse=True)
    file_context = "_".join(simContext)

    if ref_dir[-1] != "/":  # Ensure a trailing backslash to the given path
        ref_dir += "/"

    if analysis == "subClassify":
        # Reuse the results of a previous run; fails now, before anything is deleted,
        # if some of them is missing
        regions, imds = load_previous_analysis(
            project, ref_dir, contexts, correction, chrom_based
        )

    path_suffix = ""
    if correction:
        path_suffix += "_corrected"
    # Done formatting file extensions

    # Start formatting output folder structure and the log files
    output_path = ref_dir + "output/vcf_files" + path_suffix + "/single/"
    if os.path.exists(output_path) or output_type != None:
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
    os.makedirs(output_path)
    output_log_path = ref_dir + "logs/"
    if not os.path.exists(output_log_path):
        os.makedirs(output_log_path)

    time_stamp = log_stamp()

    log_file = (
        output_log_path
        + "SigProfilerClusters_"
        + project
        + "_"
        + genome
        + "_"
        + time_stamp
        + ".out"
    )
    error_file = (
        output_log_path
        + "SigProfilerClusters_"
        + project
        + "_"
        + genome
        + "_"
        + time_stamp
        + ".err"
    )
    if os.path.exists(error_file):
        os.remove(error_file)
    if os.path.exists(log_file):
        os.remove(log_file)
    sys.stderr = open(error_file, "w")
    log_out = open(log_file, "w")
    log_out.write("THIS FILE CONTAINS THE METADATA ABOUT SYSTEM AND RUNTIME\n\n\n")
    log_out.write("-------System Info-------\n")
    log_out.write(
        "Operating System Name: "
        + platform.uname()[0]
        + "\n"
        + "Nodename: "
        + platform.uname()[1]
        + "\n"
        + "Release: "
        + platform.uname()[2]
        + "\n"
        + "Version: "
        + platform.uname()[3]
        + "\n"
    )
    log_out.write("\n-------Python and Package Versions------- \n")
    log_out.write(
        "Python Version: "
        + str(platform.sys.version_info.major)
        + "."
        + str(platform.sys.version_info.minor)
        + "."
        + str(platform.sys.version_info.micro)
        + "\n"
    )
    log_out.write("SigProfilerMatrixGenerator Version: " + sig.__version__ + "\n")
    log_out.write("SigProfilerPlotting version: " + sigPlt.__version__ + "\n")
    log_out.write("matplotlib version: " + plt.__version__ + "\n")
    log_out.write("scipy version: " + scipy.__version__ + "\n")
    log_out.write("numpy version: " + np.__version__ + "\n")

    log_out.write("\n-------Vital Parameters Used for the execution -------\n")
    log_out.write(
        "Project: {}\nGenome: {}\nContext: {}\ninterdistance: {}\ninput_path: {}\noutput_type: {}\n".format(
            project, genome, simContext, interdistance, ref_dir, output_type
        )
    )
    log_out.write("\n-------Date and Time Data------- \n")
    tic = datetime.datetime.now()
    log_out.write(
        "Date and Clock time when the execution started: " + str(tic) + "\n\n\n"
    )
    log_out.write("-------Runtime Checkpoints------- \n")
    log_out.close()

    print("\n\n======================================", flush=True)
    print("Beginning SigProfilerClusters Analysis", flush=True)
    print("======================================\n\n", flush=True)
    # Log files are done begin created and the majority of the folder structures are set

    ############################################
    # Begin analysis ###########################
    ############################################
    start = time.time()
    context_ref = None
    # Only one type of analysis is specified
    if len(simContext) == 1:
        if simContext[0] == "INDEL" or simContext[0] == "ID":
            context_ref = "INDEL"
        else:
            context_ref = "SNV"

        original_vcf_path = classifyFunctions.get_input_vcf_path(ref_dir)
        vcf_path = original_vcf_path
        # Ignore hidden macOS files and subfolders
        vcf_files = [
            x
            for x in os.listdir(original_vcf_path)
            if x != ".DS_Store" and os.path.isfile(os.path.join(original_vcf_path, x))
        ]
        file_name = vcf_files[0].split(".")
        file_extension = file_name[-1]

        # Convert original input files (VCF-like files) into chromosome-based files based upon the input file format.
        # Calls upon a script used in the matrixGenerator
        if file_extension == "genome":
            convertIn.convertTxt(project, vcf_path, genome, output_path)
        else:
            if file_extension == "txt":
                snv, indel, skipped, samples = convertIn.convertTxt(
                    project, vcf_path, genome, output_path, ncbi_chrom, log_file
                )
            elif file_extension == "vcf":
                snv, indel, skipped, samples = convertIn.convertVCF(
                    project, vcf_path, genome, output_path, ncbi_chrom, log_file
                )
            elif file_extension == "maf":
                snv, indel, skipped, samples = convertIn.convertMAF(
                    project, vcf_path, genome, output_path, ncbi_chrom, log_file
                )
            elif file_extension == "tsv":
                snv, indel, skipped, samples = convertIn.convertICGC(
                    project, vcf_path, genome, output_path, ncbi_chrom, log_file
                )
            else:
                raise InvalidParameterError(
                    "The input file format of %r is not supported; the supported ones are: %s"
                    % (vcf_files[0], ", ".join(constants.INPUT_FILE_EXTENSIONS))
                )

    # More than one type of analysis is specified. Currently no longer supported
    else:
        raise InvalidParameterError(
            "Please only provide one type of context for analysis (simContext has "
            + str(len(simContext))
            + ")."
        )

    # Organize additional paths/files and corresponding suffixes
    file_context2 = inter_label
    if exome:
        file_context += "_exome"
        file_context2 += "_exome"
    output_path = (
        ref_dir
        + "output/simulations/"
        + project
        + "_intradistance_"
        + genome
        + "_"
        + file_context2
        + "/"
    )
    simulation_path_sorted = None
    original_samples = (
        ref_dir + "output/vcf_files" + path_suffix + "/single/" + context_ref + "/"
    )
    if sortSims:
        simulation_path = (
            ref_dir
            + "output/simulations/"
            + project
            + "_simulations_"
            + genome
            + "_"
            + file_context
            + "/"
        )
        simulation_path_sorted = (
            ref_dir
            + "output/simulations/"
            + project
            + "_simulations_"
            + genome
            + "_"
            + file_context
            + "_sorted/"
        )
        if os.path.exists(simulation_path_sorted):
            shutil.rmtree(simulation_path_sorted)
        os.makedirs(simulation_path_sorted)
    else:
        simulation_path = (
            ref_dir
            + "output/simulations/"
            + project
            + "_simulations_"
            + genome
            + "_"
            + file_context
            + "_sorted/"
        )
    output_path_original = (
        ref_dir
        + "output/simulations/"
        + project
        + "_intradistance_original_"
        + genome
        + "_"
        + file_context2
        + "/"
    )

    if analysis == "hotspot" and not os.path.isdir(output_path_original):
        # Results of an earlier version have "96"/"ID" in the names of the distance files
        legacy = [
            output_path_original.replace(
                "_" + genome + "_" + file_context2 + "/",
                "_" + genome + "_" + token + ("_exome" if exome else "") + "/",
            )
            for token in ("96", "ID", "INDEL")
        ]
        if any(os.path.isdir(path) for path in legacy):
            raise MissingInputError(
                "The distances of this project were saved by an earlier version, with the "
                "names '96' or 'ID' instead of 'SBS96' or 'ID83'. Run analysis=\"all\" again "
                "to recalculate them: " + output_path_original
            )

    # Checks for simulated data. If there are not simulations, then inform the user to generate simulations first before
    # performing the clusters analysis:
    if not os.path.exists(simulation_path):
        raise MissingSimulationsError(
            "There are no simulated data present for this project. Please generate simulations before running SigProfilerClusters.\n"
            "\tThe package can be installed via pip:\n\t\t\t$ pip install SigProfilerSimulator\n"
            "\n\tand used within a python3 sessions as follows:\n\t\t\t$ python3\n\t\t\t>> from SigProfilerSimulator import SigProfilerSimulator as sigSim\n"
            "\t\t\t>> sigSim.SigProfilerSimulator(project, project_path, genome, contexts=['6144'], simulations=100)\n\n"
            "\tFor a complete list of parameters, visit the github repo (https://github.com/AlexandrovLab/SigProfilerSimulator) or the documentation page (https://osf.io/usxjz/wiki/home/)"
        )
    if len(os.listdir(simulation_path)) < 100:
        raise MissingSimulationsError(
            "Please simulate a minimum of 100 simulations per sample to successfully run SigProfilerClusters."
        )

    # Collect simulation files and distribute them into the available processors for parallelization
    simulations = os.listdir(simulation_path)
    if ".DS_Store" in simulations:
        simulations.remove(".DS_Store")

    # Assign max_seed as all available or user-provided number of processors for parallelization
    if max_cpu:
        processors = max_cpu
    else:
        processors = mp.cpu_count()
    max_seed = processors

    # Calculate the mutational distances for the original samples (distance_one_file())
    # and for the simulated samples (distance_multiple_files_sims()). These series of if/else
    # statements were originally used to analyze simulations or original data separately, however,
    # we always analyze both at consecutively. The only relevant section is now output_type=='all'.
    # If output_type=='hotspot' or output_type=='subClassify', this portion of code is skipped, which
    # saves time if the downstream analysis needs to be rerun.
    if output_type != None:
        print("Calculating mutational distances...", end="", flush=True)
    if output_type == "all":
        # Remove old path iterations if they exist
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.makedirs(output_path)

        # Calculate distances for the original samples
        # Parallelize and calculate distances for the simulated samples
        # Organize the output paths (ie remove if they exist, and then create and emtpy folder)
        if os.path.exists(output_path_original):
            shutil.rmtree(output_path_original)
            os.makedirs(output_path_original)
        else:
            os.makedirs(output_path_original)
        original_samples_length = os.listdir(original_samples)
        if processors > len(original_samples_length):
            max_seed = len(original_samples_length)
        pool = mp.get_context("fork").Pool(max_seed)
        # pool_break = len(simulations)/max_seed
        samples_parallel = [[] for i in range(max_seed)]
        pool_bin = 0
        for chromFile in original_samples_length:
            if pool_bin == max_seed:
                pool_bin = 0
            samples_parallel[pool_bin].append(chromFile)
            pool_bin += 1
        results = []
        for i in range(0, len(samples_parallel), 1):
            r = pool.apply_async(
                distance_one_file,
                args=(
                    samples_parallel[i],
                    original_samples,
                    output_path_original,
                    file_context2,
                    genome,
                    centromeres,
                ),
            )
            results.append(r)
        pool.close()
        pool.join()
        for r in results:
            r.wait()
            if not r.successful():
                # Raises an error when not successful
                r.get()

        # The workers appended to the same file of each sample in the order in which they
        # finished: put its rows in a fixed order (audit finding 121)
        sort_distance_files(output_path_original)

        # Parallelize and calculate distances for the simulated samples
        if processors > len(simulations):
            max_seed = len(simulations)
        pool = mp.get_context("fork").Pool(max_seed)
        # sim_break = len(simulations)/max_seed
        simulations_parallel = [[] for i in range(max_seed)]
        sim_bin = 0
        for sim in simulations:
            if sim_bin == max_seed:
                sim_bin = 0
            simulations_parallel[sim_bin].append(sim)
            sim_bin += 1
        results = []
        for i in range(0, len(simulations_parallel), 1):
            r = pool.apply_async(
                distance_multiple_files_sims,
                args=(
                    output_path,
                    simulations_parallel[i],
                    simulation_path,
                    simulation_path_sorted,
                    file_context2,
                    genome,
                    centromeres,
                    sortSims,
                ),
            )
            results.append(r)
        pool.close()
        pool.join()
        for r in results:
            r.wait()
            if not r.successful():
                # Raises an error when not successful
                r.get()

    if output_type != None:
        print("Completed!", flush=True)  # IMD calculations have successfully completed.

    # Begin the hotspot analysis. This function will also produce the IMD distribution and spectra plots
    skipped_samples = {}
    skipped_table = ref_dir + "output/skipped_samples" + path_suffix + ".txt"
    if analysis == "hotspot" or analysis == "all":
        original = True
        firstRun = True
        signature = False
        percentage = False
        regions, imds = hotspot.hotSpotAnalysis(
            project,
            genome,
            contexts,
            simContext,
            ref_dir,
            windowSize,
            processors,
            plotIMDfigure,
            exome,
            chromLengths,
            binsDensity,
            original,
            signature,
            percentage,
            firstRun,
            clustering_vaf,
            calculateIMD,
            chrom_based,
            correction,
        )
        skipped_samples = dict(getattr(hotspot, "skipped_samplesFinal", {}))
        if regions is None:
            # No clustered mutations: not an error, but there is nothing to subclassify or plot
            print(
                "No clustered mutations were found. Skipping the subclassification and the plots."
            )
            hotspot.report_skipped_samples(skipped_samples, skipped_table, log_file)
            return

    # Allows for the clustered muutations to automactically run through the extraction.
    # This option is not recommended, but rather the extractions should be run separately after
    # the clustered analysis is completed.
    # if extraction:
    # 	print("Beginning signature extraction...")
    # 	sigs.sigProfilerExtractor("table", ref_dir+"output/extraction_clustered/", ref_dir+"output/vcf_files" + path_suffix + "/"+project+"_clustered/SNV/output/SBS/"+project+"_clustered.SBS96.all", genome, startProcess=startProcess, endProcess=endProcess, totalIterations=totalIterations)#, totalIterations=totalIterations)

    # Subclassify the clustered partition of mutations with the exception of indels. This function will also
    # produce the rainfall plots.
    if subClassify:
        if not constants.is_indel_context(contexts):
            print("Beginning subclassification of clustered mutations...", end="")
            if includedVAFs:
                classifyFunctions.pullVaf(
                    project, ref_dir, variant_caller, tumor_vaf_column, vaf_field, correction
                )
                classifyFunctions.findClustersOfClusters(
                    project,
                    chrom_based,
                    ref_dir,
                    windowSize,
                    chromLengths,
                    regions,
                    log_file,
                    genome,
                    processors,
                    imds,
                    correction,
                    includedCCFs,
                )
            elif includedCCFs:
                classifyFunctions.pullCCF(project, ref_dir, correction)
                classifyFunctions.findClustersOfClusters(
                    project,
                    chrom_based,
                    ref_dir,
                    windowSize,
                    chromLengths,
                    regions,
                    log_file,
                    genome,
                    processors,
                    imds,
                    correction,
                    includedCCFs,
                )
            else:
                classifyFunctions.findClustersOfClusters_noVAF(
                    project,
                    chrom_based,
                    ref_dir,
                    windowSize,
                    chromLengths,
                    regions,
                    log_file,
                    genome,
                    processors,
                    imds,
                    correction,
                )
            print("Completed!", flush=True)

        # Generate the rainfall plots.
        if plotRainfall:
            print("Generating a rainfall plot for all samples...", end="")
            plottingFunctions.rainfall(
                chrom_based,
                project,
                ref_dir,
                chrom_path,
                chromLengths,
                centromeres,
                contexts,
                includedVAFs,
                includedCCFs,
                correction,
                windowSize,
                bedRanges,
            )
            print("done")
            print("Subclassification of clustered mutations has finished!")

        # Generates output paths for final results saved as VCF files
        convertToVCF.generateAllPaths(ref_dir, contexts)
        convertToVCF.convertFiles(ref_dir, contexts, project, correction)

        # Add probability to all clustered mutations using same window sized-bins
        if probability:
            print("\n")
            print(
                "Calculating the probability of observing each clustered event...",
                end="",
                flush=True,
            )
            eventProbability(
                project, ref_dir, output_path, output_path_original, windowSize, k=10
            )
            print("done")

    end = time.time() - start
    hotspot.report_skipped_samples(skipped_samples, skipped_table, log_file)
    print(
        "SigProfilerClusters successfully finished! Elapsed time: "
        + str(round(end, 2))
        + " seconds."
    )


@functools.wraps(_run_analysis)
def analysis(*args, **kwargs):
    # _run_analysis redirects sys.stderr to the .err file of the run. Whatever happens
    # (success, error or exit), the original stream is given back and the file closed.
    original_stderr = sys.stderr
    try:
        return _run_analysis(*args, **kwargs)
    finally:
        redirected = sys.stderr
        sys.stderr = original_stderr
        if redirected is not original_stderr:
            redirected.close()
