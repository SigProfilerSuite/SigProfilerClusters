import argparse
from typing import List
from SigProfilerClusters import SigProfilerClusters


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")

def str2list(arg):
    return arg.split(",")

def parse_arguments_clusters(args: List[str], description: str) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=description)

    # Obligatorios
    parser.add_argument("project", help="Project name")
    parser.add_argument("genome", help="Reference genome")
    parser.add_argument("contexts", help="Mutational contexts")
    parser.add_argument("simContext", type=str2list, help="Simulated contexts (comma-separated)")
    parser.add_argument("input_path", help="Path to input directory")

    # Opcionales
    parser.add_argument("--output_type", default="all")
    parser.add_argument("--analysis", default="all")
    parser.add_argument("--interdistance", default="96")
    parser.add_argument("--exome", type=str2bool, nargs="?", const=True, default=False)
    parser.add_argument("--clustering_vaf", type=str2bool, nargs="?", const=True, default=False)
    parser.add_argument("--sortSims", type=str2bool, nargs="?", const=True, default=True)
    parser.add_argument("--extraction", type=str2bool, nargs="?", const=True, default=False)
    parser.add_argument("--correction", type=str2bool, nargs="?", const=True, default=True)
    parser.add_argument("--startProcess", type=int, default=1)
    parser.add_argument("--endProcess", type=int, default=25)
    parser.add_argument("--totalIterations", type=int, default=1000)
    parser.add_argument("--calculateIMD", type=str2bool, nargs="?", const=True, default=True)
    parser.add_argument("--chrom_based", type=str2bool, nargs="?", const=True, default=False)
    parser.add_argument("--max_cpu", type=int, default=None)
    parser.add_argument("--subClassify", type=str2bool, nargs="?", const=True, default=False)
    parser.add_argument("--variant_caller", default="standard")
    parser.add_argument("--includedVAFs", type=str2bool, nargs="?", const=True, default=True)
    parser.add_argument("--includedCCFs", type=str2bool, nargs="?", const=True, default=False)
    parser.add_argument("--windowSize", type=int, default=1000000)
    parser.add_argument("--bedRanges", default=None)
    parser.add_argument("--plotIMDfigure", type=str2bool, nargs="?", const=True, default=True)
    parser.add_argument("--plotRainfall", type=str2bool, nargs="?", const=True, default=True)
    parser.add_argument("--probability", type=str2bool, nargs="?", const=True, default=False)

    return parser.parse_args(args)


class CliController:
    def dispatch_sigProfilerClusters(self, user_args: List[str]) -> None:
        parsed_args = parse_arguments_clusters(
            user_args, "Cluster mutation analysis"
        )
        SigProfilerClusters.analysis(
            project=parsed_args.project,
            genome=parsed_args.genome,
            contexts=parsed_args.contexts,
            simContext=parsed_args.simContext,
            input_path=parsed_args.input_path,
            output_type=parsed_args.output_type,
            analysis=parsed_args.analysis,
            interdistance=parsed_args.interdistance,
            exome=parsed_args.exome,
            clustering_vaf=parsed_args.clustering_vaf,
            sortSims=parsed_args.sortSims,
            extraction=parsed_args.extraction,
            correction=parsed_args.correction,
            startProcess=parsed_args.startProcess,
            endProcess=parsed_args.endProcess,
            totalIterations=parsed_args.totalIterations,
            calculateIMD=parsed_args.calculateIMD,
            chrom_based=parsed_args.chrom_based,
            max_cpu=parsed_args.max_cpu,
            subClassify=parsed_args.subClassify,
            variant_caller=parsed_args.variant_caller,
            includedVAFs=parsed_args.includedVAFs,
            includedCCFs=parsed_args.includedCCFs,
            windowSize=parsed_args.windowSize,
            bedRanges=parsed_args.bedRanges,
            plotIMDfigure=parsed_args.plotIMDfigure,
            plotRainfall=parsed_args.plotRainfall,
            probability=parsed_args.probability,
        )