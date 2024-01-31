import argparse
import os
import shutil
import logging
from kraken2ref.src.kraken_taxonomy_report import KrakenTaxonomyReport


# get the version number from a file that is created by setuptools_scm
# when the package is installed.
try:
    from .version import version as __version__
    from .version import version_tuple
except ImportError:
    __version__ = "unknown version"
    version_tuple = (0, 0, "unknown version")

# NOTE about imports and using this script in develpoment and after production.
# The script uses absolute imports that will only work after being installed by pip.
# To run the script during development, run as follows from top level directory of the repo:
# python -m kraken_flu.cmd { OPTIONS }

def args_parser():
    """
    Command line argument parser
    """
    parser = argparse.ArgumentParser(
        description = "kraken2ref: extract reads to reference sequences from kraken2 outputs")

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='kraken2ref ' + __version__ )

    parser.add_argument(
        '-i', '--in_file',
        type = str,
        required = True,
        help = "The kraken2 taxonomy report (typically 'report.txt') to process. [str/pathlike]")

    parser.add_argument(
        '-t', '--min_read_threshold',
        type = int,
        required = False,
        default = 5,
        help = "The absolute minimum number of reads to use as threshold; taxa with fewer reads assigned to them will not be considered. [int][Default = 5]")

    return parser

def main():

    args = args_parser().parse_args()

    my_tax_report = KrakenTaxonomyReport(in_file = args.in_file, min_abs_reads = args.min_read_threshold)
    graph_meta = my_tax_report.pick_reference_taxid()

    ## for dev purposes
    for k in graph_meta.keys():
        print(f"target = {k}: {graph_meta[k]}\n\n")

if __name__ == "__main__":
    exit(main())