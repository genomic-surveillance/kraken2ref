## python imports
import os
import logging
import argparse

## import driver module
from kraken2ref.kraken2reference import KrakenProcessor

## collect version
try:
    from .version import version as __version__
    from .version import version_tuple
except ImportError:
    __version__ = "local/dev"
    version_tuple = (0, 0, "local/dev")

def collect_args():
    """Collect args
    """
    parser = argparse.ArgumentParser(
        description = "kraken2ref: extract reads to reference sequences from kraken2 outputs")

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='kraken2ref ' + __version__)

    parser.add_argument(
        '-s', '--sample_id',
        type = str,
        required = True,
        help = "Sample ID. [str]")

    parser.add_argument(
        '-i', '--in_file',
        type = str,
        required = True,
        help = "The kraken2 taxonomy report (typically 'report.txt') to process. [str/pathlike]")

    parser.add_argument(
        '-o', '--outdir',
        type = str,
        required = True,
        help = "Full path to output directory. [str/pathlike]")

    parser.add_argument(
        '-x', '--suffix',
        type = str,
        required = False,
        default="decomposed",
        help = """Suffix to add to output JSON file.
                    Format will be: <sample_id>_<suffix>.json. [str]""")

    parser.add_argument(
        '-t', '--min_read_threshold',
        type = int,
        required = False,
        default = 100,
        help = "The absolute minimum number of reads to use as threshold; taxa with fewer reads assigned to them will not be considered. [int][Default = 100]")

    parser.add_argument(
        '-m', '--poll_method',
        type = str,
        required = False,
        default = "kmeans",
        help = """Which polling method to use. [str] [Default = 'kmeans']
                    Valid choices: ['skew', 'kmeans', 'tiles']""")

    parser.add_argument(
        "-q", "--quiet",
        action = "store_true",
        required = False,
        help = "If specified, no stdderr will be run. [switch]")

    args = parser.parse_args()
    return args

def main():
    """User-facing driver function
    """
    ## collect args
    args = collect_args()

    ## explicitly set up outdir
    fixed_outdir = os.path.abspath(args.outdir)

    ## set up logging to file
    logfile = os.path.join(fixed_outdir, f"{args.sample_id}_kraken2ref.log")
    logging.basicConfig(format='%(asctime)s | %(levelname)s | %(module)s - %(funcName)s | %(message)s', level=logging.NOTSET, datefmt="%Y-%m-%d %H:%M:%S", filename = logfile)

    ## initialise driver module and analyse report
    my_processor = KrakenProcessor(args.sample_id)
    my_processor.analyse_report(input_kraken_report_file = args.in_file, input_threshold = args.min_read_threshold, input_method = args.poll_method, quiet=args.quiet)
    my_processor.write_output(prefix = fixed_outdir, suffix = args.suffix)

if __name__ == "main":
    main()

