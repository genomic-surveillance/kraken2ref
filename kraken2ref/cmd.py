import argparse
import os
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
        version='kraken2ref ' + __version__)

    parser.add_argument(
        '-s', '--sample_id',
        type = str,
        required = True,
        help = "Sample ID. [str]")

    subparsers = parser.add_subparsers(title="subcommands", help='kraken2ref sub-commands', dest='mode')
    report_parser = subparsers.add_parser("parse_report")

    report_parser.add_argument(
        '-i', '--in_file',
        type = str,
        required = True,
        help = "The kraken2 taxonomy report (typically 'report.txt') to process. [str/pathlike]")

    report_parser.add_argument(
        '-o', '--outdir',
        type = str,
        required = True,
        help = "Full path to output directory. [str/pathlike]")


    report_parser.add_argument(
        '-t', '--min_read_threshold',
        type = int,
        required = False,
        default = 100,
        help = "The absolute minimum number of reads to use as threshold; taxa with fewer reads assigned to them will not be considered. [int][Default = 100]")

    sort_read_parser = subparsers.add_parser("sort_reads")

    sort_read_parser.add_argument(
        "-fq1", "--fastq1",
        type = str,
        required = True,
        help = "First FASTQ file of paired end reads. [str/pathlike]")

    sort_read_parser.add_argument(
        "-fq2", "--fastq2",
        type = str,
        required = True,
        help = "Second FASTQ file of paired end reads. [str/pathlike]")

    sort_read_parser.add_argument(
        "-k", "--kraken_out",
        type = str,
        required = True,
        help = "Kraken2 output containing read-to-taxon mapping (typically output.kraken). [str/pathlike]")

    sort_read_parser.add_argument(
        "-r", "--ref_json",
        type = str,
        required = True,
        help = "Output JSON created by `kraken2ref parse_report`. [str/pathlike]")

    sort_read_parser.add_argument(
        "-u", "--update",
        action = "store_true",
        required = False,
        help = "Whether to update the kraken2ref JSON inplace or create a new updated copy. [switch]")

    return parser

def main():

    args = args_parser().parse_args()

    if args.mode == "parse_report":
        outdir = args.outdir
    if args.mode == "sort_reads":
        outdir = os.path.dirname(os.path.abspath(args.ref_json))

    if not os.path.exists(outdir):
            os.mkdir(outdir)
    # else:
    #     logging.warning(msg = f"CMD: OutputPathExists: The path {outdir} already exists; existing outputs may be overwritten.\n")
        # sys.stderr.write(f"OutputPathExists: The path {outdir} already exists; existing outputs may be overwritten.\n")


    logfile = os.path.join(outdir, f"{args.sample_id}_kraken2ref.log")
    logging.basicConfig(format='%(asctime)s | %(levelname)s | %(module)s - %(funcName)s | %(message)s', level=logging.NOTSET, datefmt="%Y-%m-%d %H:%M:%S", filename = logfile)

    tax_report = KrakenTaxonomyReport(sample_id = args.sample_id)
    if args.mode == "parse_report":
        tax_report.pick_reference_taxid(in_file = args.in_file, outdir = args.outdir, min_abs_reads = args.min_read_threshold)
    if args.mode == "sort_reads":
        tax_report.sort_reads_by_ref(sample_id = args.sample_id, fq1 = args.fastq1, fq2 = args.fastq2, kraken_out = args.kraken_out, ref_data = args.ref_json, update_output = args.update)

    ## for dev purposes
    # for k in graph_meta.keys():
    #     print(f"target = {k}: {graph_meta[k]}\n\n")

if __name__ == "__main__":
    exit(main())