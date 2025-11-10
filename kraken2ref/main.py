## python imports
import os, sys
import logging
import argparse

## import driver module
from kraken2ref.kraken2reference import KrakenProcessor
from kraken2ref.sort_reads import sort_reads_by_tax
from kraken2ref.dump_fastqs import dump_fastqs
import io

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

    subparsers = parser.add_subparsers(title="subcommands", help="kraken2ref sub-commands", dest="run_mode")
    report_parser = subparsers.add_parser("parse_report")

    report_parser.add_argument(
        '-i', '--in_file',
        type = str,
        required = True,
        help = "The kraken2 taxonomy report (typically 'report.txt') to process. [str/pathlike]")

    report_parser.add_argument(
        '-o', '--outdir',
        type = str,
        default=os.getcwd(),
        required = True,
        help = "Full path to output directory. [str/pathlike] (default = current working dir)")

    report_parser.add_argument(
        '-x', '--suffix',
        type = str,
        required = False,
        default="decomposed",
        help = """Suffix to add to output JSON file.
                    Format will be: <sample_id>_<suffix>.json. [str]""")

    report_parser.add_argument(
        '-t', '--min_read_threshold',
        type = int,
        required = False,
        default = 100,
        help = "The absolute minimum number of reads to use as threshold; taxa with fewer reads assigned to them will not be considered. [int][Default = 100]")

    report_parser.add_argument(
        '-m', '--poll_method',
        type = str,
        required = False,
        default = "max",
        help = """Which polling method to use. [str] [Default = 'max']
                    Valid choices: ['max', 'skew', 'kmeans', 'tiles']""")

    report_parser.add_argument(
        "-q", "--quiet",
        action = "store_true",
        required = False,
        help = "If specified, no stdderr will be run. [switch]")

    sort_reads_parser = subparsers.add_parser("sort_reads")

    sort_reads_parser.add_argument(
            '-t', '--taxon_list',
            required = False,
            type = str,
            help = "Comma-separated list of taxa to extract reads for. Eg: taxID_1,taxID_2... [str/pathlike]")

    sort_reads_parser.add_argument(
            "-k", "--kraken_out",
            type = str,
            required = True,
            help = "Path to kraken2 output file. [str/pathlike]")

    sort_reads_parser.add_argument(
            "-r", "--ref_json",
            type = str,
            required = False,
            help = "Output JSON created by `kraken2ref parse_report`. [str/pathlike]")

    sort_reads_parser.add_argument(
            "-m", "--mode",
            type = str,
            required = False,
            default="unique",
            help = "Which mode to use while sorting reads ['unique', 'tree']")

    sort_reads_parser.add_argument(
            "-c", "--condense",
            action = "store_true",
            required = False,
            help = "Whether to condense the outputs by root taxid. [switch]")

    sort_reads_parser.add_argument(
        "-u", "--update",
        action = "store_true",
        required = False,
        help = "Whether to update the kraken2ref JSON inplace or create a new updated copy. [switch]")

    sort_reads_parser.add_argument(
        '-o', '--outdir',
        type = str,
        default=os.getcwd(),
        required = False,
        help = "Full path to output directory. [str/pathlike] (default = current working dir)")


    dump_fqs_parser = subparsers.add_parser("dump_fastqs")

    dump_fqs_parser.add_argument(
        "--tax_to_readsid_path",
        type = str,
        required = True,
        help="json file containing tax to reads id (output by 'sort_to_reads' mode) [str/pathlike]"
    )

    dump_fqs_parser.add_argument(
        "-fq1", "--fastq1",
        type = str,
        required = True,
        help = "First FASTQ file of paired end reads. [str/pathlike]")

    dump_fqs_parser.add_argument(
        "-fq2", "--fastq2",
        type = str,
        required = True,
        help = "Second FASTQ file of paired end reads. [str/pathlike]")

    dump_fqs_parser.add_argument(
        '-o', '--outdir',
        type = str,
        default=os.getcwd(),
        required = False,
        help = "Full path to output directory. [str/pathlike] (default = current working dir)")

    dump_fqs_parser.add_argument(
        '--chunk_size',
        type = int,
        required = False,
        default=100_000,
        help = "Number of reads loaded into memory to process per batch [int] (default = 100_000)")

    dump_fqs_parser.add_argument(
        '--buffer_size',
        type = int,
        required = False,
        default=io.DEFAULT_BUFFER_SIZE,
        help = "buffer for writing output fq files size in bytes [int] (default=IO default buffer size)")

    dump_fqs_parser.add_argument(
        "-r", "--ref_json",
        type = str,
        required = False,
        help = "Output JSON created by `kraken2ref parse_report`. [str/pathlike]")

    dump_fqs_parser.add_argument(
        '--max_threads',
        type = int,
        required = False,
        default=1,
        help = "number of threads to provide to index [int] (default=1)")

    args = parser.parse_args()
    return args

# instantiating the decorator
#@profile
def main():
    """User-facing driver function
    """
    ## collect args
    args = collect_args()

    ## explicitly set up outdir
    fixed_outdir = os.path.abspath(args.outdir)
    if not os.path.exists(fixed_outdir):
        os.mkdir(fixed_outdir)

    ## set up logging to file
    logfile = os.path.join(fixed_outdir, f"{args.sample_id}_kraken2ref.log")
    logging.basicConfig(format='%(asctime)s | %(levelname)s | %(module)s - %(funcName)s | %(message)s', level=logging.NOTSET, datefmt="%Y-%m-%d %H:%M:%S", filename = logfile)

    if args.run_mode == "parse_report":

        ## initialise driver module and analyse report
        my_processor = KrakenProcessor(args.sample_id)
        my_processor.analyse_report(input_kraken_report_file = args.in_file, input_threshold = args.min_read_threshold, input_method = args.poll_method, quiet=args.quiet)
        my_processor.write_output(prefix = fixed_outdir, suffix = args.suffix)

    if args.run_mode == "sort_reads":
        sort_reads_by_tax(args)

    if args.run_mode == "dump_fastqs":
        dump_fastqs(args)

if __name__ == "main":
    main()