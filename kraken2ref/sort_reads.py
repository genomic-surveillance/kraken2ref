import os, sys
import json
import pandas as pd
import argparse
from Bio import SeqIO
from concurrent import futures
import datetime
import logging

def collect_args():
    """Function to collect command-lne arguments.

    Returns:
        args (Namespace): Argparse namespace
    """
    msg = f"""
        Extract reads to reference sequences from kraken2 outputs.
        This script splits a pair of fastq files based on the kraken2 output, and creates a set of <sample_id>_<taxon_id_R{1,2}.fq.
        Supports three modes: unique, tree, and condensed.
        Unique mode: Extract ONLY reads that are uniquely assigned to the specified taxon ID(s). Specify taxon IDs as follows: `-m unique -t taxon1[,taxon2,taxon3...]`.
                        NB: `-r path/to/kraken2ref.json` can be used with unique mode but is not required.
        Tree mode: Extract ALL reads in the taxonomy tree of the specified taxon ID(s). Usage: `-m tree -r path/to/kraken2ref.json`
        Condensed mode: Builds on tree mode; produce one set of fastq files per species, RATHER THAN per refernce. Usage: `-m tree -r path/to/kraken2ref.json -c`
        """
    parser = argparse.ArgumentParser(
            description = msg)

    parser.add_argument(
            '-s', '--sample_id',
            required = True,
            type = str,
            help = "Sample ID [str]")

    parser.add_argument(
            '-t', '--taxon_list',
            required = False,
            type = str,
            help = "Comma-separated list of taxa to extract reads for. Eg: taxID_1,taxID_2... [str/pathlike]")

    parser.add_argument(
            "-k", "--kraken_out",
            type = str,
            required = True,
            help = "Path to kraken2 output file. [str/pathlike]")

    parser.add_argument(
            "-fq1", "--fastq1",
            type = str,
            required = True,
            help = "First FASTQ file of paired end reads. [str/pathlike]")

    parser.add_argument(
            "-fq2", "--fastq2",
            type = str,
            required = True,
            help = "Second FASTQ file of paired end reads. [str/pathlike]")

    parser.add_argument(
            "-r", "--ref_json",
            type = str,
            required = False,
            help = "Output JSON created by `kraken2ref parse_report`. [str/pathlike]")

    parser.add_argument(
            "-m", "--mode",
            type = str,
            required = False,
            default="unique",
            help = "Which mode to use while sorting reads ['unique', 'tree']")

    parser.add_argument(
            "-c", "--condense",
            action = "store_true",
            required = False,
            help = "Whether to condense the outputs by root taxid. [switch]")

    parser.add_argument(
        "-u", "--update",
        action = "store_true",
        required = False,
        help = "Whether to update the kraken2ref JSON inplace or create a new updated copy. [switch]")

    parser.add_argument(
        '-o', '--outdir',
        type = str,
        required = False,
        help = "Full path to output directory. [str/pathlike]")


    args = parser.parse_args()

    return args

def dump_to_file(sample_id, tax_to_readids_dict, fq1, fq2, outdir):
    """Function that dumps reads to file.

    Args:
        tax_to_readids_dict (dict): Dictionary mapping of taxon IDs to their corresponding lists of read IDs
        fq1 (str/path): Path to forward fastq file
        fq2 (str/path): Path to reverse fastq file
        outdir (str/path): Path to output directory
    """
    ## load in fastq files as dictionaries for constant-time lookup
    fq1_dict = SeqIO.index(fq1, "fastq")
    fq2_dict = SeqIO.index(fq2, "fastq")

    ## check if read IDs have slash notations
    chosen_ref = list(fq1_dict.keys())[0]
    if chosen_ref[-2:] == "/1":
        slashes = True
    else:
        slashes = False

    def fq_write_wrapper(output_taxid, sample_id, tax_to_readids_dict, fq1_dict, fq2_dict, outdir):
        """Function that wraps around actual file I/O, run in parallel

        Args:
            output_taxid (str): Taxonomic ID to write
            tax_to_readids_dict (dict): Dictionary mapping of taxon IDs to their corresponding lists of read IDs
            fq1_dict (_IndexedSeqFileDict): Dictionary of forward fastq records
            fq2_dict (_IndexedSeqFileDict): Dictionary of reverse fastq records
            outdir (str/path): Path to output directory
        """
        ## initiate counter
        written = 0

        ## initialise files to write to
        R1_file = open(os.path.join(outdir, f"{sample_id}_{output_taxid}_R1.fq"), "w")
        R2_file = open(os.path.join(outdir, f"{sample_id}_{output_taxid}_R2.fq"), "w")

        ## iterate over read ids in list and dump to files
        for read_id in tax_to_readids_dict[output_taxid]:
            if slashes:
                R1_file.write(fq1_dict[read_id+"/1"].format("fastq"))
                R2_file.write(fq2_dict[read_id+"/2"].format("fastq"))
                written += 1
            else:
                R1_file.write(fq1_dict[read_id].format("fastq"))
                R2_file.write(fq2_dict[read_id].format("fastq"))
                written += 1

    ## initialise parallel function calls
    with futures.ThreadPoolExecutor(max_workers=4) as executor:
        ## get list of functions for execution
        functions = [executor.submit(fq_write_wrapper(taxid, sample_id, tax_to_readids_dict, fq1_dict, fq2_dict, outdir)) for taxid in list(tax_to_readids_dict.keys())]
        ## make main wait until functions conclude
        futures.wait(functions)


def sort_reads(sample_id: str, kraken_output: str, mode: str, fastq1: str, fastq2: str, ref_json_file: str, outdir: str, update_output: bool, condense: bool = False, taxon_list: list = None):
    """Control flow of taking args and producing output fastq files

    Args:
        sample_id (str): Sample ID
        kraken_output (str/path): Path to kraken2 output file containing taxid to read_id mappings
        mode (str): Which mode to sort in ["unique", "tree"]
        fastq1 (str/path): Path to forward fastq file
        fastq2 (str/path): Path to reverse fastq file
        condense (bool, optional): If mode == tree, whether to condense reads by species instead of reference. Defaults to False.
        update_output (bool, optional): Whether to update ref_json inplace. Defaults to True.
        taxon_list (str, optional): If mode == unique, list of taxids to extract reads for eg. taxid1,taxid2,taxid3. Defaults to None.
        ref_json_file (str/path, optional): Path to ref_json file produced by kraken2ref. Defaults to None.
    """

    ## time for logging
    NOW = f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S}"

    ## read in kraken output file, populate dict {taxid: [readid1, readid2...]}
    file_handle = open(kraken_output, "r")
    tax_to_read_ids = {}
    read_count = 0
    classified_reads_count = 0
    for line in file_handle.readlines():
        decomp = line.split("\t")
        if decomp[0] == "C":
            if decomp[2] in tax_to_read_ids.keys():
                tax_to_read_ids[decomp[2]].append(decomp[1])
            else:
                tax_to_read_ids[decomp[2]] = [decomp[1]]
            classified_reads_count += 1
        read_count += 1

    logging.debug(f"Found {read_count} read pairs.")
    logging.debug(f"Of which {classified_reads_count} are classified.")

    ############################
    #                          #
    ####    UNIQUE MODE    #####
    #                          #
    ############################

    ## Extract only reads uniquely assigned to specified taxa.
    if mode == "unique":
        if not taxon_list:
            sys.stderr.write("No taxa provided for mode: unique. Exiting.")
            sys.exit(0)
        ## generate list of taxids from CLI arg
        taxids_to_extract = [str(i) for i in taxon_list.split(",")]

        ## generate dict {taxid1: [readid1, readid2...], taxid2: [readid11, readid12...]}
        umode_tax_to_reads = {k: v for k, v in tax_to_read_ids.items() for k in taxids_to_extract}
        ## generate dict {{taxid1: num_reads1, taxid2: num_reads2}}
        umode_numreads_per_taxon = {k: len(v) for k, v in umode_tax_to_reads.items()}

        ## dump to file
        dump_to_file(sample_id, umode_tax_to_reads, fastq1, fastq2, outdir)
        file_read_counts = umode_numreads_per_taxon
        reads_written = []
        ## collect reads that were written for summary logging
        for k, v in umode_tax_to_reads.items():
            reads_written.extend(v)

    ############################
    #                          #
    ####     TREE MODE     #####
    #                          #
    ############################

    ## Extract reads assigned to chosen reference and also all reads in subtree from which ref was chosen
    if mode == "tree":
        ref_json = json.load(open(ref_json_file))
        if len(ref_json["metadata"]["selected"]) == 0:
            sys.stderr.write(f"No FASTQ files to generate for sample: {sample_id}: no reference taxids selected.")
            sys.exit(0)

        ## get selected refs from ref_json
        selected_refs = [str(i) for i  in ref_json["metadata"]["selected"]]

        ## initialise dict {taxid1: [readid1, readid2...], taxid2: [readid11, readid12...]}
        tmode_tax_to_reads = {k: [] for k in selected_refs}

        ## set up condense_mode dict {parent1: [taxid1, taxid2...]} in case needed
        cmode_parent_to_refs = {}

        ## populate dict {taxid1: [readid1, readid2...], taxid2: [readid11, readid12...]}
        for ref in selected_refs:
            data = ref_json["outputs"][ref]
            parent = data["source_taxid"]
            if parent in cmode_parent_to_refs.keys():
                cmode_parent_to_refs[parent].append(ref)
            else:
                cmode_parent_to_refs[parent] = [ref]
            taxa_in_this_ref = data["all_taxa"]
            for taxon in taxa_in_this_ref:
                try:
                    tmode_tax_to_reads[ref].extend(tax_to_read_ids[str(taxon)])
                except KeyError as ke:
                    pass

        ## generate dict {{taxid1: num_reads1, taxid2: num_reads2}}
        tmode_numreads_per_taxon = {k: len(v) for k, v in tmode_tax_to_reads.items()}

        ## if not condense, dump to file now
        if not condense:
            file_read_counts = tmode_numreads_per_taxon
            reads_written = []
            ## collect reads that were written for summary logging
            for k, v in tmode_tax_to_reads.items():
                reads_written.extend(v)
            dump_to_file(sample_id, tmode_tax_to_reads, fastq1, fastq2, outdir)

        #############################
        #                           #
        #####  CONDENSED MODE   #####
        #                           #
        #############################

        ## if condense == True, make one filepair per parent (really useful for flu)
        if condense:

            ## initialise dict {parend_taxid1: [readid1, readid2...], parent_taxid2: [readid11, readid12...]}
            cmode_tax_to_reads = {k: set() for k in cmode_parent_to_refs.keys()}
            for k, v in cmode_parent_to_refs.items():
                for leaf_tax in v:
                    cmode_tax_to_reads[k].update(tmode_tax_to_reads[leaf_tax])

            ## generate dict {{taxid1: num_reads1, taxid2: num_reads2}}
            cmode_numreads_per_taxon = {k: len(v) for k, v in cmode_tax_to_reads.items()}

            file_read_counts = cmode_numreads_per_taxon
            reads_written = []
            ## collect reads that were written for summary logging
            for k, v in cmode_tax_to_reads.items():
                reads_written.extend(v)
            dump_to_file(sample_id, cmode_tax_to_reads, fastq1, fastq2, outdir)

    ## populate summary dict
    summary = {
        "info": {
                "time_updated": str(NOW),
                "sort_mode": mode,
                "total_input_reads": read_count,
                "total_classified_reads": classified_reads_count,
                "unclasified_reads": read_count - classified_reads_count,
                },
        "per_taxon": file_read_counts
    }

    ## log info related to condense
    if not condense:
        summary["condense"] = {"condense_by_root": False}
    else:
        summary["condense"] = {
                                "condense_by_root": True,
                                "condense_info": cmode_parent_to_refs
                            }

    if mode == "unique" and not ref_json_file:
        summary["sample_id"] = sample_id
        with open(f"{sample_id}_sort_reads_summary.json", "w") as sort_reads_json:
            json.dump(summary, sort_reads_json, indent=4)
    else:
        ## update ref_json or write summary to file
        with open(ref_json_file, "r+") as data_json:
            data = json.load(data_json)
            data["metadata"]["summary"] = summary
            if update_output:
                data_json.seek(0)
                json.dump(data, data_json, indent=4)
            else:
                new_json_path = ref_json_file.replace("_decomposed.json", "_updated_decomposed.json")
                with open(new_json_path, "w") as new_json:
                    json.dump(data, new_json, indent=4)
    logging.info(f"Wrote {len(file_read_counts.keys())} file-pairs at path {outdir}.\n\n")

    ## get a list of unwritten reads and dump to file as tsv
    full_kraken_out = pd.read_csv(kraken_output, sep = "\t", header = None)
    unwritten = full_kraken_out[~full_kraken_out[1].isin(reads_written)]
    unwritten.to_csv(os.path.join(outdir, f"{sample_id}_unwritten_reads.txt"), sep = "\t", header = False, index = False)

def sort_reads_by_tax(args):
    """Driver function.
    """

    args.mode = args.mode.lower()

    ## map failing conditions to error message
    error_map = {
        args.mode == "tree" and args.taxon_list: f"Cannot use mode: tree with taxon list...\n",
        args.mode == "unique" and not args.taxon_list: f"No taxon Ids provides for unique mode...\n",
        args.mode != "tree" and args.condense: f"Cannot condense outputs when not using mode: tree...\n",
        args.mode == "unique" and (not args.outdir and not args.ref_json): f"Either provide a JSON produced by kraken2ref or provide a valid outdir...\n"
    }

    ## check failing conditions
    for k, v in error_map.items():
        if k:
            sys.stderr.write(v)
            sys.exit(0)

    ## populate args into main() namespace (bit redundant)
    sample_id = args.sample_id
    kraken_output = args.kraken_out
    mode = args.mode
    if not args.mode:
        logging.debug("No sorting mode provided, defaulting to mode: tree.")
        mode = "tree"
    fastq1 = args.fastq1
    fastq2 = args.fastq2

    condense = args.condense
    update_output = args.update
    taxon_list = args.taxon_list
    ref_json_file = args.ref_json

    ## set up logfile
    if ref_json_file:
        full_path_to_ref_json = os.path.abspath(ref_json_file)
        fixed_outdir = os.path.dirname(full_path_to_ref_json)
        logfile = os.path.join(fixed_outdir, f"{sample_id}_kraken2ref.log")
    else:
        fixed_outdir = os.path.abspath(args.outdir)
        logfile = os.path.join(fixed_outdir, f"{sample_id}_sort_reads.log")
        full_path_to_ref_json = None
    logging.basicConfig(format='%(asctime)s | %(levelname)s | %(module)s - %(funcName)s | %(message)s', level=logging.NOTSET, datefmt="%Y-%m-%d %H:%M:%S", filename = logfile)


    ## run sort_reads
    sort_reads(sample_id=sample_id, kraken_output=kraken_output, mode=mode, fastq1=fastq1, fastq2=fastq2, condense=condense, update_output=update_output, taxon_list=taxon_list, ref_json_file=full_path_to_ref_json, outdir=fixed_outdir)

if __name__ == "__main__":
    sort_reads_by_tax(args = collect_args())
