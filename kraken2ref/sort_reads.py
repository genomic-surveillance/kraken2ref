import os, sys
import json
import pandas as pd
from Bio import SeqIO
import datetime
import logging

def sort_reads(sample_id: str, kraken_output: str, mode: str,
        ref_json_file: str, outdir: str, update_output: bool,
        condense: bool = False, taxon_list: list = None,):
    """
    Control flow of taking args and producing output fastq files

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

    def load_reads_written(tax_to_reads):
        """
        Collects and returns a list of reads written for summary logging.

        Parameters:
            tax_to_reads (dict): A dictionary where keys are taxonomy IDs
                         and values are lists of reads that were written at those levels

        Returns:
            reads_written (list): A list of all reads that were written, regardless of taxonomy level

        Examples:
            >>> load_reads_written({'tax1': ['read1', 'read2'], 'tax2': ['read3']})
            ['read1', 'read2', 'read3']
        """
        reads_written = []
        ## collect reads that were written for summary logging
        for k, v in tax_to_reads.items():
            reads_written.extend(v)
        return reads_written

    def write_out_json(sample_id: str, outdir: str,tax_to_reads: dict):
        """
        Writes the taxonomy-to-reads dictionary to a JSON file.

        Parameters:
            sample_id (str): The ID of the sample being written.
            outdir (str): The directory where the output files should be written.
            tax_to_reads (dict): A dictionary where keys are taxonomy IDs and values are lists of reads that were written at those levels.

        Returns:
            None

        Examples:
            >>> write_out_json('sample1', '/path/to/output/directory', {'tax1': ['read1', 'read2'], 'tax2': ['read3']})
            > output file written to /path/to/output/directory/sample1_tax_to_reads.json
        """
        # write tax_to_reads json file
        json_out_path = f"{outdir}/{sample_id}_tax_to_reads.json"
        tax_json_out = open(json_out_path, "w")
        json_content_str = json.dumps(tax_to_reads, indent=4)
        tax_json_out.write(json_content_str)
        sys.stdout.write(f"> output file written to {json_out_path}")

    def compute_numreads_per_taxon(tax_to_reads):
        """
        Computes and returns a dictionary where keys are taxonomy IDs and values are the number of reads associated with each ID.

        Parameters:
            tax_to_reads (dict): A dictionary where keys are taxonomy IDs and values are lists of reads that were written at those levels.

        Returns:
            dict: A dictionary where keys are taxonomy IDs and values are the number of reads associated with each ID.

        Examples:
            >>> compute_numreads_per_taxon({'tax1': ['read1', 'read2'], 'tax2': ['read3']})
            {'tax1': 2, 'tax2': 1}
        """
        return {k: len(v) for k, v in tax_to_reads.items()}
    
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
        tax_to_reads = {k: v for k, v in tax_to_read_ids.items() for k in taxids_to_extract}
        ## generate dict {{taxid1: num_reads1, taxid2: num_reads2}}
        numreads_per_taxon = compute_numreads_per_taxon(tax_to_reads)

        
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
        numreads_per_taxon = compute_numreads_per_taxon(tmode_tax_to_reads)

        ## if not condense, dump to file now
        if not condense:
            #file_read_counts = tmode_numreads_per_taxon
            ## collect reads that were written for summary logging
            reads_written = load_reads_written(tmode_tax_to_reads)
            write_out_json(sample_id, outdir,tmode_tax_to_reads)

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
                
            # convert set to list, json is not happy to dump files containing sets
            for k, v in cmode_tax_to_reads.items():
                cmode_tax_to_reads[k] = list(v)

            ## generate dict {{taxid1: num_reads1, taxid2: num_reads2}}
            numreads_per_taxon = compute_numreads_per_taxon(cmode_tax_to_reads)

            reads_written = load_reads_written(cmode_tax_to_reads)
            write_out_json(sample_id, outdir, cmode_tax_to_reads)

    ## populate summary dict
    summary = {
        "info": {
                "time_updated": str(NOW),
                "sort_mode": mode,
                "total_input_reads": read_count,
                "total_classified_reads": classified_reads_count,
                "unclasified_reads": read_count - classified_reads_count,
                },
        "per_taxon": numreads_per_taxon
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

    logging.info(f"Wrote {len(numreads_per_taxon.keys())} file-pairs at path {outdir}.\n\n")

    ## get a list of unwritten reads and dump to file as tsv
    full_kraken_out = pd.read_csv(kraken_output, sep = "\t", header = None)
    unwritten = full_kraken_out[~full_kraken_out[1].isin(reads_written)]
    unwritten.to_csv(os.path.join(outdir, f"{sample_id}_unwritten_reads.txt"), sep = "\t", header = False, index = False)


# instantiating the decorator
#@profile
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

    condense = args.condense
    update_output = args.update
    taxon_list = args.taxon_list
    ref_json_file = args.ref_json

    absolute_outdir = os.path.abspath(args.outdir)

    ## set up logfile
    if ref_json_file:
        full_path_to_ref_json = os.path.abspath(ref_json_file)
        log_outdir = os.path.dirname(full_path_to_ref_json)
        logfile = os.path.join(log_outdir, f"{sample_id}_kraken2ref.log")
    else:
        logfile = os.path.join(absolute_outdir, f"{sample_id}_sort_reads.log")
        full_path_to_ref_json = None
    logging.basicConfig(format='%(asctime)s | %(levelname)s | %(module)s - %(funcName)s | %(message)s', level=logging.NOTSET, datefmt="%Y-%m-%d %H:%M:%S", filename = logfile)


    ## run sort_reads
    sort_reads(
        sample_id=sample_id,
        kraken_output=kraken_output,
        mode=mode,
        condense=condense,
        update_output=update_output,
        taxon_list=taxon_list,
        ref_json_file=full_path_to_ref_json,
        outdir=absolute_outdir)

