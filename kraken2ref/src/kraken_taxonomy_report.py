import pandas as pd
import os, sys
# from cached_property import cached_property
import logging
import json
import datetime

from kraken2ref.src.graph_functions import build_graph, find_valid_graphs
from kraken2ref.src.sort_reads_by_ref import write_fastq

try:
    from kraken2ref.version import version as __version__
except ImportError as ie:
    __version__ = "dev/test"

class KrakenTaxonomyReport():
    """
    This class handles a kraken2 taxonomy report output file.

    Parameters:
        in_file: str, required
            path to the kraken2 taxonomy report file (should be called "report.txt")

        min_abs_reads: int, optional, defaults to 5
            minimum absolute number of reads that are directly assigned to a taxon. A taxon with
            fewer reads will be discarded.

    """

    def __init__(self, sample_id: str):

        NOW = f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S}"
        self.sample_id = sample_id

        self.metadata = {
                            "k2r_version": __version__,
                            "sample": self.sample_id,
                            "timestamp": NOW
                        }

        logging.info(f"\nkraken2ref version = {__version__}\nSTARTED = {NOW}\nSample: {sample_id}")


    def read_kraken_report(self, kraken_report_file:str):
        """Read in kraken2 report and produce outputs used to build and find paths in taxonomy graph.

        Args:
            kraken_report (str/pathlike): Path to input kraken2 taxonomy report.

        Returns:
            all_nodes_list (list): List of nodes (from species level - "S" - onward only)
                Each node is represented as a tuple, where:
                    node[0] = index of that node in the kraken report
                    node[1] = taxon level of that node ("S"/"S1"/etc)
            data_dict (dict): Dictionary representation of kraken report (from species level - "S" - onward only)
                The contents of the data_dict are:
                    key = node (node is a tuple as described above)
                    value = tuple, where
                        value[0] = number of reads assigned to that node
                        value[1] = taxonomy ID of that node
        """
        ## read in kraken report and collect lists of data needed
        kraken_report = pd.read_csv(kraken_report_file, sep = "\t", header = None)

        num_hits = list(kraken_report[2])
        num__cumulative_hits = list(kraken_report[1])
        tax_levels = list(kraken_report[3])
        tax_ids = list(kraken_report[4])

        keys = list(zip(kraken_report.index, tax_levels))
        vals = list(zip(num_hits, tax_ids, num__cumulative_hits))

        ## construct data dict
        data_dict = dict(zip(keys, vals))

        ## iterate over keys in data dict and collect nodes for each graph in report
        all_node_lists = []
        for k in data_dict.keys():
            if k[1] == "S":
                all_node_lists.append([])
            if all_node_lists:
                if "S" in k[1]:
                    all_node_lists[-1].append(k)

        if len(all_node_lists) == 0:
            logging.critical(msg = f"NoDataFoundError: No Data in Report: File {kraken_report_file} does not contain any usable data.\n")
            sys.stderr.write(f"NoDataFoundError: No Data in Report: File {kraken_report_file} does not contain any usable data.\n")
            sys.exit(0)

        return all_node_lists, data_dict

    def pick_reference_taxid(self, in_file: str, outdir: str, min_abs_reads: int = 5):
        """Build all graphs contained in the kraken2 taxonoic report;
            From each graph, identify nodes that are assigned more reads
                than the threshold (passing nodes);
            Find all paths leading to passing nodes;
            Collect all paths found in the kraken2 report;
            Collect the taxon IDs for each node in each path

        Returns:
            dump_dict (dict): A dictionary with data that can be dumped to file.
                The contents of dump_dict are:
                    key = tax ID chosen for each path
                    value = list where:
                        value[0] = list of tax_ids for the path leading to this key
                        value[1] = list of node (ie the path) leading to this key
        """
        logging.info(msg = f"CMD: kraken2r -s {self.sample_id} parse_report \n\t\t-i {in_file} \n\t\t-o {outdir} \n\t\t-t {min_abs_reads} \n\t\t")

        if not os.path.isfile(in_file):
            logging.critical(msg = f"FileNotFoundError: Missing Input: Path {in_file} does not exist or is not a file.\n")
            raise FileNotFoundError(f"Missing Input: Path {in_file} does not exist or is not a file.\n")


        self.in_file = in_file
        self.threshold = min_abs_reads
        self.outdir = outdir
        self.graphs = []
        self.metadata["threshold"] = self.threshold
        self.all_node_lists, self.data_dict = self.read_kraken_report(self.in_file)
        if len(self.all_node_lists) == 0:
            logging.warning(msg=f"The given report: {self.in_file} does not contain any graphs. Sample ID: {self.sample_id} will have no outputs. Exiting...")
            sys.stderr.write("The given report: {self.in_file} does not contain any graphs. Sample ID: {self.sample_id} will have no outputs. Exiting...")
            sys.exit(0)
        for node_list in self.all_node_lists:
            self.graphs.append(build_graph(node_list))

        graph_meta_dict = {}
        for graph in self.graphs:
            graph_meta_dict.update(find_valid_graphs(graph, self.data_dict, self.threshold))

        self.graph_meta = graph_meta_dict

        selected_refs = list(graph_meta_dict.keys())
        if len(selected_refs) == 0:
            logging.warning(msg=f"There are no selected references in sample {self.sample_id}.")
        self.metadata["selected"] = selected_refs

        to_json = {
            "metadata": self.metadata,
            "outputs": self.graph_meta
        }

        with open(os.path.join(self.outdir, self.sample_id+"_decomposed.json"), "w") as outfile:
            json.dump(to_json, outfile, indent=4)
        logging.info(msg = f"Output written to {self.sample_id}_decomposed.json at location {outdir}\n\n")

    def sort_reads_by_ref(self, sample_id: str, fq1: str, fq2:str, kraken_out:str, update_output:bool = True, ref_data = None, condense = False):
        logging.info(msg = f"CMD: kraken2r -s {self.sample_id} ref_sort_reads \n\t\t-fq1 {fq1} \n\t\t-fq2 {fq2} \n\t\t-k {kraken_out} \n\t\t-r {ref_data} \n\t\t-u {update_output} \n\t\t-c {condense}\n\n")
        self.summary = write_fastq(sample_id, fq1, fq2, kraken_out, update_output, ref_data, condense)

