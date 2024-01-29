import re
import pandas as pd
import os
from cached_property import cached_property
import logging

from kraken2ref.src.graph_functions import build_graph, get_graph_endpoints

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)

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

    def __init__(self, in_file: str, min_abs_reads: int = 5):
        self.in_file = in_file
        self.threshold = min_abs_reads

        if not os.path.isfile( self.in_file ):
            raise ValueError(f'path { in_file} does not exist or is not a file')


    def read_kraken_report(self, kraken_report):
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
        kraken_report = pd.read_csv(kraken_report, sep = "\t", header = None)

        num_hits = list(kraken_report[2])
        tax_levels = list(kraken_report[3])
        tax_ids = list(kraken_report[4])

        keys = list(zip(kraken_report.index, tax_levels))
        vals = list(zip(num_hits, tax_ids))

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

        return all_node_lists, data_dict

    def pick_reference_taxid(self):
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
        self.graphs = []
        self.all_node_lists, self.data_dict = self.read_kraken_report(self.in_file)
        for node_list in self.all_node_lists:
            self.graphs.append(build_graph(node_list))

        graph_meta_dict = get_graph_endpoints(graphs=self.graphs, data_dict=self.data_dict, threshold=self.threshold)
        self.graph_meta = graph_meta_dict

        return graph_meta_dict

