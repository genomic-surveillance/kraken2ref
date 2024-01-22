import re
import pandas as pd
import os.path
from cached_property import cached_property
import logging

from kraken2ref.src.graph_functions import build_graph, get_end_points

logging.basicConfig( format='%(asctime)s %(message)s', level=logging.DEBUG )

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

    def __init__( self, in_file: str, min_abs_reads: int = 5 ):
        self.in_file = in_file
        self.threshold = min_abs_reads

        if not os.path.isfile( self.in_file ):
            raise ValueError(f'path { in_file} does not exist or is not a file')


    def read_kraken_report(self, kraken_report):
        kraken_report = pd.read_csv(kraken_report, sep = "\t", header = None)

        keys = list(zip(kraken_report.index, kraken_report[3]))
        vals = list(zip(list(kraken_report[2]), list(kraken_report[4])))

        data_dict = dict(zip(keys, vals))
        all_node_lists = []
        for k in data_dict.keys():
            if k[1] == "S":
                all_node_lists.append([])
            if all_node_lists:
                if "S" in k[1]:
                    all_node_lists[-1].append(k)

        return all_node_lists, data_dict

    @cached_property
    def pick_reference_taxid(self):
        graphs = []
        all_node_lists, data_dict = self.read_kraken_report(self.in_file)
        for node_list in all_node_lists:
            graphs.append(build_graph(node_list))

        all_paths_lists = []
        all_taxids_lists = []
        for graph in graphs:
            this_graph_paths, this_graph_tax_lists = get_end_points(graph, data_dict, self.threshold)
            # print(this_graph_paths)
            print("\n\n\n")
            for path in this_graph_paths:
                all_paths_lists.append(path)
            for tax_id_list in this_graph_tax_lists:
                all_taxids_lists.append(tax_id_list)

        dump_dict = {all_taxids_lists[i][-1]: [all_taxids_lists[i], all_paths_lists[i]] for i in range(len(all_taxids_lists))}
        # print("\n\n\n\nDUMP DICT = ")
        # for k, [t, p] in dump_dict.items():
        #     print(f"Filename = {k}\ttaxa = {t}\tpath = {p}")

        # for i in range(len(all_paths_lists)):
        #     print(f"{i}. path = {all_paths_lists[i]}\ttax_ids_list = {all_taxids_lists[i]}")

        return dump_dict

    @cached_property
    def taxonomy_graph( self ):
        """
        A graph representation of the taxonomy with number of reads that are assigned directly to each node
        """
        # TODO
        pass