## python imports
import os
import sys
import json
import logging
import datetime
import pandas as pd

## import package modules
from kraken2ref.taxonomytree import TaxonomyTree
from kraken2ref.poll import Poll

## retrieve version
try:
    from kraken2ref.version import version as __version__
except ImportError as ie:
    __version__ = "dev/test"

class KrakenProcessor:
    """User-facing driver class to analyse kraken2 taxonomic report
    """
    def __init__(self, sample_id: str):
        """Initialiser

        Args:
            sample_id (str): Sample ID
        """

        ## collect and setup metadata info
        NOW = f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S}"
        self.sample_id = sample_id

        self.metadata = {
                            "k2r_version": __version__,
                            "sample": self.sample_id,
                            "timestamp": NOW
                        }

        logging.info(f"\nkraken2ref version = {__version__}\nSTARTED = {NOW}\nSample: {sample_id}\n")

    def read_report(self, kraken_report_file: str):
        """Function to ingest kraken2 taxonomic report

        Args:
            kraken_report_file (str/path): Path to kraken2 taxonomic report

        Returns:
            kraken_report (pd.DataFrame): Dataframe of kraken2 taxonomic report
            has_minimizer_info (bool): Whether the input report contains minimizer data
        """

        logging.info("Reading report...")

        ## read in kraken report and collect lists of data needed
        self.report_file = kraken_report_file
        kraken_report = pd.read_csv(kraken_report_file, sep = "\t", header = None)
        try:
            kraken_report.columns = ["pct_comp", "cumulative_num_reads", "unique_num_reads", "cumulative_minimizers", "unique_minimizers", "taxon_level", "taxon_id", "desc_name"]
            has_minimizer_info = True
        except Exception as e:
            kraken_report.columns = ["pct_comp", "cumulative_num_reads", "unique_num_reads", "taxon_level", "taxon_id", "desc_name"]
            has_minimizer_info = False

        return kraken_report, has_minimizer_info

    def find_node_lists(self, kraken_report: pd.DataFrame, has_minimizer_info: bool):
        """Collect node-lists from kraken report and populate data dictionary

        Args:
            kraken_report (pd.DataFrame): Dataframe of kraken2 taxonomic report
            has_minimizer_info (bool): Whether the input report contains minimizer data

        Returns:
            all_nodes_list (list(list)): List of nodes (from species level - "S" - onward only)
                Each node is represented as a tuple, where:
                    node[0] = index of that node in the kraken report
                    node[1] = taxon level of that node ("S"/"S1"/etc)
            data_dict (dict): Dictionary representation of kraken report (from species level - "S" - onward only)
                                The contents of the data_dict are:
                                    key = node (node is a tuple as described above)
                                    value = tuple, where
                                        value[0] = cumulative number of reads assigned to that node and its children nodes
                                        value[1] = number of reads *uniquely* assigned to that node
                                        value[2] = taxonomy ID of that node
                                        if has_minimizer_info:
                                            value[3] = cumulative number of minimizers assigned to that node
                                            value[4] = number of unique minimizers assigned to that node
        """

        logging.info("Collecting node lists and data...")

        ## collect all nodes in kraken report
        nodes = list(zip(kraken_report.index, kraken_report["taxon_level"]))

        ## populate data
        if has_minimizer_info:
            data = list(zip(kraken_report["cumulative_num_reads"], kraken_report[ "unique_num_reads"], kraken_report["taxon_id"], kraken_report["cumulative_minimizers"], kraken_report["unique_minimizers"]))
        else:
            data = list(zip(kraken_report["cumulative_num_reads"], kraken_report[ "unique_num_reads"], kraken_report["taxon_id"]))

        ## construct data_dict
        data_dict = dict(zip(nodes, data))

        ## populate node_lists, identifying trees in the report
        all_node_lists = []
        for k in data_dict.keys():
            if k[1] == "S":
                all_node_lists.append([])
            if all_node_lists:
                if "S" in k[1]:
                    all_node_lists[-1].append(k)

        ## handle no data exception
        if len(all_node_lists) == 0:
            logging.critical(msg = f"NoDataFoundError: No Data in Report: File {self.report_file} does not contain any usable data.\n")
            sys.stderr.write(f"NoDataFoundError: No Data in Report: File {self.report_file} does not contain any usable data.\n")
            sys.exit(0)

        return all_node_lists, data_dict

    def analyse_report(self, input_kraken_report_file: str, input_threshold: int, input_method: str, quiet: bool = False):
        """Driver function to read kraken report, collect data from it, and analyse it to pick references

        Args:
            input_kraken_report_file (str): Path to kraken2 taxonomic report
            input_threshold (int): Minimum number of reads required to pass a leaf node as valid
            input_method (str): Polling method to apply ["kmeans", "tiles", "skew"]
            quiet (bool, optional): _description_. Defaults to False.

        Returns:
            tree_meta_out: Dictionary containing output info. This is written to the output JSON file.
        """

        logging.info("Analysing report...")

        ## add threshold value to metadata
        self.metadata["threshold"] = input_threshold

        ## get kraken report
        kraken_report, has_minimizer_info = self.read_report(kraken_report_file = input_kraken_report_file)
        self._log_to_stderr(f"{has_minimizer_info = }\n", quiet)

        ## generate node_lists and data_dict
        all_node_lists, data_dict = self.find_node_lists(kraken_report = kraken_report, has_minimizer_info = has_minimizer_info)

        ## iterate over node_lists
        ### create trees from each node list
        ### if tree complexity == 0: add to list of simple trees
        ### if not, then add subtrees to the list of simple trees

        ## polling functions expect simple trees
        simple_trees = []
        for node_list in all_node_lists:
            ## make tree from node list
            tree = TaxonomyTree(nodes = node_list)

            ## tree complexity is 0, it is ready to be analysed
            if tree.complexity == 0:
                self._log_to_stderr(f"Adding tree rooted at {tree.root} to simple_trees.\n", quiet)
                logging.debug(f"Adding tree rooted at {tree.root} to simple_trees.\n")
                simple_trees.append(tree)

            ## if not, all subtrees in this tree are considered instead
            else:
                num_simple_trees_in_tree = len(tree.subgraphs)
                logging.debug(f"Tree rooted at {tree.root} contains {num_simple_trees_in_tree} simple sub-trees.")
                self._log_to_stderr(f"Tree rooted at {tree.root} contains {num_simple_trees_in_tree} simple sub-trees.\n", quiet)
                simple_trees.extend([TaxonomyTree(tree = subtree) for subtree in tree.subgraphs])

        logging.info("---")

        ## initialise output dict
        self.tree_meta_out = {}

        ## iterate over trees in list simple_trees
        ### poll each tree and collect outputs
        for simple_tree in simple_trees:
            simple_tree_root = simple_tree.root

            ## loging
            logging.debug(f"Now polling simple tree rooted at {simple_tree_root}")
            logging.debug(f"Subterminals are: {simple_tree.subterminal_nodes}")
            self._log_to_stderr(f"Now polling simple tree rooted at {simple_tree_root}\n", quiet)
            self._log_to_stderr(f"Subterminals are: {simple_tree.subterminal_nodes}\n", quiet)

            ## create a poll object
            poll = Poll(taxonomy_tree = simple_tree, data_dict = data_dict, threshold = input_threshold)

            ## logging
            logging.debug(f"{poll.singleton = }")
            logging.debug(f"Number of leaves = {len(poll.leaves)}")
            self._log_to_stderr(f"{poll.singleton = }\n", quiet)
            self._log_to_stderr(f"Number of leaves = {len(poll.leaves)}\n", quiet)

            ## if only one valid leaf node in tree, update output dict with it
            if poll.singleton:
                self.tree_meta_out.update(self._update_tree_meta(filt_leaves = poll.valid_leaves, simple_source_tree = simple_tree, data_dict = data_dict, parent_selected = poll.parent_selected))

                ## logging
                logging.debug(f"Only one valid node found. Not polling.\n")
                self._log_to_stderr(f"Only one valid node found. Not polling.\n\n", quiet)

                continue

            ## if no valid leaf nodes found, check if parent is valid
            if poll.parent_selected:
                ## logging
                logging.debug(f"No leaf nodes passed the minimum read threshold: {input_threshold}.")
                self._log_to_stderr(f"No leaf nodes passed the minimum read threshold: {input_threshold}.\n", quiet)

                ## if parent is valid, add info to output dict
                if poll.valid_parent:
                    self.tree_meta_out.update(self._update_tree_meta(filt_leaves = poll.valid_subterminals, simple_source_tree = simple_tree, data_dict = data_dict, parent_selected = poll.parent_selected))

                    ## logging
                    logging.debug(f"Jumped up 1 level, selected parent node: {poll.valid_subterminals}.\n")
                    self._log_to_stderr(f"Jumped up 1 level, selected parent node: {poll.valid_subterminals}.\n\n", quiet)

                    continue

                ## if no valid leaves and parent not valid, skip
                else:
                    ## logging
                    logging.debug(f"Jumped up 1 level, parent node did not pass threshold. No data passed on from this graph, rooted at {simple_tree_root}.\n")
                    self._log_to_stderr(f"Jumped up 1 level, parent node did not pass threshold. No data passed on from this graph, rooted at {simple_tree_root}.\n\n", quiet)

                    continue

            ## if >1 valid leaf nodes found, run polling
            else:
                poll.poll_leaves(method = input_method)

                ## logging
                logging.debug(f"Valid leaves = {poll.valid_leaves}\n")
                self._log_to_stderr(f"Valid leaves = {poll.valid_leaves}\n", quiet)

                ## if polling returned no leaf nodes, handle and skip
                if len(poll.filt_leaves) == 0:
                    logging.debug(f"No leaf nodes passed polling")
                    self._log_to_stderr(f"No leaf nodes passed polling.\n\n")
                    continue

                ## update output dict sith polling results
                self.tree_meta_out.update(self._update_tree_meta(filt_leaves = poll.filt_leaves, simple_source_tree = simple_tree, data_dict = data_dict, parent_selected = poll.parent_selected))

                ## logging
                logging.debug(f"Valid leaves = {poll.valid_leaves}\n")
                self._log_to_stderr(f"{poll.filt_leaves = }\n\n", quiet)

        return self.tree_meta_out

    def _update_tree_meta(self, filt_leaves: list, simple_source_tree: TaxonomyTree, data_dict: dict, parent_selected: bool):
        """Function to encapsulate output dict updtes

        Args:
            filt_leaves (list): List of leaf nodes that passed the polling
            simple_source_tree (TaxonomyTree): Tree for which updates are being written
            data_dict (dict): data_dict as described above
            parent_selected (bool): Whether the filt_leaf provided is a parent (non-terminal) node

        Returns:
            tree_meta_chunk: The chunk of the output dict generated for the input leaf nodes
        """
        for filt_leaf in filt_leaves:
            simple_tree_root = simple_source_tree.root ## root node
            simple_tree_root_taxid = data_dict[simple_tree_root][2] ## root taxID
            meta_key = data_dict[filt_leaf][2] ## taxID
            simple_tree_idx = simple_tree_root[0] ## tree_idx
            all_taxa_in_simple_tree = [data_dict[node][2] for node in simple_source_tree.nodes] ## all taxIDs in this simple (sub)tree
            path_to_filt_leaf = simple_source_tree.find_all_paths(graph = simple_source_tree.graph, source = simple_tree_root, target = filt_leaf)[0] ## path as nodes
            path_as_taxids = [data_dict[i][2] for i in path_to_filt_leaf] ## path as taxIDs
            tree_meta_chunk =  {
                                meta_key :
                                            {
                                                "graph_idx": simple_tree_idx,
                                                "source": simple_tree_root,
                                                "source_taxid": simple_tree_root_taxid,
                                                "target": filt_leaf,
                                                "parent_selected": parent_selected,
                                                "all_taxa": all_taxa_in_simple_tree,
                                                "path": path_to_filt_leaf,
                                                "path_as_taxids": path_as_taxids,
                                            }
                            }

        return tree_meta_chunk

    def write_output(self, prefix:str = os.getcwd(), suffix:str = "decomposed"):
        """Function to write output to file. Separated and not run in fn analyse report to allow programmatic use of this module

        Args:
            prefix (str/path, optional): Path to output directory. Defaults to os.getcwd().
            suffix (str/path, optional): Suffix to append to sample name in output filename. Defaults to "decomposed".
        """
        logging.info("Writing output JSON...")

        ## clean up suffix
        if suffix[0] not in [".", "_"]:
            suffix = "_"+suffix
        if ".json" in suffix:
            suffix.replace(".json", "")

        ## explicitly set up filenames
        outfile_name = self.sample_id + suffix + ".json"
        outfile_path = os.path.join(prefix, outfile_name)

        ## collect and organise data to be output
        selected_refs = list(self.tree_meta_out.keys())
        if len(selected_refs) == 0:
            logging.warning(f"No suitable references found in sample: {self.sample_id}.")
            sys.exit(0)

        self.metadata["selected"] = selected_refs
        to_file = {
                    "metadata": self.metadata,
                    "outputs": self.tree_meta_out
                    }

        with open(outfile_path, "w") as outfile:
            json.dump(to_file, outfile, indent=4)
        logging.info(msg = f"Output written to {self.sample_id}_decomposed.json at location {prefix}\n\n")

    def _log_to_stderr(self, message: str, quiet: bool):
        """Helper function to write to stderr based on quiet mode

        Args:
            message (str): The message to write to stderr
            quiet (bool): Whether quiet mode is active
        """
        if quiet:
            return
        sys.stderr.write(message)

