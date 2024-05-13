import logging, sys
import numpy as np
import scipy.stats as sts
from sklearn.cluster import KMeans

from kraken2ref.taxonomytree import TaxonomyTree

class Poll:
    """Class to encapsulate polling functions and handle related exceptions/standard behaviours
    """
    def __init__(self, taxonomy_tree: TaxonomyTree, data_dict:dict, threshold: int):
        """Initialiser

        Args:
            taxonomy_tree (TaxonomyTree): Input TaxonomyTree to analyse
            data_dict (dict): data_dict containing data for the input tree
            threshold (int): Minimum read threshold to use as cutoff for leaf node validity
        """
        ## set up attributes
        self.leaves = taxonomy_tree.leaf_nodes
        self.class_method = []

        self.nodes_to_poll = None
        self.dist = None
        self.prob_dist = None
        self.parent_selected = False
        self.valid_parent = False
        self.singleton = False

        ## find the number of valid leaf nodes
        self.valid_leaves = [i for i in self.leaves if data_dict[i][1] > threshold]

        ## if no valid leaf nodes, jump up one level
        if len(self.valid_leaves) == 0:
            self.parent_selected = True
            self.valid_subterminals = [i for i in taxonomy_tree.subterminal_nodes if data_dict[i][0] > threshold]
            ## if subterminal node is valid, return this info
            if len(self.valid_subterminals) > 0 and len(self.leaves) > 1:
                self.valid_parent = True
                valid_leaves_dict = {i: data_dict[i] for i in self.leaves}
                max_freq = max(list(valid_leaves_dict.values()))
                self.max_leaves = [k for k, v in valid_leaves_dict.items() if v == max_freq]

        ## if only one valid leaf node, mark singleton attribute as true: no polling needed
        if len(self.valid_leaves) == 1:
            self.singleton = True

        ## if >1 valid leaf node, set up polling
        if len(self.valid_leaves) > 1:
            filt_data_dict = {k: data_dict[k] for k in data_dict.keys() if k in self.valid_leaves}
            nodes, freq_dist = [[k for k in filt_data_dict.keys()], [filt_data_dict[k][1] for k in filt_data_dict.keys()]]

            self.nodes_to_poll = nodes
            self.dist = freq_dist
            self.prob_dist = [i/sum(self.dist) for i in self.dist]

    def step_thru(self):
        """Function that steps forward through a frequency distribution and finds the first inflection

        Returns:
            idxs_to_return: Index locations of retained frequencies to map back to X-axis values and retrieve them
        """
        steps = sorted(self.dist)
        max_step = 0
        for i in range(len(steps) - 1):
            j = i+1
            step = steps[j] - steps[i]
            # print(steps[i], steps[j], step)
            if step > max_step:
                max_step = step
                i += 1
            elif step == 0 or step == max_step:
                pass
            else:
                break_point = i
                break

        try:
            assert break_point
        except UnboundLocalError as ule:
            break_point = None

        filt_freqs = steps[break_point:]
        idxs_to_return = [self.dist.index(freq) for freq in filt_freqs]
        return idxs_to_return

    def step_thru_back(self):
        """Function that steps backward through reverse-sorted list of frequencies and finds the last inflection point.

        Returns:
            idxs_to_return: Index locations of retained frequencies to map back to X-axis values and retrieve them
        """
        steps = sorted(self.dist, reverse=True)
        # max_step = 100000000000000000
        max_step = 65_535
        for i in range(len(steps) - 1):
            j = i+1
            step = steps[i] - steps[j]
            if step > max_step:
                break_point = j
                break
            elif step == 0:
                pass
            else:
                max_step = step
                i += 1

        try:
            assert break_point
        except UnboundLocalError as ule:
            break_point = None

        filt_freqs = steps[:break_point]
        idxs_to_return = [self.dist.index(freq) for freq in filt_freqs]
        return idxs_to_return

    def poll_with_skew(self):
        """Apply polling functions to end nodes, treating data as a frequency distribution of hits v/s end nodes

        Returns:
            filt_leaves (list): List of indexed end nodes retained after filtration
            surprise_postfilter (float): Entropy of frequency distribution after filtration
        """

        freq_dist = self.dist
        prob_dist = self.prob_dist
        logging.debug(f"Frequency Distribution: {freq_dist}")
        logging.debug(f"Probability Distribution: {prob_dist}")

        if len(prob_dist) < 8:
            padding = [0]*(8-len(prob_dist))
            prob_dist.extend(padding)

        skew_test = sts.skewtest(prob_dist)
        if skew_test.pvalue < 0.005:
            logging.info("Mode = MAX")
            mode = "max"
        if 0.005 < skew_test.pvalue < 0.05:
            logging.info("Mode = STEP")
            mode = "step"
        if 0.05 < skew_test.pvalue:
            logging.info("Mode = CONSERVATIVE")
            mode = "conservative"
        self.class_method.append(mode)
        filt_leaves = []
        if mode == "conservative":
            ## using step_thru to step forward in ascending sorted dist
            idxs_to_keep = self.step_thru()
            for idx in idxs_to_keep:
                filt_leaves.append(self.nodes_to_poll[idx])

        if mode == "step":
            ## using step_thru_back to step forward in **descending** sorted dist
            idxs_to_keep = self.step_thru_back()
            for idx in idxs_to_keep:
                filt_leaves.append(self.nodes_to_poll[idx])

        ## caveman method
        if mode == "max":
            max_freq = max(freq_dist)
            max_node = self.nodes_to_poll[freq_dist.index(max_freq)]
            filt_leaves.append(max_node)

        logging.debug(f"Selected nodes: {filt_leaves}")
        filt_prob_dist = []
        for filt_leaf in filt_leaves:
            idx = self.nodes_to_poll.index(filt_leaf)
            filt_prob_dist.append(self.prob_dist[idx])

        surprise_postfilter = sts.entropy(filt_prob_dist)
        logging.debug(f"Post-filter Entropy: {surprise_postfilter}")

        return filt_leaves, surprise_postfilter

    def poll_with_tiles(self):
        """Apply polling using quantile-based outlier analysis

        Returns:
            filt_leaves (list): List of indexed end nodes retained after filtration
            surprise_postfilter (float): Entropy of frequency distribution after filtration
        """

        freq_dist = self.dist
        q1 = np.quantile(freq_dist, 0.25)
        q3 = np.quantile(freq_dist, 0.75)
        iqr = q3 - q1
        left_fence = q1 - (1.5 * iqr)
        right_fence = q3 + (1.5 * iqr)

        logging.debug(f"Lower quartile fence = {left_fence}")
        logging.debug(f"Lower quartile fence = {right_fence}")

        filt_leaves = []
        filt_prob_dist = []
        for i, freq in enumerate(freq_dist):
            if freq > right_fence:
                filt_leaves.append(self.nodes_to_poll[i])
                filt_prob_dist.append(self.prob_dist[i])

        logging.debug(f"Selected nodes: {filt_leaves}")
        surprise_postfilter = sts.entropy(filt_prob_dist)
        logging.debug(f"Post-filter Entropy: {surprise_postfilter}")

        return filt_leaves, surprise_postfilter

    def poll_with_kmeans(self):
        """Apply polling using kmeans-based outlier analysis

        Returns:
            filt_leaves (list): List of indexed end nodes retained after filtration
            surprise_postfilter (float): Entropy of frequency distribution after filtration

        """

        freq_arr = np.array(self.dist).reshape(-1, 1)
        kmeans = KMeans(n_clusters = 1)
        kmeans.fit(freq_arr)
        distances = kmeans.transform(freq_arr)
        sorted_idx = np.argsort(distances.ravel())[::-1][:5]
        med = np.median(freq_arr)
        filt_leaves = []
        filt_prob_dist = []
        for idx in sorted_idx:
            if self.dist[idx] > med:
                filt_leaves.append(self.nodes_to_poll[idx])
                filt_prob_dist.append(self.prob_dist[idx])

        surprise_postfilter = sts.entropy(filt_prob_dist)
        logging.debug(f"Post-filter Entropy: {surprise_postfilter}")

        return filt_leaves, surprise_postfilter

    def poll_leaves(self, method: str = "kmeans"):
        """Driver function to run polling

        Args:
            method (str): Polling method to use. Defaults to "kmeans"
        """

        ## if singleton: return
        if self.singleton:
            logging.info("Only one non-spurious leaf-node found. NOT POLLING.")
            self.class_method.append("singleton")
            self.filt_leaves = self.valid_leaves
            self.surprise_prefilter = 0
            self.postfilter_surprise = 0
            return

        ## if no valid leaves, and parent is valid, return parent
        ## if parent not valid, return None
        if self.parent_selected:
            self.class_method.append("singleton")
            if self.valid_parent:
                self.filt_leaves = self.valid_subterminals
                self.surprise_prefilter = 0
                self.postfilter_surprise = 0
            else:
                self.filt_leaves = None
                self.surprise_prefilter = 0
                self.postfilter_surprise = 0
                return

        ## ensure method is valid
        method = method.lower()
        self.class_method.append(method)
        try:
            assert method in ["skew", "tiles", "kmeans"]
        except AssertionError as ae:
            logging.warning(f"Invalid polling method specified. Defaulting to 'kmeans'")
            method = "kmeans"

        logging.debug(f"Frequency Distribution: {self.dist}")
        logging.debug(f"Probability Distribution: {self.prob_dist}")
        self.surprise_prefilter = sts.entropy(self.prob_dist)
        logging.debug(f"Pre-filter Entropy: {self.surprise_prefilter}")

        if method == "skew":
            self.filt_leaves, self.postfilter_surprise = self.poll_with_skew()

        if method == "tiles":
            self.filt_leaves, self.postfilter_surprise = self.poll_with_tiles()

        if method == "kmeans":
            self.filt_leaves, self.postfilter_surprise = self.poll_with_kmeans()

