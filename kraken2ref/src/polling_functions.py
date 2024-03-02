import logging
import scipy.stats as sts

def step_thru(freq_dist):
    """Function that steps forward through a frequency distribution and finds the first inflection

    Args:
        freq_dist (list): List of frequencies to be used

    Returns:
        idxs_to_return: Index locations of retained frequencies to map back to X-axis values and retrieve them
    """
    steps = sorted(freq_dist)
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
            # print(steps[i], steps[j])
            break

    try:
        assert break_point
    except UnboundLocalError as ule:
        break_point = None

    filt_freqs = steps[break_point:]
    idxs_to_return = [freq_dist.index(freq) for freq in filt_freqs]
    return idxs_to_return

def step_thru_back(freq_dist):
    """Function that steps backward through reverse-sorted list of frequencies and finds the last inflection point.

    Args:
        freq_dist (list): List of frequencies to be used

    Returns:
        idxs_to_return: Index locations of retained frequencies to map back to X-axis values and retrieve them
    """
    steps = sorted(freq_dist, reverse=True)
    max_step = 100000000000000000
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
    idxs_to_return = [freq_dist.index(freq) for freq in filt_freqs]
    return idxs_to_return

def poll_leaves(end_nodes, data_dict, parent_selected = False):
    """Apply polling functions to end nodes, treating data as a frequency distribution of hits v/s end nodes

    Args:
        end_nodes (list): List of indexed end nodes to examine
        data_dict (dict): Dictionary representation of kraken2 report

    Returns:
        filtered_end_nodes (list): List of indexed end nodes retained after filtration
        surprise_prefilter (float): Entropy of frequency distribution prior to filtration
        surprise_postfilter (float): Entropy of frequency distribution after filtration
    """

    if len(end_nodes) == 1:
        logging.info("Only one non-spurious leaf-node found. NOT POLLING.")
        return end_nodes, 0, 0

    logging.info("Now polling...")

    filt_data_dict = {k: data_dict[k] for k in data_dict.keys() if k in end_nodes}

    filtered_end_nodes = []
    if parent_selected:
        nodes, freq_dist = [[k for k in filt_data_dict.keys()], [filt_data_dict[k][2] for k in filt_data_dict.keys()]]
    else:
        nodes, freq_dist = [[k for k in filt_data_dict.keys()], [filt_data_dict[k][0] for k in filt_data_dict.keys()]]
    prob_dist = [i/sum(freq_dist) for i in freq_dist]
    logging.debug(f"Frequency Distribution: {freq_dist}")
    logging.debug(f"Probability Distribution: {prob_dist}")
    surprise_prefilter = sts.entropy(prob_dist)
    logging.debug(f"Pre-filter Entropy: {surprise_prefilter}")

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

    max_freq = max(freq_dist)
    max_node = nodes[freq_dist.index(max_freq)]

    if mode == "conservative":
        ## we actually CAN think of this as the left half of a "normal" distribution
        ## or essentially a "sorted" version of a normal distribution
        ## not using this approach for now, but retained for consideration

        # std_dev = np.std(freq_dist)
        # include_range = 2*std_dev
        # for i in range(len(freq_dist)):
        #     if (max_freq-include_range) <= freq_dist[i] <= max_freq:
        #         filtered_end_nodes.append(nodes[i])

        ## using step_thru to step forward in ascending sorted dist
        idxs_to_keep = step_thru(freq_dist)
        for idx in idxs_to_keep:
            filtered_end_nodes.append(nodes[idx])

    if mode == "step":
        ## using step_thru_back to step forward in **descending** sorted dist
        idxs_to_keep = step_thru_back(freq_dist)
        for idx in idxs_to_keep:
            filtered_end_nodes.append(nodes[idx])

    ## caveman method
    if mode == "max":
        filtered_end_nodes.append(max_node)

    logging.debug(f"Selected nodes: {filtered_end_nodes}")
    filt_freqs = [filt_data_dict[k][0] for k in filtered_end_nodes]
    filt_prob_dist = [i/sum(filt_freqs) for i in filt_freqs]
    surprise_postfilter = sts.entropy(filt_prob_dist)
    logging.debug(f"Post-filter Entropy: {surprise_postfilter}")

    return filtered_end_nodes, surprise_prefilter, surprise_postfilter
