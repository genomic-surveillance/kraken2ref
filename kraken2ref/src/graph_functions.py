import logging
from kraken2ref.src.taxonlevel import TaxonLevel, find_parent
from kraken2ref.src.polling_functions import poll_leaves

def get_graph_complexity(graph):
    """Get a measure of graph "complexity" based on multiplicity of nodes at each level.
        Complexity is:
            0 IF graph has multiplicity at only the leaf level
            1 IF graph has multiplicity at multiple levels AND ALL PATHS END AT TAXON LEVEL
            2 IF graph has multiplicity at multiple levels AND THERE ARE MULTIPLE TERMINAL TAXON LEVELS


    Args:
        graph (dict): A dictionary representation of the graph, where:
                    key = Source node
                    value = Target node
                    (Nodes are represented as described elsewhere)

    Returns:
        complexity (int): Integer description of graph complexity.
    """
    terminal_levels = set()
    terminal_nodes = []
    for k, v in graph.items():
        if len(v) == 0:
            terminal_nodes.append(k)
            terminal_levels.add(k[1])
    if len(terminal_levels) > 1:
        complexity = 2
    else:
        breadth_indicator_level = (TaxonLevel(list(terminal_levels)[0]) - 1).lvl
        indicator_nodes = [node for node in graph.keys() if node[1] == breadth_indicator_level]
        if len(indicator_nodes) > 1:
            complexity = 1
        else:
            complexity = 0

    return complexity

def build_graph(indexed_nodes):
    """Function to build a graph representation from a list of indexed nodes.

    Args:
        indexed_nodes (list(tuples)): List of nodes; each node is represented as a tuple,
                where:
                    node[0] = index of that node in the original data
                    node[1] = taxon level of that node ("S"/"S1"/etc)

    Returns:
        graph (dict): A dictionary representation of the graph, where:
                    key = Source node
                    value = Target node
                    (Nodes are represented as described above)
    """

    ## initialise required data structures for graph building
    graph = {}

    plain_levels = []  ## e.g. ["S", "S1", "S2", "S3", "S3"]

    ## levels_as_taxa is used to backtrack
    ## good to have it computed once as will probably be needed often
    levels_as_taxa = []  ## same as above except as list of TaxonLevel objects

    ## iterate over indexed_nodes list
    ## (looks like: [(11,"S"), (12,"S1)", (13,"S2"), (14,"S3"), (15,"S3")] )
    ## to populate data structures
    for indexed_node in indexed_nodes:
        graph[indexed_node] = []
        plain_levels.append(indexed_node[1])
        levels_as_taxa.append(TaxonLevel(indexed_node[1]))

    ## set up counters for list trawling
    start = 0
    # offset = indexed_nodes[0][0]

    ## trawl over list to find edges in graph
    while start < len(plain_levels)-1:
        end = start + 1
        left_node = TaxonLevel(plain_levels[start])
        right_node = TaxonLevel(plain_levels[end])

        ## if left node less than right node: that's an edge
        if left_node < right_node:
            ## generalise this, do not depend on offset
            graph[indexed_nodes[start]].append(indexed_nodes[end])
            # graph[(start+offset,plain_levels[start])].append((end+offset, plain_levels[end]))

        ## if left node equals right node: oops! find the closest parent
        if left_node >= right_node:
            sublist = levels_as_taxa[:end]
            parent = right_node - 1
            parent_idx = find_parent(sublist, parent)
            ## generalise this, do not depend on offset
            graph[indexed_nodes[parent_idx]].append(indexed_nodes[end])
            # graph[(parent_idx+offset, parent.lvl)].append((end+offset, plain_levels[end]))
        start += 1

    return graph

def find_all_paths(graph, source, target, path=[]):
    """Generic function to find all paths from given source node to given target node.

    Args:
        graph (dict): A dictionary representation of the graph, where:
                        key = Source node
                        value = Target node
        source (hashable object): Source node from which to start
        target (hashable object): Target node to which paths should lead
        path (list, optional): A known path from given source to given target [Default = []]

    Returns:
        list(lists): A list of paths (where each path is a list)
    """
    ## initialise path by adding source
    path = path + [source]

    ## handle edge cases
    if source == target:
        return [path]
    if source not in graph:
        return []

    ## paths is a list in case of multiple paths from source to target
    paths = []

    ## iterate over graph dict (i.e. edges)
    for node in graph[source]:
        if node not in path:
            ## recursively find all "sub-paths" and update paths list
            newpaths = find_all_paths(graph, node, target, path)
            for newpath in newpaths:
                paths.append(newpath)

    return paths

def split_graph_properly(graph, split_at):
    """Split a given graph at specified level. For example, splitting the following graph at "S1":
                                (0, S)
                                   |
                                   |
                        (1, S1)--------(2, S1)   <----- Splitting at this level
                            |             |
                            |             |
                    (3, S2)---(4, S2)  (5, S2)

        Gives:
                (0, S)                  (0, S)
                   |                       |
                   |                       |
                (1, S1)         AND     (2, S1)
                   |                       |
                   |                       |
           (3, S2)---(4, S2)            (5, S2)


    Args:
        graph (dict): A dictionary representation of the graph, where:
                    key = Source node
                    value = Target node
                    (Nodes are represented as described elsewhere)
        split_at (str): Taxon level to split graph at

    Returns:
        subgraphs (list): List of subgraphs found after splitting intput graph at given level
    """
    subgraphs = []
    root = sorted(list(graph.keys()))[0]

    split_points = [k for k in graph.keys() if k[1] == split_at]
    if len(split_points) <= 1:
        logging.debug(msg=f"Can't split at {split_at}")
        return graph

    for split_pt in split_points:
        path_to_split_pt = find_all_paths(graph, root, split_pt)[0]
        if len(graph[split_pt]) != 0:
            post_split_leaves = graph[split_pt]
            path_to_split_pt.extend(post_split_leaves)
            for node in post_split_leaves:
                if len(graph[node]) != 0:
                    path_to_split_pt.extend(graph[node])
        final_path = sorted(path_to_split_pt)

        subgraphs.append(build_graph(final_path))

    return subgraphs

def decompose_graph(graph, known_subgraphs = []):
    """Split a complex graph into simple graphs.
        Here, a simple graph is defined as one where there is only one node at all levels except the lowest one.

    Args:
        graph (dict): A dictionary representation of the graph, where:
                    key = Source node
                    value = Target node
                    (Nodes are represented as described elsewhere)
        known_subgraphs (list, optional): List of known subgraphs to output. Defaults to [].

    Returns:
        subgraphs (list(dicts)): List of simple subgraphs decomposed from input graph
    """
    ## collect known subgraphs
    subgraphs = known_subgraphs

    ## get a list of taxon levels and find the highest level with count > 1
    levels_list = [k[1] for k in graph.keys()]
    for i in sorted(levels_list):
        if levels_list.count(i) > 1:
            multiplicity_level = i
            break

    ## dump graphs split at level identified above into temporary list
    ## for each subgraph check if it is simple, if not, call this function on it again
    tmp_subgraphs = split_graph_properly(graph, multiplicity_level)
    for g in tmp_subgraphs:
        if get_graph_complexity(g) == 0:
            subgraphs.append(g)
        else:
            logging.debug(f"Decomposing {g.keys()} further")
            ## recursion is cool!
            ## (but only when it works)
            cmp_subgraphs = decompose_graph(g, subgraphs)

    return subgraphs

def find_valid_graphs(graph, data_dict, threshold):
    """Given a graph, decopose into simple subgraphs if needed,
        then analyse it to find passing nodes and update the output dictionary

    Args:
        graph (dict): A dictionary representation of the graph, where:
                    key = Source node
                    value = Target node
                    (Nodes are represented as described elsewhere)
        data_dict (dict): Dictionary representation of kraken2 report
                    key = Indexed node, made of the kraken report index and the taxon level (eg (11, "S1"))
                    value = List where:
                            value[0] = Uniquely assigned reads for that entry
                            value[1] = Taxonomic ID of that entry
                            value[2] = Cumulative assigned reads, i.e. sum(value[0] for this entry + sum(value[0] of each child entry))
        threshold (int): Minimum number of reads required to consider an entry as non-spurious

    Returns:
        graph_meta (dict): A dict containing various info about processed graphs, where:
                                    key = target node
                                    value = dictionary containing info about path to target node
    """

    ## initialise output dict
    graph_meta = {}

    ## collect info for logging
    root = sorted(list(graph.keys()))[0]
    root_taxid = data_dict[root][1]

    ## check graph complexity
    graph_complexity = get_graph_complexity(graph)

    ## if graph is simple, pass as is into subgraphs
    ## else decompose graph into simple subgraphs
    if graph_complexity == 0:
        subgraphs = [graph]
    else:
        subgraphs = decompose_graph(graph)

    logging.info(msg=f"Found {len(subgraphs)} subgraphs in graph rooted at {root}, taxonomic ID: {root_taxid}.")

    chosen_subgraphs = []
    ## iterate over each subgraph
    ## check if there are any passing leaf nodes
    ## if yes, run polling on the passing leaf nodes and update output dict
    for subgraph in subgraphs:
        maximum_level = sorted(list(subgraph.keys()))[-1][1]
        terminals = [key for key in subgraph.keys() if key[1] == maximum_level and data_dict[key][0] >= threshold]
        if len(terminals) > 0:
            chosen_subgraphs.append(subgraph)
            taxa = []
            for term in terminals:
                taxa.append(data_dict[term][1])
            logging.debug(msg=f"Selected valid terminals {terminals}, taxonomic ID(s): {taxa} in graph rooted at {root}, root taxonomic ID: {root_taxid}.")


            graph_meta.update(update_graph_meta(subgraph, data_dict, terminals, root, root_taxid))

        ## if no passing leaf nodes, check for passing parent nodes
        ## in this case, check based on the cumulative num_reads for that parent
        ## i.e. if num_reads for parent and all its children > threshold, pass it
        else:
            penultimate = (TaxonLevel(maximum_level)-1).lvl
            penultimate_node = None
            ## find the parent node and check if it is valid
            for key in subgraph.keys():
                if key[1] == penultimate and data_dict[key][2] >= threshold:
                    penultimate_node = key
                    break

            if penultimate_node is not None:
                logging.debug(msg=f"{penultimate_node = }")
                pen_taxon = data_dict[penultimate_node][1]
                children = subgraph[penultimate_node]

                ## importantnly, check that the parent has > 1 children
                ## this catches cases where no leaves pass but their parents,often higher level taxa
                ### are passed leading to data duplication at sort_reads step
                if len(children) > 1:
                    chosen_subgraphs.append(subgraph)
                    logging.debug(msg=f"Selected valid parent node(s) {penultimate_node}, taxonomic ID(s): {pen_taxon} in graph rooted at {root}, root taxonomic ID: {root_taxid}")

                    ## run polling and update output using passing parent nodes
                    graph_meta.update(update_graph_meta(subgraph, data_dict, [penultimate_node], root, root_taxid, True))

                else:
                    logging.debug(msg=f"REJECTED {penultimate_node}, taxonomic ID(s): {pen_taxon} in graph rooted at {root}, root taxonomic ID: {root_taxid}")

            ## explicitly skip graphs where no leaves pass, nor parents
            ## this was somehow needed explicitly: weird
            else:
                continue

    return graph_meta

def update_graph_meta(valid_subgraph, data_dict, terminals, root, root_taxid, parent_selected = False):
    """Function to silo off:
        - polling
        - creating an entry for each valid result

    Args:
        valid_subgraph (dict): A dictionary representation of the graph, where:
                    key = Source node
                    value = Target node
                    (Nodes are represented as described elsewhere)
        data_dict (dict): Dictionary representation of kraken2 report
                    key = Indexed node, made of the kraken report index and the taxon level (eg (11, "S1"))
                    value = List where:
                            value[0] = Uniquely assigned reads for that entry
                            value[1] = Taxonomic ID of that entry
                            value[2] = Cumulative assigned reads, i.e. sum(value[0] for this entry + sum(value[0] of each child entry))
        terminals (list(tuples)): Passing terminal nodes in valid_subgraph
        root (tuple): Root node of valid_subgraph, used for logging.
        root_taxid (int): Taxonomic ID of root, used for logging
        parent_selected (bool, optional): Whether the terminals are leaf nodes or parent nodes. Defaults to False.

    Returns:
        graph_meta (dict): A dict containing various info about processed graphs, where:
                                    key = target node
                                    value = dictionary containing info about path to target nod
    """
    graph_meta = {}
    filt_end_nodes, pre_surprise, post_surprise = poll_leaves(end_nodes=terminals, data_dict=data_dict, parent_selected=parent_selected)
    if parent_selected:
        logging.info(msg=f"AFTER POLLING: Selected valid penultimate terminals {filt_end_nodes} in graph rooted at {root}, root taxonomic ID: {root_taxid}.")
        logging.info(msg=f"Entropy metrics | pre-poll: {pre_surprise} | post-poll: {post_surprise}")
    else:
        logging.info(msg=f"AFTER POLLING: Selected valid terminals {filt_end_nodes} in graph rooted at {root}, root taxonomic ID: {root_taxid}.")
        logging.info(msg=f"Entropy metrics | pre-poll: {pre_surprise} | post-poll: {post_surprise}")

    all_taxids_in_subgraph = [data_dict[key][1] for key in valid_subgraph.keys()]
    for end_node in set(filt_end_nodes):
        path_found = find_all_paths(valid_subgraph, root, end_node)[0]
        path_as_taxids = [data_dict[node][1] for node in path_found]
        graph_meta[data_dict[end_node][1]] = {
                                                "graph_idx": root[0],
                                                "source": root,
                                                "target": end_node,
                                                "parent_selected": parent_selected,
                                                "all_taxa": all_taxids_in_subgraph,
                                                "path": path_found,
                                                "path_as_taxids": path_as_taxids,
                                            }

    return graph_meta
