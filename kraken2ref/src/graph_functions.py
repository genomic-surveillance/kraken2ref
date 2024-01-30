import sys
from kraken2ref.src.taxonlevel import TaxonLevel, find_parent

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
    offset = indexed_nodes[0][0]

    ## trawl over list to find edges in graph
    while start < len(plain_levels)-1:
        end = start + 1
        left_node = TaxonLevel(plain_levels[start])
        right_node = TaxonLevel(plain_levels[end])

        ## if left node less than right node: that's an edge
        if left_node < right_node:
            graph[(start+offset,plain_levels[start])].append((end+offset, plain_levels[end]))

        ## if left node equals right node: oops! find the closest parent
        if left_node >= right_node:
            sublist = levels_as_taxa[:end]
            parent = right_node - 1
            parent_idx = find_parent(sublist, parent)
            graph[(parent_idx+offset, parent.lvl)].append((end+offset, plain_levels[end]))
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

    ## paths is a list in case of multiple paths from ource to target
    paths = []

    ## iterate over graph dict (i.e. edges)
    for node in graph[source]:
        if node not in path:
            ## recursively find all "sub-paths" and update paths list
            newpaths = find_all_paths(graph, node, target, path)
            for newpath in newpaths:
                paths.append(newpath)

    return paths

def get_graph_endpoints(graphs, data_dict, threshold):
    """Operate on list of graphs to collect key information

    Args:
        graphs (list(dicts)): A list of graphs to collect info on; each graph is a dict describing
                                edges between nodes in that graph
        data_dict (dict): Dictionary representation of kraken2 report
        threshold (int): Minimum number of read needed to retain a leaf node

    Returns:
        graph_meta (dict(dicts)): A dict containing various info about processed graphs, where:
                                    key = target node
                                    value = dictionary containing info about path to target node
    """
    graph_meta = {}
    for idx, graph in enumerate(graphs):
        all_taxids_in_graph = [data_dict[key][1] for key in graph.keys()]
        graph_nodes = [TaxonLevel(key) for (idx, key) in graph.keys()]
        maximum = max(graph_nodes)
        for (idx, lvl) in graph.keys():
            if lvl == "S":
                root = (idx, lvl)
                break
        end_nodes = [key for key in graph.keys() if key[1] == maximum.lvl and data_dict[key][0] >= threshold]
        ## handle exception where end_nodes is empty
        if len(end_nodes) == 0:
            ## create a unique "undetermined target" entry with no path; use parent as filename
            ## jump up 1 level to get parent
            parent_level = maximum - 1

            ## find parent node
            for (idx,lvl) in graph.keys():
                if lvl == parent_level.lvl:
                    parent_node = (idx,lvl)
                    break

            ## write error: this is important
            sys.stderr.write(f"\n\nNoSuitableTargetError: No leaf nodes found suitable, reverting to pre-selected reference for node: {parent_node}, taxonomic ID: {data_dict[parent_node][1]}\n\n")

            ## populate graph_meta
            graph_meta["parent_selected_"+str(idx)] = {"graph_idx": idx,
                                    "source": root,
                                    "all_taxa": all_taxids_in_graph,
                                    "path": None,
                                    "path_as_taxids": None,
                                    "taxid_filename": data_dict[parent_node][1]
                                    }
        else:
            ##TODO: add polling function here, apply to end nodes list
            for end_node in end_nodes:
                path_found = find_all_paths(graph, root, end_node)[0]
                path_as_taxids = [data_dict[node][1] for node in path_found]
                graph_meta[end_node] = {"graph_idx": idx,
                                        "source": root,
                                        "all_taxa": all_taxids_in_graph,
                                        "path": path_found,
                                        "path_as_taxids": path_as_taxids,
                                        "taxid_filename": path_as_taxids[-1]
                                        }

    return graph_meta

def poll_leaves(end_nodes, data_dict):
    filtered_end_nodes = []

    return filtered_end_nodes