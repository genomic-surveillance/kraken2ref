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

    ## populate required data structures for graph building
    graph = {}
    plain_levels = []
    levels_as_taxa = []
    for indexed_node in indexed_nodes:
        graph[indexed_node] = []
        plain_levels.append(indexed_node[1])
        levels_as_taxa.append(TaxonLevel(indexed_node[1]))

    start = 0
    offset = indexed_nodes[0][0]
    while start < len(plain_levels)-1:
        end = start + 1
        left_node = TaxonLevel(plain_levels[start])
        right_node = TaxonLevel(plain_levels[end])

        ## ["S", "S1", "S2", "S3", "S3"]
        ## [(11,"S"), (12,"S1)", (13,"S2"), (14,"S3"), (15,"S3")]

        ## if left node less than right node: that's a path
        if left_node < right_node:
            graph[(start+offset,plain_levels[start])].append((end+offset, plain_levels[end]))

        ## if lef node equals right node: oops! find the closest parent
        if left_node >= right_node:
            sublist = levels_as_taxa[:end]
            parent = right_node - 1
            parent_idx = find_parent(sublist, parent)
            graph[(parent_idx+offset, parent.lvl)].append((end+offset, plain_levels[end]))
        start += 1

    return graph

def find_all_paths(graph, source, target, path=[]):
    """Function to find all paths from given source node to given target node.

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
    path = path + [source]
    if source == target:
        return [path]
    if source not in graph:
        return []
    paths = []
    for node in graph[source]:
        if node not in path:
            newpaths = find_all_paths(graph, node, target, path)
            for newpath in newpaths:
                paths.append(newpath)

    return paths

def get_end_points(graph, data_dict, threshold):
    paths = []
    path_as_tax_ids = []
    graph_nodes = [TaxonLevel(key) for (idx, key) in graph.keys()]
    maximum = max(graph_nodes)
    for (idx, lvl) in graph.keys():
        if lvl == "S":
            root = (idx, lvl)
            break
    end_nodes = [key for key in graph.keys() if key[1] == maximum.lvl]
    all_taxids_in_graph = [data_dict[key][1] for key in graph.keys()]
    for end_node in end_nodes:
        if data_dict[end_node][0] >= threshold:
            path_found = find_all_paths(graph, root, end_node)[0]
            tax_list = [data_dict[node][1] for node in path_found]
            # print("paths found =", path_found)
            # print(tax_list)
            paths.append(path_found)
            path_as_tax_ids.append(tax_list)
    # print("all paths in this graph = ", paths)
    # print("all tax_lists in this graph = ", tax_lists)
    return paths, path_as_tax_ids, all_taxids_in_graph


def get_graph_endpoints(graphs, data_dict, threshold):
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
        for end_node in end_nodes:
            path_found = find_all_paths(graph, root, end_node)[0]
            path_as_taxids = [data_dict[node][1] for node in path_found]
            graph_meta[end_node] = {"graph_idx": idx,
                                    "source": root,
                                    "all_taxa": all_taxids_in_graph,
                                    "path": path_found,
                                    "path_as_taxids": path_as_taxids,
                                    "taxid_filename": path_found[-1]
                                    }
    return graph_meta