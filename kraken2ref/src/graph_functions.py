from kraken2ref.src.taxonlevel import TaxonLevel, find_parent

def build_graph(indexed_nodes):

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

def find_all_paths(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return [path]
    if start not in graph:
        return []
    paths = []
    for node in graph[start]:
        if node not in path:
            newpaths = find_all_paths(graph, node, end, path)
            for newpath in newpaths:
                paths.append(newpath)

    return paths

def get_end_points(graph, data_dict, threshold):
    paths = []
    tax_lists = []
    graph_nodes = [TaxonLevel(key) for (idx, key) in graph.keys()]
    maximum = max(graph_nodes)
    for (idx, lvl) in graph.keys():
        if lvl == "S":
            root = (idx, lvl)
            break
    end_nodes = [key for key in graph.keys() if key[1] == maximum.lvl]
    for end_node in end_nodes:
        if data_dict[end_node][0] >= threshold:
            path_found = find_all_paths(graph, root, end_node)[0]
            tax_list = [data_dict[node][1] for node in path_found]
            # print("paths found =", path_found)
            # print(tax_list)
            paths.append(path_found)
            tax_lists.append(tax_list)
    # print("all paths in this graph = ", paths)
    # print("all tax_lists in this graph = ", tax_lists)
    return paths, tax_lists