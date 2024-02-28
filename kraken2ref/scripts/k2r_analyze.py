from kraken2ref.src.graph_functions import build_graph, split_graph
from kraken2ref.src.kraken_taxonomy_report import KrakenTaxonomyReport
from kraken2ref.src.taxonlevel import TaxonLevel
from kraken2ref.src.polling_functions import poll_leaves
import logging

expect = [ # With threshold = 350
    3464915, # Segment 4 H1 - S4
    3312104, # Segment 1 - S3
    3673473, # Segment 2 - S3
    3573866, 3793226, # Segment 3 - S3
    3284012, # Segment 5 - S3
    3645382, # Segment 7 - S3
    3121653, # Segment 6 N1 - S3
    3332919, # Segment 8 - S3
    114727, # PARENT_SELECTED - S2
]
threshold = 350
my_report = KrakenTaxonomyReport("new_logic")
all_nodes, data_dict = my_report.read_kraken_report(kraken_report_file="/Users/bd8/Developer/kraken2ref/kraken_reports/48367_2_1_Flu_H1N1_100K_P1.report.txt")
graphs = []
for node_list in all_nodes:
    graphs.append(build_graph(node_list))

# complex_graphs = []
# for graph in graphs:
#     root = min(graph.keys())
#     analysis = {root: {"terminal": [], "graph": None} }
#     terminal_nodes = set()
#     for k, v in graph.items():
#         # print(f"{k}: {v}")
#         if len(v) == 0:
#             analysis[root]["terminal"].append(k)
#             terminal_nodes.add(k[1])
#     print(f"{analysis = }")
#     if len(terminal_nodes) > 1:
#         print(f"OH NO!!! GRAPH ROOTED AT {root} IS COMPLICATED ON SO MANY LEVELS: {terminal_nodes}!!!")
#         analysis[root]["graph"] = graph
#         complex_graphs.append()

# print(complex_graphs)


def get_graph_complexity(graph):
    root = min(graph.keys())
    terminal_levels = set()
    terminal_nodes = []
    for k, v in graph.items():
        if len(v) == 0:
            terminal_nodes.append(k)
            terminal_levels.add(k[1])
    if len(terminal_levels) > 1:
        print(f"OH NO!!! GRAPH ROOTED AT {root} IS COMPLICATED ON SO MANY LEVELS: {terminal_levels}!!!")
        return True, terminal_nodes, terminal_levels
    else:
        return False, terminal_nodes, list(terminal_levels)[0]

multilevel_graph = []
simple_graphs = []
broad_graphs = []
for graph in graphs:
    graph_is_complex, terminal_nodes, terminal_levels = get_graph_complexity(graph)
    if graph_is_complex:
        multilevel_graph.append(graph)
    else:
        breadth_indicator_level = (TaxonLevel(terminal_levels) - 1).lvl
        # print(f"{breadth_indicator_level = }")
        indicator_nodes = [node for node in graph.keys() if node[1] == breadth_indicator_level]
        # print(f"{indicator_nodes = }")
        if len(indicator_nodes) > 1:
            broad_graphs.append(graph)
        else:
            simple_graphs.append(graph)

def process_simple_graph(graph, data_dict, threshold):
    ## graph_entry looks like: {root: [graph, terminal_nodes]}
    root = min(sorted(graph.keys()))
    root_tax = data_dict[root][1]
    root_cumulative_hits = data_dict[root][2]
    if root_cumulative_hits < threshold:
        logging.info(msg=f"Graph rooted at {root} (taxonomic ID: {root_tax}) does not contain enough reads ({root_cumulative_hits}).")
        print(f"Graph rooted at {root} (taxonomic ID: {root_tax}) does not contain enough reads ({root_cumulative_hits}).")
        return None
    terminals = [key for key in graph.keys() if len(graph[key]) == 0 and data_dict[key][0] >= threshold]
    print(f"{terminals = }")
    if len(terminals) == 0:
        logging.info(msg=f"Graph rooted at {root} (taxonomic ID: {root_tax}) contains no suitable targets.")
        print(f"Graph rooted at {root} (taxonomic ID: {root_tax}) contains no suitable targets.")
        return None
    if len(terminals) == 1:
        return list(graph.keys())
    else:
        filt_end_nodes, pre_surprise, post_surprise = poll_leaves(end_nodes=terminals, data_dict=data_dict)
        return filt_end_nodes

def process_broad_graph(graph, data_dict, threshold, split_at):
    terminals = [key for key in graph.keys() if len(graph[key]) == 0 and data_dict[key][0] >= threshold]
    parent_level = (TaxonLevel(terminals[0][1]) - 1).lvl
    parent_nodes = [k for k in graph.key() if k[1] == parent_level]
    to_remove = []
    for parent in parent_nodes:
        if data_dict[parent][2] < threshold:
            if data_dict[parent][0] < threshold:
                to_remove.append(parent)
                to_remove.extend(graph[parent])
                logging.info(msg=f"Removing sub-tree under {parent} (taxonomic ID: {data_dict[parent][1]}) due to insufficient reads {data_dict[parent][2]}.")
    for node_to_remove in to_remove:
        del graph[node_to_remove]
        for v in graph.values():
            try:
                v.remove(node_to_remove)
            except ValueError as ve:
                pass

    subgraphs = split_graph(graph, split_at)


    return


for g in simple_graphs:
    x = process_simple_graph(g, data_dict, 350)
    print(x)
