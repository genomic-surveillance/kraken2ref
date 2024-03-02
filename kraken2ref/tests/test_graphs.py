from kraken2ref.src.graph_functions import build_graph, split_graph_properly, find_all_paths

def __init__():
    test_input = [(11,"S"), (12,"S1"), (13,"S2"), (14,"S3"), (15,"S4"), (16,"S3"), (17,"S2"), (18,"S3"), (19,"S3")]
    test_output = {(11, 'S'): [(12, 'S1')], (12, 'S1'): [(13, 'S2'), (17, 'S2')], (13, 'S2'): [(14, 'S3'), (16, 'S3')], (14, 'S3'): [(15, 'S4')], (15, 'S4'): [], (16, 'S3'): [], (17, 'S2'): [(18, 'S3'), (19, 'S3')], (18, 'S3'): [], (19, 'S3'): []}
    test_split_output = [{(11, 'S'): [(12, 'S1')], (12, 'S1'): [(13, 'S2')], (13, 'S2'): [(14, 'S3'), (16, 'S3')], (14, 'S3'): [(15, 'S4')], (15, 'S4'): [], (16, 'S3'): []},
                            {(11, 'S'): [(12, 'S1')], (12, 'S1'): [(17, 'S2')], (17, 'S2'): [(18, 'S3'), (19, 'S3')], (18, 'S3'): [], (19, 'S3'): []}]
    return test_input, test_output, test_split_output

def test_graph():
    test_input, test_output, test_split_output = __init__()
    assert build_graph(test_input) == test_output

def test_splitting():
    test_input, test_output, test_split_output = __init__()
    graph = build_graph(test_input)
    split = split_graph_properly(graph, "S2")
    assert split == test_split_output

def test_path_find():
    test_input, test_output, test_split_output = __init__()
    graph = build_graph(test_input)
    path = find_all_paths(graph, (11, "S"), (19, "S3"))
    assert path[0] == [(11, 'S'), (12, 'S1'), (17, 'S2'), (19, 'S3')]

