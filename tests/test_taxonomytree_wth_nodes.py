from kraken2ref.taxonomytree import TaxonomyTree

def __init__():
    test_input = [(11,"S"), (12,"S1"), (13,"S2"), (14,"S3"), (15,"S4"), (16,"S3"), (17,"S2"), (18,"S3")]
    test_output = {(11, 'S'): [(12, 'S1')],
                   (12, 'S1'): [(13, 'S2'), (17, 'S2')],
                   (13, 'S2'): [(14, 'S3'), (16, 'S3')],
                   (14, 'S3'): [(15, 'S4')],
                   (15, 'S4'): [],
                   (16, 'S3'): [],
                   (17, 'S2'): [(18, 'S3')],
                   (18, 'S3'): [],
                   }
    test_split_output = [
                        {(11, 'S'): [(12, 'S1')],
                          (12, 'S1'): [(13, 'S2')],
                          (13, 'S2'): [(14, 'S3'),(16, 'S3')],
                          (14, 'S3'): [(15, 'S4')],
                          (15, 'S4'): [],
                          (16, 'S3'): []},

                        {(11, 'S'): [(12, 'S1')],
                         (12, 'S1'): [(17, 'S2')],
                         (17, 'S2'): [(18, 'S3'), ],
                         (18, 'S3'): [],
                         }
                        ]
    test_decomposed_output = [
                                {(11, 'S'): [(12, 'S1')],
                                (12, 'S1'): [(13, 'S2')],
                                (13, 'S2'): [(14, 'S3')],
                                (14, 'S3'): [(15, 'S4')],
                                (15, 'S4'): []},

                                {(11, 'S'): [(12, 'S1')],
                                (12, 'S1'): [(13, 'S2')],
                                (13, 'S2'): [(16, 'S3')],
                                (16, 'S3'): []},

                                {(11, 'S'): [(12, 'S1')],
                                (12, 'S1'): [(17, 'S2')],
                                (17, 'S2'): [(18, 'S3'), ],
                                (18, 'S3'): [],
                                }
                            ]
    return test_input, test_output, test_split_output, test_decomposed_output


def test_taxonomy_tree_with_nodes():
    test_input, test_output, test_split_output, test_decomposed_output  = __init__()

    test_tree_from_nodes = TaxonomyTree(nodes = test_input)

    ## check attributes
    assert test_tree_from_nodes.nodes == test_input
    assert test_tree_from_nodes.max_lvl == "S4"
    assert test_tree_from_nodes.root == (11, "S")
    assert test_tree_from_nodes.root_idx == 11
    assert sorted(test_tree_from_nodes.leaf_lvls) == ["S3", "S4"]
    assert sorted(test_tree_from_nodes.leaf_nodes) == [(15, 'S4'), (16, 'S3'), (18, 'S3')]
    assert sorted(test_tree_from_nodes.subterminal_lvls) == ["S2", "S3"]
    assert sorted(test_tree_from_nodes.subterminal_nodes) == [(13, 'S2'), (14, 'S3'), (17, 'S2')]

    ## check functions
    ## test build_from
    assert test_tree_from_nodes.graph == test_output
    ## check split_at
    assert test_tree_from_nodes.split_tree(test_tree_from_nodes.graph, "S2") == test_split_output
    ## check decompose_tree
    # test_tree_from_nodes.decompose_tree(test_tree_from_nodes.graph)
    assert test_tree_from_nodes.subgraphs == test_decomposed_output
    ## check get_tree_complexity
    assert test_tree_from_nodes.complexity == 2

