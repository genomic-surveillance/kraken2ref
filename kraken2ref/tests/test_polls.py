from kraken2ref.src.polling_functions import poll_leaves

def __init__():
    end_nodes = [(11, "S3"), (18, "S3"), (19, "S3"), (22, "S3"), (31, "S3"), (68, "S3")]
    data_dict_1 = {(11, "S3"): (45, "Tax1"),
                (18, "S3"): (52, "Tax2"),
                (19, "S3"): (11, "Tax3"),
                (22, "S3"): (5, "Tax4"),
                (31, "S3"): (88, "Tax5"),
                (68, "S3"): (32, "Tax6")
                }

    data_dict_2 = {(11, "S3"): (4, "Tax1"),
                (18, "S3"): (5, "Tax2"),
                (19, "S3"): (11, "Tax3"),
                (22, "S3"): (5, "Tax4"),
                (31, "S3"): (88, "Tax5"),
                (68, "S3"): (63, "Tax6")
                }

    end_nodes_3 = [(11, "S3"), (18, "S3"), (19, "S3"), (22, "S3"), (31, "S3"), (68, "S3"), (69, "S3"), (70, "S3")]
    data_dict_3 = {(11, "S3"): (916211, "Tax1"),
                (18, "S3"): (718948, "Tax2"),
                (19, "S3"): (523656, "Tax3"),
                (22, "S3"): (491852, "Tax4"),
                (31, "S3"): (360002, "Tax5"),
                (68, "S3"): (308869, "Tax6"),
                (69, "S3"): (235529, "Tax7"),
                (70, "S3"): (201223, "Tax8")
                }

    end_nodes_4 = [(11, "S3"), (18, "S3"), (19, "S3")]
    data_dict_4 = {(11, "S3"): (916211, "Tax1"),
                (18, "S3"): (718948, "Tax2"),
                (19, "S3"): (523656, "Tax3")
                }

    expected_outputs = {"data_dict_1": [end_nodes, data_dict_1, [(11, 'S3'), (18, 'S3'), (31, 'S3'), (68, 'S3')]],
                         "data_dict_2": [end_nodes, data_dict_2, [(31, 'S3'), (68, 'S3')]],
                         "data_dict_3": [end_nodes_3, data_dict_3, [(11, 'S3'), (18, 'S3'), (19, 'S3'), (22, 'S3'), (31, 'S3'), (68, 'S3')]],
                         "data_dict_4": [end_nodes_4, data_dict_4, [(11, 'S3'), (18, 'S3'), (19, 'S3')]]}

    return expected_outputs

def test_polling():
    expected_outputs = __init__()

    for k, [n, d, e] in expected_outputs.items():
        filt_nodes, prefilt_surprise, postfilt_surprise = poll_leaves(n, d)
        try:
            assert e == sorted(filt_nodes)
        except AssertionError as ae:
            print(k)
            print(f"for {k}: found {filt_nodes}, expected {e}")

