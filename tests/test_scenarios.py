import pytest
from kraken2ref.kraken2reference import KrakenProcessor

def test_exit():

    def empty_report():
        empty_proc = KrakenProcessor("empty_report")
        empty_proc.analyse_report(input_kraken_report_file="tests/artificial_reports/empty_file.report.txt", input_threshold=100, input_method="max", quiet=True)

    def no_data_report():
        no_data_proc = KrakenProcessor("no_data_report")
        no_data_proc.analyse_report(input_kraken_report_file="tests/artificial_reports/no_data.report.txt", input_threshold=100, input_method="max", quiet=True)


    with pytest.raises(SystemExit) as pytest_wrapped_e:
        empty_report()
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        no_data_report()
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0

def test_singleton_report():
    singleton_proc = KrakenProcessor("singleton_report")
    singleton_proc.analyse_report(input_kraken_report_file="tests/artificial_reports/scov2_clean.report.txt", input_threshold=100, input_method="max", quiet=True)

    assert 2697049 in singleton_proc.tree_meta_out.keys(), "Wrong taxID appears to be selected, should be 2697049"
    assert singleton_proc.tree_meta_out[2697049]["graph_idx"] == 13, "Graph linked to wrong line in file, graph idx should be 13"
    assert singleton_proc.tree_meta_out[2697049]["source_taxid"] == 694009, "Wrong parent reported, should be 694009"
    assert singleton_proc.tree_meta_out[2697049]["path"] == [(13, 'S'), (14, 'S1')], "Wrong path reported, should be [(13, 'S'), (14, 'S1')]"

def test_complex_tree_report():
    complex_tree_report_proc = KrakenProcessor("complex_tree")
    complex_tree_report_proc.analyse_report(input_kraken_report_file="tests/artificial_reports/fluA_clean.report.txt", input_threshold=100, input_method="max", quiet=True)

    assert sorted(complex_tree_report_proc.tree_meta_out.keys()) == [3149250, 3149251, 3149253, 3185065, 3219756, 3219759, 3285787, 3290795, 3369048, 3514037, 3847902], \
    "Wrong taxIDs appear to be selected, should be [3149250, 3149251, 3149253, 3185065, 3219756, 3219759, 3285787, 3290795, 3369048, 3514037, 3847902]"

def test_l1_tree_report():
    l1_tree_report_max_proc = KrakenProcessor("l1_tree_max")
    l1_tree_report_max_proc.analyse_report(input_kraken_report_file="tests/artificial_reports/adenovirus_clean.report.txt", input_threshold=100, input_method="max", quiet=True)

    assert sorted(l1_tree_report_max_proc.tree_meta_out.keys()) == [10519, 28285], "Wrong taxIDs appear to be selected in mode MAX, should be [10519, 28285]"

    l1_tree_report_skew_proc = KrakenProcessor("l1_tree_kmeans")
    l1_tree_report_skew_proc.analyse_report(input_kraken_report_file="tests/artificial_reports/adenovirus_clean.report.txt", input_threshold=100, input_method="skew", quiet=True)

    assert sorted(l1_tree_report_skew_proc.tree_meta_out.keys()) == [10519, 10533, 28285], "Wrong taxIDs appear to be selected in mode MAX, should be [10519, 28285, 28285]"
