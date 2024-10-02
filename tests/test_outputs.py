import json
from kraken2ref.kraken2reference import KrakenProcessor

def test_singleton_output():
    singleton_proc = KrakenProcessor("singleton_report")
    singleton_proc.analyse_report(input_kraken_report_file="tests/artificial_reports/scov2_clean.report.txt", input_threshold=100, input_method="max", quiet=True)
    singleton_proc.write_output(prefix=".")

    outdata = json.load(open("singleton_report_decomposed.json", "r"))
    assert len(outdata["metadata"]["selected"]) == 1, "Wrong number of refs selected"
    assert outdata["metadata"]["selected"] == [2697049], "wrong ref selected"

def test_complex_tree_output():
    complex_tree_report_proc = KrakenProcessor("complex_tree")
    complex_tree_report_proc.analyse_report(input_kraken_report_file="tests/artificial_reports/fluA_clean.report.txt", input_threshold=100, input_method="max", quiet=True)
    complex_tree_report_proc.write_output(prefix=".")

    outdata = json.load(open("complex_tree_decomposed.json", "r"))
    assert len(outdata["metadata"]["selected"]) == 11, "Wrong number of refs selected"

    complex_tree_report_proc = KrakenProcessor("complex_tree_kmeans")
    complex_tree_report_proc.analyse_report(input_kraken_report_file="tests/artificial_reports/fluA_clean.report.txt", input_threshold=100, input_method="kmeans", quiet=True)
    complex_tree_report_proc.write_output(prefix=".")

    outdata_kmeans = json.load(open("complex_tree_kmeans_decomposed.json", "r"))
    assert len(outdata_kmeans["metadata"]["selected"]) == 12, "Wrong number of refs selected" ## one extra segment 1 ref should be selected AND written with kmeans method

test_complex_tree_output()