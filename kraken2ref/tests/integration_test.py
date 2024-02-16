import json
import os
from kraken2ref.src.kraken_taxonomy_report import KrakenTaxonomyReport

def __init__():
    RUN_DIR = os.getcwd()
    test_tax_report = KrakenTaxonomyReport(sample_id="int_test")

    assert test_tax_report
    return test_tax_report



def test_basic():

    test_tax_report = __init__()
    test_tax_report.pick_reference_taxid(in_file="kraken2ref/tests/fixtures/report.txt", outdir="./", min_abs_reads=3)

    tmp = open("int_test_decomposed.json")
    data = json.load(tmp)
    assert sorted(list(data.keys())) == sorted(['metadata', 'outputs'])
    assert sorted(list(data["outputs"].keys())) == sorted(['335341', '518987', '11250', '12814', '2697049', 'parent_selected_694009_T'])

def test_split():

    test_tax_report = __init__()
    test_tax_report.pick_reference_taxid(in_file="kraken2ref/tests/fixtures/report.txt", outdir="./", min_abs_reads=3, split_at = "S2")

    tmp = open("int_test_decomposed.json")
    data = json.load(tmp)
    assert list(data.keys()) == ['metadata', 'outputs']
    assert list(data["outputs"].keys()) == ['335341', '641809', 'parent_selected_102797', 'parent_selected_114728',
                                            'parent_selected_102793', 'parent_selected_102796', 'parent_selected_102800', 'parent_selected_119212',
                                            'parent_selected_222769', 'parent_selected_140020', 'parent_selected_119218', 'parent_selected_11320',
                                            'parent_selected_384619', 'parent_selected_119211', 'parent_selected_184006', 'parent_selected_437442',
                                            '518987', 'parent_selected_11520', '11250', '12814', '2697049', 'parent_selected_694009_T']
    return

