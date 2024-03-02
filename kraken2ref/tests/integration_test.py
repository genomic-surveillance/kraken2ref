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
    test_tax_report.pick_reference_taxid(in_file="kraken2ref/tests/fixtures/report.txt", outdir="./", min_abs_reads=25)

    tmp = open("int_test_decomposed.json")
    data = json.load(tmp)
    print(sorted(list(data["outputs"].keys())))
    assert sorted(list(data.keys())) == sorted(['metadata', 'outputs'])
    assert sorted([int(i) for i in data["outputs"].keys()]) == sorted([335341, 114727, 518987, 11250, 12814, 2697049])


