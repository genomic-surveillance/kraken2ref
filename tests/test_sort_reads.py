import json

from kraken2ref import sort_reads
from kraken2ref.kraken2reference import KrakenProcessor

def test_sort_basic():
    sort_basic_proc = KrakenProcessor("test_sort_basic")
    sort_basic_proc.analyse_report(input_kraken_report_file="tests/artificial_reports/adenovirus_clean.report.txt", input_threshold=100, input_method="max", quiet=True)
    sort_basic_proc.write_output(prefix = ".", suffix="decomposed")

    sort_reads.sort_reads(sample_id="test_sort_basic",
                          kraken_output="tests/test_set/with_fluseg_db/testset_50x.kraken.output",
                          mode="tree",
                          ref_json_file="test_sort_basic_decomposed.json",
                          outdir=".",
                          update_output=False,
                          condense=False)

    outdata = json.load(open("test_sort_basic_tax_to_reads.json", "r"))
    assert sorted(outdata.keys()) == ["10519", "28285"], "Wrong taxIDs appear to be in output"
    assert len(outdata["10519"]) == 5900, "Wrong number of reads found"
    assert len(outdata["28285"]) == 17950, "Wrong number of reads found"

    with open("test_sort_basic_unwritten_reads.txt", "r") as unwritten:
        num_lines = len(unwritten.readlines())
    assert num_lines == 74019



