import json
from Bio import SeqIO
from kraken2ref import sort_reads
from kraken2ref.dump_fastqs import dump_to_files
from kraken2ref.kraken2reference import KrakenProcessor

def test_dump_basic(tmp_path):
    dump_basic_proc = KrakenProcessor("test_dump_basic")
    dump_basic_proc.analyse_report(input_kraken_report_file="tests/artificial_reports/adenovirus_clean.report.txt", input_threshold=100, input_method="max", quiet=True)
    dump_basic_proc.write_output(prefix = tmp_path, suffix="decomposed")

    sort_reads.sort_reads(sample_id="test_dump_basic",
                          kraken_output="tests/test_set/with_fluseg_db/testset_50x.kraken.output",
                          mode="tree",
                          ref_json_file=f"{tmp_path}/test_dump_basic_decomposed.json",
                          outdir=tmp_path,
                          update_output=False,
                          condense=False)

    tax_to_read_ids_here = json.load(open(f"{tmp_path}/test_dump_basic_tax_to_reads.json", "r"))
    dump_to_files(sample_id="test_dump_basic",
        tax_to_readids_dict=tax_to_read_ids_here,
        fq1="tests/test_set/with_fluseg_db/testset_50x_1.fq",
        fq2="tests/test_set/with_fluseg_db/testset_50x_2.fq",
        outdir=tmp_path)

    files_expected = {f"{tmp_path}/test_dump_basic_10519_R1.fq": 5900, f"{tmp_path}/test_dump_basic_10519_R2.fq": 5900, f"{tmp_path}/test_dump_basic_28285_R1.fq": 17950, f"{tmp_path}/test_dump_basic_28285_R2.fq": 17950}
    for filepath, num_reads_expected in files_expected.items():
        recs = list(SeqIO.parse(filepath, "fastq"))
        assert len(recs) == num_reads_expected, f"File {filepath} does not contain the expected number of reads {num_reads_expected}"

