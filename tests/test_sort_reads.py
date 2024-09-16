from kraken2ref import sort_reads
import os, json
import pandas as pd

from kraken2ref.kraken2reference import KrakenProcessor

def __init__():
    my_processor = KrakenProcessor("test_sort")
    my_processor.analyse_report(input_kraken_report_file="tests/test_set/with_fluseg_db/testset_50x.uniq.report.txt", input_threshold=100, input_method="kmeans", quiet=True)
    my_processor.write_output(prefix = ".", suffix="decomposed")


def test_tree():
    __init__()
    sort_reads.sort_reads(
            sample_id="test_sort",
            kraken_output="tests/test_set/with_fluseg_db/testset_50x.kraken.output",
            mode="tree",
            condense=False,
            update_output=False,
            taxon_list=None,
            ref_json_file="test_sort_decomposed.json",
            outdir=".")

    with open("./test_sort_updated_decomposed.json", "r") as data_json:
        data = json.load(data_json)
        per_taxon = data["metadata"]["summary"]["per_taxon"]

        numreads_df = open("tests/test_set/with_fluseg_db/tree_numreads.txt")
        r1_dict = {}
        r2_dict = {}
        for line in numreads_df.readlines():
            decomp = line.split("\t")
            if "R1" in decomp[0]:
                tax = decomp[0].replace("_R1.fq", "")
                r1_dict[tax] = int(decomp[1])
            if "R2" in decomp[0]:
                tax = decomp[0].replace("_R2.fq", "")
                r2_dict[tax] = int(decomp[1])

        for k, v in per_taxon.items():
            assert v == r1_dict[k]
            assert v == r2_dict[k]

    # os.system("rm ./*.fq")
    os.system("rm ./test_sort_*")


# def test_unique():
#     __init__()
#     sort_reads.sort_reads(sample_id="test_sort",
#         kraken_output="tests/test_set/with_fluseg_db/testset_50x.kraken.output",
#         mode="unique",
#         fastq1="tests/test_set/with_fluseg_db/testset_50x.class_seqs_1.fq",
#         fastq2="tests/test_set/with_fluseg_db/testset_50x.class_seqs_2.fq",
#         condense=False,
#         update_output=False,
#         taxon_list="2697049",
#         ref_json_file="test_sort_decomposed.json",
#         outdir=".")

#     r1_reads = int(os.popen("grep ^@ 2697049_R1.fq | wc -l").read().strip())
#     print(r1_reads)
#     r2_reads = int(os.popen("grep ^@ 2697049_R2.fq | wc -l").read().strip())
#     assert r1_reads == 40394
#     assert r2_reads == 40394

#     os.system("rm ./*.fq")
#     os.system("rm ./test_sort_*")

def test_condensed():
    __init__()
    sort_reads.sort_reads(sample_id="test_sort",
        kraken_output="tests/test_set/with_fluseg_db/testset_50x.kraken.output",
        mode="tree",
        condense=True,
        update_output=False,
        taxon_list=None,
        ref_json_file="test_sort_decomposed.json",
        outdir=".")


    with open("./test_sort_updated_decomposed.json", "r") as data_json:
        data = json.load(data_json)
        per_taxon = data["metadata"]["summary"]["per_taxon"]

        numreads_df = open("tests/test_set/with_fluseg_db/condensed_numreads.txt")
        r1_dict = {}
        r2_dict = {}
        for line in numreads_df.readlines():
            decomp = line.split("\t")
            if "R1" in decomp[0]:
                tax = decomp[0].replace("_R1.fq", "")
                r1_dict[tax] = int(decomp[1])
            if "R2" in decomp[0]:
                tax = decomp[0].replace("_R2.fq", "")
                r2_dict[tax] = int(decomp[1])

    for k, v in per_taxon.items():
        assert v == r1_dict[k]
        assert v == r2_dict[k]

    os.system("rm ./*.fq")
    os.system("rm ./test_sort_*")
