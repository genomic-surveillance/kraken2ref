import sys
import json
import datetime
import pandas as pd
from Bio import SeqIO

def write_fastq(fq1, fq2, kraken_out, update_output, ref_data):

    NOW = datetime.datetime.now()
    fq1_dict = SeqIO.index(fq1, "fastq")
    fq2_dict = SeqIO.index(fq2, "fastq")
    k = list(fq1_dict.keys())[0]
    if k[-2:] == "/1":
        slashes = True
    else:
        slashes = False

    read_data = pd.read_csv(kraken_out, sep = "\t", header = None)
    num_all = len(read_data)
    read_data = read_data[read_data[0] == "C"]
    num_class = len(read_data)
    print(len(read_data))
    ## keys are taxid; values are reads assigned to the key
    read_dict = {}
    for i, row in read_data.iterrows():
        if str(row[2]) not in read_dict.keys():
            read_dict[str(row[2])] = [row[1]]
        else:
            read_dict[str(row[2])].append(row[1])

    ## keys are selected refs; values are taxa included in this ref, and empty list for reads
    ref_map = {k: [set(v["all_taxa"]), []] for k, v in json.load(open(ref_data))["outputs"].items()}

    for tax_id in read_dict.keys():
        for ref_tax, vals in ref_map.items():
            if tax_id in vals[0]:
                vals[1].extend(read_dict[tax_id])

    to_remove = []
    for ref_tax, [v, reads] in ref_map.items():
        if len(reads) == 0:
            to_remove.append(ref_tax)
    for tax in to_remove:
        del ref_map[tax]

    file_handle_map = {k: [open(k+"_R1.fq", "w"), open(k+"_R2.fq", "w")] for k in ref_map.keys()}
    file_read_counts = {k: 0 for k in ref_map.keys()}
    wrote = 0
    for ref_tax, [v, reads] in ref_map.items():
        if len(reads) == 0:
            continue
        if slashes:
            for read in reads:
                r1_key = read+"/1"
                r2_key = read+"/2"
                file_handle_map[ref_tax][0].write(fq1_dict[r1_key].format("fastq"))
                file_handle_map[ref_tax][1].write(fq2_dict[r2_key].format("fastq"))
                file_read_counts[ref_tax] += 1
                wrote += 1
        else:
            for read in reads:
                file_handle_map[ref_tax][0].write(fq1_dict[read].format("fastq"))
                file_handle_map[ref_tax][1].write(fq2_dict[read].format("fastq"))
                wrote += 1
                file_read_counts[ref_tax] += 1

    # for k in file_handle_map.keys():
    #     k[0].close()
    #     k[1].close()

    summary = {
        "info": {
            "time_updated": str(NOW),
            "total_input_reads": num_all,
            "total_classified_reads": num_class,
            "non_unique_writes": wrote
                    },
        "per_taxon": file_read_counts
    }


    with open(ref_data, "r+") as data_json:
        data = json.load(data_json)
        data["metadata"]["summary"] = summary
        if update_output:
            data_json.seek(0)
            json.dump(data, data_json, indent=4)
        else:
            new_json_path = ref_data.replace("_decomposed.json", "_updated_decomposed.json")
            with open(new_json_path, "w") as new_json:
                json.dump(data, new_json, indent=4)

    sys.stderr.write(f"Wrote {wrote} non-unique read-pairs to {len(file_read_counts.keys())} file-pairs.\n\n")

    return summary