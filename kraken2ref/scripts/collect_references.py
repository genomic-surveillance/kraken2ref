import sys
import json
import argparse
import pandas as pd
from Bio import SeqIO

## read in inputs
def get_args():
    parser = argparse.ArgumentParser("Script to convert kraken2ref JSON to TSV report.")
    parser.add_argument("-i", "--input_json", type = str, required = True, help = "The JSON file produced by kraken2ref.")
    parser.add_argument("-l", "--library", type = str, required = True, help = "The kraken2 taxonomic report.")
    parser.add_argument("-m", "--name_map", type = str, required = True, help = "The kraken2 DB prelim_map.")


    return parser

def collect_references(in_file, library_fa, map_file):

    data_json = json.load(open(in_file))
    metadata = data_json["metadata"]
    outputs = data_json["outputs"]
    selected = [i for i in metadata["selected"]]
    is_condensed = metadata["summary"]["condense"]["condense_by_root"]
    if is_condensed:
        condense_info = metadata["summary"]["condense"]["condense_info"]

    # selected = [3464915, 3312104, 3673473, 3573866, 3793226, 3284012, 3645382, 3121653, 3332919, 114727]
    # condense_info = {2955291: ["3464915", "3312104", "3673473", "3573866", "3793226", "3284012", "3645382", "3121653", "3332919", "114727"]}
    # is_condensed = True

    prelim_map = pd.read_csv(map_file, sep = "\t", header = None)
    prelim_map_dict = dict(zip(prelim_map[2], prelim_map[1]))
    print(prelim_map_dict)
    output_map = {}
    not_found = []
    library = SeqIO.index(library_fa, "fasta")
    print(list(library.keys()))
    for taxid in selected:
        try:
            ref_name = prelim_map_dict[taxid]
            output_map[taxid] = library[ref_name]
        except KeyError as ke:
            not_found.append(taxid)
    sys.stderr.write(f"{not_found = }\n")

    if is_condensed:
        for k, v in condense_info.items():
            print(k)
            with open(f"ref_{str(k)}.fa", "w") as r:
                for i in v:
                    try:
                        r.write(output_map[int(i)].format("fasta"))
                    except KeyError as ke:
                        continue
    else:
        for k, v in output_map.items():
            with open(f"ref_{str(k)}.fa", "w") as r:
                if i not in not_found:
                    r.write(output_map[i].format("fasta"))
            r.close()

def main():
    args = get_args().parse_args()
    collect_references(args.input_json, args.library, args.name_map)


if __name__ == "__main__":
    main()









# library_file = "/Users/bd8/Developer/kraken2ref/scratch/collect_refs/Zydu0XGTj3.fna"
# library = SeqIO.index(library_file, "fasta")

# prelim_map = pd.read_csv("/Users/bd8/Developer/kraken2ref/scratch/collect_refs/taxid.map.txt", sep = "\t", header = None)
# prelim_map_dict = dict(zip(prelim_map[2], prelim_map[1]))
# print(prelim_map_dict[3597048])
# x = library["kraken:taxid|3464915|gb|HE584756"]
# print(x)
# selected = [3464915, 3312104, 3673473, 3573866, 3793226, 3284012, 3645382, 3121653, 3332919, 114727]
# condensed = {2955291: ["3464915", "3312104", "3673473", "3573866", "3793226", "3284012", "3645382", "3121653", "3332919", "114727"]}

# not_found = []
# for i in selected:
#     try:
#         ref_name = prelim_map_dict[i]
#         print(library[ref_name].format("fasta"))
#     except KeyError as ke:
#         not_found.append(i)
# print(f"{not_found = }")

def main():
    print("collecting refs")
