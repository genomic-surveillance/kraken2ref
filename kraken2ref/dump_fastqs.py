import io, os, sys
import json
import logging
from collections import defaultdict
from Bio import SeqIO

def dump_to_files(sample_id, tax_to_readids_dict, fq1, fq2, outdir,
        buffer_size=io.DEFAULT_BUFFER_SIZE):

    # Build readid -> taxid map (fast lookup)
    readids_to_taxids: Dict[str, Set[str]] = {}
    for taxid, read_list in tax_to_readids_dict.items():        
        for rid in read_list:
            taxdict = readids_to_taxids.get(rid)
            if taxdict is None:
                readids_to_taxids[rid] = {taxid}
            else:
                readids_to_taxids[rid].add(taxid)

    # Open all output file handles once
    outputs = {
        taxid: (
            open(os.path.join(outdir, f"{sample_id}_{taxid}_R1.fq"), "w", buffering=buffer_size),
            open(os.path.join(outdir, f"{sample_id}_{taxid}_R2.fq"), "w", buffering=buffer_size)
        )
        for taxid in tax_to_readids_dict
    }

    it1 = SeqIO.parse(fq1, "fastq")
    it2 = SeqIO.parse(fq2, "fastq")

    # Pass through input fastq, funnelling read pairs to appropriate output file(s)
    for r1, r2 in zip(it1, it2):
        rid = r1.id
        if rid.endswith("/1"):
            rid = rid[:-2]

        taxids = readids_to_taxids.get(rid)

        if taxids is None:
            continue

        for taxid in taxids:
            R1, R2 = outputs[taxid]
            R1.write(r1.format("fastq"))
            R2.write(r2.format("fastq"))

    # Close everything
    for R1, R2 in outputs.values():
        R1.close()
        R2.close()

def dump_fastqs(args):
    sample_id = args.sample_id
    json_tax_to_readsid_path = args.tax_to_readsid_path
    fq1 = args.fastq1
    fq2 = args.fastq2
    outdir = args.outdir
    buffer_size = args.buffer_size

    # load tax to reads id dictionary
    tax_to_readids_dict = json.load(open(json_tax_to_readsid_path, "r"))

    ## Check if output directory exists and create if not
    absolute_outdir = os.path.abspath(outdir)

    if not os.path.exists(absolute_outdir):
        os.makedirs(absolute_outdir)

    dump_to_files(
        sample_id,
        tax_to_readids_dict,
        fq1,fq2,absolute_outdir,
        buffer_size=buffer_size
    )

