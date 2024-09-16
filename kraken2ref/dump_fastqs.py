import io
import json
import os
from Bio import SeqIO
from concurrent import futures
import gc
import sys

def dump_to_file_index(sample_id, tax_to_readids_dict, fq1, fq2, outdir, 
                        max_threads=1, buffer_size=io.DEFAULT_BUFFER_SIZE):
    """Function that dumps reads to file.

    Args:
        tax_to_readids_dict (dict): Dictionary mapping of taxon IDs to their corresponding lists of read IDs
        fq1 (str/path): Path to forward fastq file
        fq2 (str/path): Path to reverse fastq file
        outdir (str/path): Path to output directory
    """

    def fq_write_wrapper(output_taxid, sample_id, tax_to_readids_dict, fq1_dict, fq2_dict, outdir, slashes):
        """Function that wraps around actual file I/O, run in parallel

        Args:
            output_taxid (str): Taxonomic ID to write
            tax_to_readids_dict (dict): Dictionary mapping of taxon IDs to their corresponding lists of read IDs
            fq1_dict (_IndexedSeqFileDict): Dictionary of forward fastq records
            fq2_dict (_IndexedSeqFileDict): Dictionary of reverse fastq records
            outdir (str/path): Path to output directory
        """

        ## initialise files to write to
        R1_file = open(os.path.join(outdir, f"{sample_id}_{output_taxid}_R1.fq"), "w", buffering=buffer_size)
        R2_file = open(os.path.join(outdir, f"{sample_id}_{output_taxid}_R2.fq"), "w", buffering=buffer_size)

        ## iterate over read ids in list and dump to files
        for read_id in tax_to_readids_dict[output_taxid]:
            if slashes:
                try:
                    R1_file.write(fq1_dict[read_id+"/1"].format("fastq"))
                    R2_file.write(fq2_dict[read_id+"/2"].format("fastq"))

                except KeyError as ke: 
                    # if read_id not present, skip it (necessary if splitted fq are processed ) 
                    continue

            else:
                try:
                    R1_file.write(fq1_dict[read_id].format("fastq"))
                    R2_file.write(fq2_dict[read_id].format("fastq"))

                except KeyError as ke: 
                    # if read_id not present, skip it (necessary if fq is splitted) 
                    continue

    # -------------------------------------------------------------#
    ## load in fastq files as dictionaries for constant-time lookup
    fq1_dict = SeqIO.index(fq1, "fastq")
    fq2_dict = SeqIO.index(fq2, "fastq")

    ## check if read IDs have slash notations
    chosen_ref = list(fq1_dict.keys())[0]

    if chosen_ref[-2:] == "/1":
        slashes = True
    else:
        slashes = False

    ## initialise parallel function calls
    with futures.ThreadPoolExecutor(max_workers=max_threads) as executor:
        ## get list of functions for execution
        functions = [executor.submit(fq_write_wrapper(taxid, sample_id, tax_to_readids_dict, fq1_dict, fq2_dict, outdir, slashes)) for taxid in list(tax_to_readids_dict.keys())]
        ## make main wait until functions conclude
        futures.wait(functions)

def dump_to_file_chunks(sample_id, tax_to_readids_dict, fq1, fq2,
    outdir, chunk_size=10_000, buffer_size=io.DEFAULT_BUFFER_SIZE):
    """Function that dumps reads to file in chunks.

    Args:
        tax_to_readids_dict (dict): Dictionary mapping of taxon IDs to their corresponding lists of read IDs
        fq1 (str/path): Path to forward fastq file
        fq2 (str/path): Path to reverse fastq file
        outdir (str/path): Path to output directory
        chunk_size (int): Number of records to process in each batch (default is 1000)
    """

    ## Create output files for each taxid
    output_files = {}
    for taxid in tax_to_readids_dict.keys():
        R1_file = open(os.path.join(outdir, f"{sample_id}_{taxid}_R1.fq"), "w", buffering=buffer_size)
        R2_file = open(os.path.join(outdir, f"{sample_id}_{taxid}_R2.fq"), "w", buffering=buffer_size)
        output_files[taxid] = (R1_file, R2_file)

    ## Use iterators to process the fastq files
    fq1_iter = SeqIO.parse(fq1, "fastq")

    ## Check if read IDs have slash notations
    first_record = next(fq1_iter)
    slashes = first_record.id.endswith("/1")

    ## Reinitialize iterators after peeking
    fq1_iter = SeqIO.parse(fq1, "fastq")
    fq2_iter = SeqIO.parse(fq2, "fastq")

    ## Process files in chunks
    sys.stdout.write(f"@ processing chunks (chunk_size = {chunk_size})")
    c = 1
    while True:

        chunk1 = list(next(fq1_iter, None) for _ in range(chunk_size))
        chunk2 = list(next(fq2_iter, None) for _ in range(chunk_size))
        print(f"  > {c}", end="\r")
        ## Break the loop if no more records to process
        if not chunk1 or not chunk2:
            break

        ## Filter out None values in case the iterators run out at different times
        chunk1 = [record for record in chunk1 if record is not None]
        chunk2 = [record for record in chunk2 if record is not None]

        # break loop if only none values are present on chunk1
        if len(chunk1) == 0:
            break

        ## Iterate over both chunks simultaneously
        for record1, record2 in zip(chunk1, chunk2):
            read_id = record1.id.rstrip("/1") if slashes else record1.id

            ## Check which taxid (if any) the read belongs to and write to corresponding file
            for taxid, read_ids in tax_to_readids_dict.items():
                R1_file, R2_file = output_files[taxid]

                ## DEBUG ##
                # uncomment the line below to print buffer moving checks
                #print(R1_file.buffer.tell(), "buffer,", R1_file.buffer.raw.tell(), "raw")
                ###########

                if read_id in read_ids:
                    R1_file.write(record1.format("fastq"))
                    R2_file.write(record2.format("fastq"))
                    break  # No need to check further once found


        ## Explicitly clear the chunk variables and run garbage collection
        del chunk1, chunk2
        gc.collect()
        c +=1
    # print output files generated
    sys.stdout.write(":> Output files:")
    for taxid in output_files.keys():
        out_file = os.path.join(outdir, f"{sample_id}_{taxid}_R*.fq")
        sys.stdout.write(f"  -> {out_file}")


    ## Close all output files
    for R1_file, R2_file in output_files.values():
        R1_file.close()
        R2_file.close()

    gc.collect()  # Final garbage collection


def dump_fastqs(args):
    
    # inputs
    sample_id = args.sample_id
    json_tax_to_readsid_path = args.tax_to_readsid_path
    fq1 = args.fastq1
    fq2 = args.fastq2
    outdir = args.outdir
    max_threads=args.max_threads
    fq_load_mode = args.fq_load_mode
    chunk_size = args.chunk_size
    buffer_size = args.buffer_size

    # load tax to reads id dictionary
    tax_to_readids_dict = json.load(open(json_tax_to_readsid_path, "r"))

    # check if modes are valid -- /
    ACCEPTED_MODES = ["full", "chunks"]

    try:
        assert(fq_load_mode in ACCEPTED_MODES)
    except(AssertionError):
        err_msg = f"'{fq_load_mode}' is an invalid mode (valid ones: {ACCEPTED_MODES}"
        sys.stdout.write(err_msg)
    # --------------------------- #

    ## Check if output directory exists and create if not
    absolute_outdir = os.path.abspath(outdir)

    if not os.path.exists(absolute_outdir):
        os.makedirs(absolute_outdir)

    sys.stdout.write(f"> {fq_load_mode} mode selected ")
    # dump files according to mode
    if fq_load_mode == "full":
        dump_to_file_index(
            sample_id,tax_to_readids_dict,
            fq1,fq2,absolute_outdir,
            max_threads=max_threads,
            buffer_size=buffer_size
        )
    if fq_load_mode == "chunks":
        dump_to_file_chunks(
            sample_id,tax_to_readids_dict,
            fq1,fq2,absolute_outdir,
            chunk_size=chunk_size,
            buffer_size=buffer_size
        )
