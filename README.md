# kraken2ref
This tool takes the output files from a kraken2 run and performs the task of extracting reads per reference sequence. This is similar to the functionality of the publicly available krakentools script [extract_kraken_reads.py](https://github.com/jenniferlu717/KrakenTools/blob/master/extract_kraken_reads.py) but our tool implements our specific decision logic for selecting reference genomes from the kraken output, which requires analyses of the branches of the taxonomy tree and uses the number of hits (reads) assigned directly at the various levels of nodes ( such as "species level 1 (S1)", "species level 2 (S2)", etc.).

The public tool is simply a python script that traverses the kraken output file, the kraken taxonomy report file and the original read input files in pure python, which is essentially what we recreate here but with the added decision logic of reference genome selection. While it would be technically possible to add this decision logic into the existing utility, it is better for testing purposes to re-implement the functionality from the ground up.

## Installation
Install with pip. You will probably want to create a venv for this first.
```shell
pip install kraken2ref@git+ssh://git@gitlab.internal.sanger.ac.uk/malariagen1/misc_utils/kraken2ref.git
```

## Usage
Once installed, run with the following command:

TODO

```shell
kraken2ref 
```

Get help:
```shell
kraken2ref -h
```

## Details about the algorithm
