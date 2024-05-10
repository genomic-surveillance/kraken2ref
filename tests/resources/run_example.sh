#!/bin/bash

set -e
set -u
set -o pipefail

kraken2ref -s test_cli parse_report \
            -i tests/test_set/with_refseq_db/testset_50x.report.txt \
            -o test_cli \
            -t 100 \
            -x decomposed \
            -q

kraken2ref -s test_cli sort_reads \
            -k tests/test_set/with_refseq_db/testset_50x.kraken.output \
            -fq1 tests/test_set/with_refseq_db/testset_50x.class_seqs_1.fq \
            -fq2 tests/test_set/with_refseq_db/testset_50x.class_seqs_2.fq \
            -r test_cli/test_cli_decomposed.json \
            -m tree -u