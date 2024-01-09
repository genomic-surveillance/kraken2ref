import pytest
import os.path
from importlib_resources import files

from kraken2ref.src.kraken_taxonomy_report import KrakenTaxonomyReport

FIXTURE_DIR = files('kraken2ref.tests.fixtures')
NCBI_TEST_FILE = FIXTURE_DIR.joinpath(os.path.join('all_ncbi_flu_download','ncbi_flu.fna'))

def test_init():
    ktr = KrakenTaxonomyReport( in_file= FIXTURE_DIR.joinpath(os.path.join('report.txt')))
    assert ktr
    