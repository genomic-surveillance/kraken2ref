from kraken2ref.src.taxonlevel import TaxonLevel

def test_taxonlevel():

    one = TaxonLevel("S10")
    two = TaxonLevel("S3")

    assert one > two
    assert one != two
    assert one - two == 7
    assert one + two == TaxonLevel("S13")
    assert two.to(one) == ["S4", "S5", "S6", "S7", "S8", "S9"]