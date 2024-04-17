from kraken2ref.kraken2reference import KrakenProcessor

my_processor = KrakenProcessor("test_k2r")
my_processor.analyse_report(input_kraken_report_file="tests/test_set/with_fluseg_db/testset_50x.uniq.report.txt", input_threshold=100, input_method="kmeans")
my_processor.write_output(suffix="decomposed")
