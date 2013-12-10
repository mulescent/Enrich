from __future__ import print_function
from selection import Selection
import sys
import json
import time

logfile_fmt = "Enrich__%Y_%m_%d__%H_%M_%S.log"

if __name__ == "__main__":
	logfile_name = time.strftime(logfile_fmt)
	config = json.load(open(sys.argv[1]))
	test = Selection(config)
	test.enable_logging(open(logfile_name, "w"))
	test.get_counts()
	test.calc_frequencies()
	test.calc_ratios()
	test.calc_enrichments()
#	test.count_mutations()
	test.write_enrichments("/Users/afrubin/Code/forAlan/LS710/output")
#	for lib in test.libraries:
#		for bc, count in lib.orphan_barcodes().iteritems():
#			print(bc, count, sep="\t")
