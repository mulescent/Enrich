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
	test.count()
	test.write_enrichments("/Users/afrubin/Code/seqlib_testing/output")
