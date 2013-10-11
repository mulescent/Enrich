from __future__ import print_function
import seqlib.basic
import seqlib.barcode
import seqlib.overlap
import json
from sys import argv


def choose_seqlib(config):
	if 'barcodes' in config:
		return seqlib.barcode.BarcodeSeqLib(config)
	elif 'overlap' in config:
		return seqlib.overlap.OverlapSeqLib(config)
	else:
		return seqlib.basic.BasicSeqLib(config)


config = json.load(open(argv[1], "U"))

for i, lib in enumerate(config['libraries']):
	seq = choose_seqlib(lib)
	seq.enable_logging(open("seqlib_log_%d.txt" % (i + 1), "w"))
	seq.count()
	seq.count_mutations()
	seq.print_variants(file=open("seqlib_variants_%d.txt" % (i + 1), "w"))
	seq.print_mutations(file=open("seqlib_mutations_%d.txt" % (i + 1), "w"))

