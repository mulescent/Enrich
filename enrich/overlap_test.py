from __future__ import print_function
import seqlib.overlap
import json
from sys import argv


config = json.load(open(argv[1], "U"))

test = seqlib.overlap.OverlapSeqLib(config['libraries'][0])
fwd = ("", "AAAAGGGCCCC", "###########")
rev = ("", "TTCCCAAAA", "##########")

print(test.fuse_reads(fwd, rev))