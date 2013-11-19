from seqlib import SeqLib
from enrich_error import EnrichError
from fastq_util import *


class PureBarcodeSeqLib(SeqLib):
	def __init__(self, config):
		SeqLib.__init__(self, config)
	
	