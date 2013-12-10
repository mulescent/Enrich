from __future__ import print_function
import re
from seqlib import SeqLib
from enrich_error import EnrichError
from fastq_util import read_fastq, check_fastq
import pandas as pd

# debugging
from sys import stdout, stderr


class BarcodeSeqLib(SeqLib):
    def __init__(self, config, parent=True):
        if parent:
            SeqLib.__init__(self, config)
        try:
            if 'map file' in config['barcodes']:
                self.barcode_map = BarcodeMap(config['barcodes']['map file'])
            else:
                self.barcode_map = None

            if 'forward' in config['fastq'] and 'reverse' in config['fastq']:
                raise EnrichError("Multiple FASTQ files specified", self.name)
            elif 'forward' in config['fastq']:
                self.reads = config['fastq']['forward']
                self.reverse_reads = False
            elif 'reverse' in config['fastq']:
                self.reads = config['fastq']['reverse']
                self.reverse_reads = True
            else:
                raise KeyError("'forward' or 'reverse'")

            if 'start' in config['fastq']:
                self.bc_start = config['fastq']['start']
            else:
                self.bc_start = 1
            if 'length' in config['fastq']:
                self.bc_length = config['fastq']['length']
            else:
                self.bc_length = 2147483647 # longer than any read... for now

            self.set_filters(config, {'min quality' : 0,
                                      'avg quality' : 0,
                                      'chastity' : False,
                                      'max mutations' : len(self.wt_dna)})
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, self.name)

        try:
            check_fastq(self.reads)
        except IOError as fqerr:
            raise EnrichError("FASTQ file error: %s" % fqerr, self.name)

        self.counts['barcodes'] = None


    def count(self):
        self.counts['barcodes'] = dict()

        # flags for verbose output of filtered reads
        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        # count all the barcodes
        for fq in read_fastq(self.reads):
            fq.trim_length(self.bc_length, start=self.bc_start)
            if self.reverse_reads:
                fq.reverse()

            for key in filter_flags:
                filter_flags[key] = False

            # filter the barcode based on specified quality settings
            if self.filters['chastity']:
                if not fq.is_chaste():
                    self.filter_stats['chastity'] += 1
                    filter_flags['chastity'] = True
            if self.filters['min quality'] > 0:
                if fq.min_quality() < self.filters['min quality']:
                    self.filter_stats['min quality'] += 1
                    filter_flags['min quality'] = True
            if self.filters['avg quality'] > 0:
                if fq.mean_quality() < self.filters['avg quality']:
                    self.filter_stats['avg quality'] += 1
                    filter_flags['avg quality'] = True
            if any(filter_flags.values()): # failed quality filtering
                self.filter_stats['total'] += 1
                if self.verbose:
                    self.report_filtered_read(fq, filter_flags)
            else: # passed quality filtering
                try:
                    self.counts['barcodes'][fq.sequence.upper()] += 1
                except KeyError:
                    self.counts['barcodes'][fq.sequence.upper()] = 1

        self.counts['barcodes'] = \
                pd.DataFrame.from_dict(self.counts['barcodes'], 
                                       orient="index", dtype="int32")
        if len(self.counts['barcodes']) == 0:
            raise EnrichError("Failed to count barcodes", self.name)
        self.counts['barcodes'].columns = ['count']

