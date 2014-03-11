from __future__ import print_function
import re
from seqlib import SeqLib
from enrich_error import EnrichError
from fqread import read_fastq, check_fastq
import pandas as pd

# debugging
from sys import stdout, stderr


class BarcodeSeqLib(SeqLib):
    """
    Class for count data from barcoded sequencing libraries. Designed for 
    barcode-only quantification or as a parent class for 
    :py:class:`~seqlib.barcodevariant.BarcodeVariantSeqLib`. Creating a 
    :py:class:`~seqlib.barcode.BarcodeSeqLib` requires a valid *config* 
    object with a ``"barcodes"`` entry (this entry can be empty).

    Example config file for a :py:class:`~seqlib.barcode.BarcodeSeqLib`:

    .. literalinclude:: config_examples/barcode.json

    :download:`Download this JSON file <config_examples/barcode.json>`

    The ``"fastq"`` config entry can contain one read file, with the key 
    ``"forward"`` or ``"reverse"``. If the read file is ``"reverse"``, all 
    barcodes will be reverse-complemented before being counted. The 
    ``"fastq"`` entry can contain optional values ``"start"`` and 
    ``"length"``, which will be used to trim the barcodes before counting. 
    Bases are counted starting at 1.

    The ``"min count"`` entry in ``"barcodes"`` is used to filter out 
    low-abundance barcodes (likely the result of technical artifacts). 
    Barcodes that occur less than ``"min count"`` times are removed from the 
    dataset. Setting the ``"min count"`` option appropriately can 
    dramatically improve execution time and reduce memory usage.

    .. note:: The :py:class:`~seqlib.barcodevariant.BarcodeVariantSeqLib` \
    class implements an alternative method for removing artifactual barcodes \
    that may be more appropriate for users of that module.
    """
    def __init__(self, config, barcodevariant=False):
        self.barcodevariant = barcodevariant
        if not self.barcodevariant:
            SeqLib.__init__(self, config)
        try:
            if 'forward' in config['fastq'] and 'reverse' in config['fastq']:
                raise EnrichError("Multiple FASTQ files specified", self.name)
            elif 'forward' in config['fastq']:
                self.reads = config['fastq']['forward']
                self.revcomp_reads = False
            elif 'reverse' in config['fastq']:
                self.reads = config['fastq']['reverse']
                self.revcomp_reads = True
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

            if 'min count' in config['barcodes']:
                self.min_count = config['barcodes']['min count']
            else:
                self.min_count = 0

            self.set_filters(config['filters'], {'min quality' : 0,
                                      'avg quality' : 0,
                                      'chastity' : False})
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, self.name)

        try:
            check_fastq(self.reads)
        except IOError as fqerr:
            raise EnrichError("FASTQ file error: %s" % fqerr, self.name)

        self.counts['barcodes'] = None


    def calculate(self):
        """
        Reads the forward or reverse FASTQ file (reverse reads are 
        reverse-complemented), performs quality-based filtering, and counts 
        the barcodes.
        """
        self.counts['barcodes'] = dict()

        # flags for verbose output of filtered reads
        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        # count all the barcodes
        for fq in read_fastq(self.reads):
            fq.trim_length(self.bc_length, start=self.bc_start)
            if self.revcomp_reads:
                fq.revcomp()

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
        self.counts['barcodes'] = \
                self.counts['barcodes'][self.counts['barcodes']['count'] \
                    > self.min_count]

        logging.info("Counted %d barcodes (%d unique) [%s]" % \
                (self.counts['barcodes']['count'].sum(), len(self.counts['barcodes'].index), self.name))
        if not self.barcodevariant:
            self.report_filter_stats()
