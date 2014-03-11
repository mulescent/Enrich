from variant import VariantSeqLib
from enrich_error import EnrichError
from fqread import read_fastq, check_fastq
import pandas as pd
import logging


class BasicSeqLib(VariantSeqLib):
    """
    Class for count data from sequencing libraries with a single read for 
    each variant. Creating a :py:class:`BasicSeqLib` requires a valid 
    *config* object, usually from a ``.json`` configuration file.

    Example config file for a :py:class:`BasicSeqLib`:

    .. literalinclude:: config_examples/basic.json

    :download:`Download this JSON file <config_examples/basic.json>`

    The ``"fastq"`` config entry can contain one read file, with the key 
    ``"forward"`` or ``"reverse"``. If the read file is ``"reverse"``, all 
    reads will be reverse-complemented before being compared to the wild type 
    sequence.
    """
    def __init__(self, config):
        VariantSeqLib.__init__(self, config)
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

            self.set_filters(config['filters'], {'min quality' : 0,
                                      'avg quality' : 0,
                                      'chastity' : False,
                                      'max mutations' : len(self.wt_dna)})
        except KeyError as key:
            raise EnrichError("missing required config value: %s" % key, self.name)

        try:
            check_fastq(self.reads)
        except IOError as fqerr:
            raise EnrichError("FASTQ file error: %s" % fqerr, self.name)


    def calculate(self):
        """
        Reads the forward or reverse FASTQ file (reverse reads are reverse-complemented),
        performs quality-based filtering, and counts the variants.
        """
        self.counts['variants'] = dict()

        # flags for verbose output of filtered reads
        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        for fq in read_fastq(self.reads):
            if self.revcomp_reads:
                fq.revcomp()

            for key in filter_flags:
                filter_flags[key] = False

            # filter the read based on specified quality settings
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
            if not any(filter_flags.values()): # passed quality filtering
                mutations = self.count_variant(fq.sequence)
                if mutations is None: # read has too many mutations
                    self.filter_stats['max mutations'] += 1
                    filter_flags['max mutations'] = True
            if any(filter_flags.values()):
                self.filter_stats['total'] += 1
                if self.verbose:
                    self.report_filtered_read(fq, filter_flags)

        self.counts['variants'] = \
                pd.DataFrame.from_dict(self.counts['variants'], 
                                       orient="index", dtype="int32")
        if len(self.counts['variants']) == 0:
            raise EnrichError("Failed to count variants", self.name)
        self.counts['variants'].columns = ['count']

        logging.info("Counted %d variants (%d unique) [%s]" % \
                (self.counts['variants']['count'].sum(), len(self.counts['variants'].index), self.name))
        if self.aligner is not None:
            logging.info("Aligned %d variants [%s]" % (self.aligner.calls, self.name))
        self.report_filter_stats()
