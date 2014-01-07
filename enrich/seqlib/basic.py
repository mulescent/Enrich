from variant import VariantSeqLib
from enrich_error import EnrichError
from fastq_util import read_fastq, check_fastq


class BasicSeqLib(VariantSeqLib):
    """
    Class for count data from sequencing libraries with a single read for each variant.
    Creating a :py:class:`BasicSeqLib` requires a valid *config* .
    """
    def __init__(self, config):
        VariantSeqLib.__init__(self, config)
        try:
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


    def count(self):
        """
        Reads the forward or reverse FASTQ file (reverse reads are reverse-complemented),
        performs quality-based filtering, and 
        counts the variants.
        """
        self.counts['variants'] = dict()

        # flags for verbose output of filtered reads
        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        for fq in read_fastq(self.reads):
            if self.reverse_reads:
                fq.reverse()

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
                if mutations is None: # fused read has too many mutations
                    self.filter_stats['max mutations'] += 1
                    filter_flags['max mutations'] = True
            if any(filter_flags.values()):
                self.filter_stats['total'] += 1
                if self.verbose:
                    self.report_filtered_read(self.log, fq, filter_flags)

        self.counts['variants'] = \
                pd.DataFrame.from_dict(self.counts['variants'], 
                                       orient="index", dtype="int32")
        if len(self.counts['variants']) == 0:
            raise EnrichError("Failed to count variants", self.name)
        self.counts['variants'].columns = ['count']


