from seqlib import SeqLib
from enrich_error import EnrichError
from fastq_util import *


class BasicSeqLib(SeqLib):
    def __init__(self, config):
        SeqLib.__init__(self, config)
        try:
            if 'forward' in config['fastq'] and 'reverse' in config['fastq']:
                raise EnrichError("Multiple FASTQ files specified")
            elif 'forward' in config['fastq']:
                if check_fastq_extension(config['fastq']['forward']):
                    self.reads = config['fastq']['forward']
                    self.reverse_reads = False
            elif 'reverse' in config['fastq']:
                if check_fastq_extension(config['fastq']['reverse']):
                    self.reads = config['fastq']['reverse']
                    self.reverse_reads = True
            else:
                raise KeyError("'forward' or 'reverse'")
        except KeyError as key:
            raise EnrichError("Missing required config value: %s" % key)
        except ValueError as value:
            raise EnrichError("Count not convert config value: %s" % value)
        self.set_filters(config, {'min quality' : 0,
                                  'avg quality' : 0,
                                  'max mutations' : len(self.wt_dna)})


    def count(self):
        # flags for verbose output of filtered reads
        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        for fq in read_fastq(self.reads):
            if self.reverse_reads:
                fq = reverse_fastq(fq)

            for key in filter_flags:
                filter_flags[key] = False

            # filter the read based on specified quality settings
            if self.filters['min quality'] > 0:
                if fastq_min_quality(fq) < self.filters['min quality']:
                    self.filter_stats['min quality'] += 1
                    filter_flags['min quality'] = True
            if self.filters['avg quality'] > 0:
                if fastq_average_quality(fq) < self.filters['avg quality']:
                    self.filter_stats['avg quality'] += 1
                    filter_flags['avg quality'] = True
            if not any(filter_flags.values()): # passed quality filtering
                mutations = self.count_variant(fq[1])
                if mutations is None: # fused read has too many mutations
                    self.filter_stats['max mutations'] += 1
                    filter_flags['max mutations'] = True
            if any(filter_flags.values()):
                self.filter_stats['total'] += 1
                if self.verbose:
                    self.report_filtered_read(fq, filter_flags)



