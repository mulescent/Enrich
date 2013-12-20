from __future__ import print_function
from sys import stderr



from variant import VariantSeqLib
from enrich_error import EnrichError
from fastq_util import read_fastq_multi, check_fastq, FQRead
import pandas as pd


class OverlapSeqLib(VariantSeqLib):
    def __init__(self, config):
        VariantSeqLib.__init__(self, config)
        try:
            self.forward = config['fastq']['forward']
            self.reverse = config['fastq']['reverse']

            self.fwd_start = int(config['overlap']['forward start'])
            self.rev_start = int(config['overlap']['reverse start'])
            self.overlap_length = int(config['overlap']['length'])
            self.trim = config['overlap']['overlap only']
            self.max_overlap_mismatches = int(config['overlap']
                                                    ['max mismatches'])

            if 'fuser failure' in config['filters']:
                raise EnrichError("'fuser failure' is not user-configurable", self.name)
            self.set_filters(config, {'remove unresolvable' : False, 
                                      'min quality' : 0,
                                      'avg quality' : 0,
                                      'max mutations' : len(self.wt_dna),
                                      'chastity' : False,
                                      'fuser failure' : True})
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, 
                              self.name)
        except ValueError as value:
            raise EnrichError("Invalid parameter value %s" % value, self.name)

        try:
            check_fastq(self.forward)
            check_fastq(self.reverse)
        except IOError as fqerr:
            raise EnrichError("FASTQ file error: %s" % fqerr, self.name)


    def fuse_reads(self, fwd, rev):
        rev.reverse()

        rev_extra_start = len(rev) - self.rev_start + 1
        fwd_end = self.fwd_start + self.overlap_length - 1
        fused = FQRead(header=fwd.header, 
                       sequence="A",
                       header2=fwd.header2,
                       quality="#",
                       qbase=fwd.qbase)
        fused.sequence = list(fwd.sequence[:fwd_end] + \
                                     rev.sequence[rev_extra_start:])
        fused.quality = fwd.quality[:fwd_end] + \
                                     rev.quality[rev_extra_start:]

        mismatches = 0
        for i in xrange(self.overlap_length):
            a = self.fwd_start - 1 + i
            b = len(rev) - self.rev_start - self.overlap_length + i + 1
            if fwd.sequence[a] == rev.sequence[b]:
                # take the highest quality value
                if rev.quality[b] > fwd.quality[a]:
                    fused.quality[a] = rev.quality[b]
            else:
                mismatches += 1
                if fwd.quality[a] == rev.quality[b]:
                    fused.sequence[a] = 'X' # unresolvable
                elif rev.quality[b] > fwd.quality[a]:
                    fused.sequence[a] = rev.sequence[b]
                    fused.quality[a] = rev.quality[b]
                else:
                    pass # overlap region already same as fwd
        if mismatches > self.max_overlap_mismatches:
            return None # fusing failed

        fused.sequence = "".join(fused.sequence)
        if self.trim:
            fused.trim_length(self.fwd_start, self.overlap_length)
        return fused


    def count(self):
        self.counts['variants'] = dict()

        # flags for verbose output of filtered reads
        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        for fwd, rev in read_fastq_multi([self.forward, self.reverse]):
            for key in filter_flags:
                filter_flags[key] = False

            # filter the read based on specified quality settings
            if self.filters['chastity']:
                if not fwd.is_chaste():
                    filter_flags['chastity'] = True
                    if self.verbose:
                        self.report_filtered_read(fwd, filter_flags)
                if not rev.is_chaste():
                    filter_flags['chastity'] = True
                    if self.verbose:
                        self.report_filtered_read(rev, filter_flags)
                if filter_flags['chastity']:
                    self.filter_stats['chastity'] += 1
                    self.filter_stats['total'] += 1
                    continue
            fused = self.fuse_reads(fwd, rev)
            if fused is None: # fuser failed
                self.filter_stats['fuser failure'] += 1
                self.filter_stats['total'] += 1
                filter_flags['fuser failure'] = True
                if self.verbose:
                    self.report_filtered_read(fwd, filter_flags)
                    self.report_filtered_read(rev, filter_flags)
            else:
                if self.filters['remove unresolvable']:
                    if 'X' in fused.sequence:
                        self.filter_stats['remove unresolvable'] += 1
                        filter_flags['remove unresolvable'] = True
                if self.filters['min quality'] > 0:
                    if fused.min_quality() < self.filters['min quality']:
                        self.filter_stats['min quality'] += 1
                        filter_flags['min quality'] = True
                if self.filters['avg quality'] > 0:
                    if fused.mean_quality() < self.filters['avg quality']:
                        self.filter_stats['avg quality'] += 1
                        filter_flags['avg quality'] = True
                if not any(filter_flags.values()): # passed quality filtering
                    mutations = self.count_variant(fused.sequence)
                    if mutations is None: # fused read has too many mutations
                        self.filter_stats['max mutations'] += 1
                        filter_flags['max mutations'] = True
                if any(filter_flags.values()):
                    self.filter_stats['total'] += 1
                    if self.verbose:
                        self.report_filtered_read(fused, filter_flags)

        self.counts['variants'] = \
                pd.DataFrame.from_dict(self.counts['variants'], 
                                       orient="index", dtype="int32")
        if len(self.counts['variants']) == 0:
            raise EnrichError("Failed to count variants", self.name)
        self.counts['variants'].columns = ['count']


