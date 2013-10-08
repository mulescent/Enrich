from seqlib import SeqLib
from enrich_error import EnrichError
from fastq_util import *


class OverlapSeqLib(SeqLib):
    def __init__(self, config):
        SeqLib.__init__(self, config)
        try:
            if check_fastq_extension(config['fastq']['forward']):
                self.forward = config['fastq']['forward']
            if check_fastq_extension(config['fastq']['reverse']):
                self.reverse = config['fastq']['reverse']
            self.fwd_start = int(config['overlap']['forward start'])
            self.rev_start = int(config['overlap']['reverse start'])
            self.overlap_length = int(config['overlap']['length'])
            self.trim = config['overlap']['overlap only']
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key)
        except ValueError as value:
            raise EnrichError("Count not convert config value: %s" % value)
        self.set_filters(config, {'remove unresolvable' : False, 
                                  'min quality' : 0,
                                  'avg quality' : 0,
                                  'max mutations' : len(self.wt_dna)})


    def fuse_reads(self, fwd, rev):
        rev = reverse_fastq(rev)
        fwd_quality = fastq_quality(fwd)
        rev_quality = fastq_quality(rev)

        rev_extra_start = self.rev_start - 1 + self.overlap_length + 1
        fused_seq = list(fwd[1] + rev[1][rev_extra_start:])
        fused_quality = fwd_quality + rev_quality[rev_extra_start:]

        consecutive_mismatches = 0
        for i in xrange(self.overlap_length):
            i_f = self.fwd_start - 1 + i
            i_r = self.rev_start - 1 + i
            if fwd[1][i_f] == rev[1][i_r]:
                consecutive_mismatches = 0
                # take the highest quality value
                if rev_quality[i_r] > fwd_quality[i_f]:
                    fused_quality[i_f] = rev_quality[i_r]
            else: # mismatch
                if fwd_quality[i_f] == rev_quality[i_r]:
                    fused_seq[i_f] = 'X' # unresolvable
                elif rev_quality[i_f] > fwd_quality[i_r]:
                    fused_seq[i_f] = rev[1][i_r]
                    fused_quality[i_f] = rev_quality[i_r]
                else:
                    pass # overlap region already same as fwd
                consecutive_mismatches += 1
                if consecutive_mismatches > 1:
                    matrix = needleman_wunsch(fwd[1], rev[1])
                    # FUSE USING THE ALIGNER
                    break

        fused = (fwd[0], fused_seq, 
                 fastq_quality_reconvert(fused_quality))
        if self.trim:
            fused = trim_fastq_length(fused, self.fwd_start, 
                                      self.overlap_length)
        return fused


    def count(self):
        # flags for verbose output of filtered reads
        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        for fwd, rev in read_fastq_multi([self.forward, self.reverse]):
            for key in filter_flags:
                filter_flags[key] = False

            # filter the read based on specified quality settings
            fused = fuse_reads(fwd, rev)
            if self.filters['unresolvable']:
                if 'X' in fused[1]:
                    self.filter_stats['unresolvable'] += 1
                    filter_flags['unresolvable'] = True
            if self.filters['min quality'] > 0:
                if fastq_min_quality(fused) < self.filters['min quality']:
                    self.filter_stats['min quality'] += 1
                    filter_flags['min quality'] = True
            if self.filters['avg quality'] > 0:
                if fastq_average_quality(fq) < self.filters['avg quality']:
                    self.filter_stats['avg quality'] += 1
                    filter_flags['avg quality'] = True
            if not any(filter_flags.values()): # passed quality filtering
                mutations = self.count_variant(fused[1])
                if mutations is None: # fused read has too many mutations
                    self.filter_stats['max mutations'] += 1
                    filter_flags['max mutations'] = True
            if any(filter_flags.values()):
                self.filter_stats['total'] += 1
                if self.verbose:
                    self.report_filtered_read(fused, filter_flags)
