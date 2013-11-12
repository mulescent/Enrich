from seqlib import SeqLib
from enrich_error import EnrichError
from fastq_util import *


class OverlapSeqLib(SeqLib):
    def __init__(self, config):
        SeqLib.__init__(self, config)
        self.libtype = "overlap"
        try:
            if check_fastq_extension(config['fastq']['forward']):
                self.forward = config['fastq']['forward']
            if check_fastq_extension(config['fastq']['reverse']):
                self.reverse = config['fastq']['reverse']
            self.fwd_start = int(config['overlap']['forward start'])
            self.rev_start = int(config['overlap']['reverse start'])
            self.overlap_length = int(config['overlap']['length'])
            self.trim = config['overlap']['overlap only']
            self.set_filters(config, {'remove unresolvable' : False, 
                                      'min quality' : 0,
                                      'avg quality' : 0,
                                      'max mutations' : len(self.wt_dna),
                                      'chastity' : False,
                                      'remove overlap indels' : True})
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key)
        except ValueError as value:
            raise EnrichError("Count not convert config value: %s" % value)


    def fuse_reads(self, fwd, rev):
        rev = reverse_fastq(rev)
        fwd_quality = fastq_quality(fwd)
        rev_quality = fastq_quality(rev)

        rev_extra_start = len(rev[1]) - self.rev_start + 1
        fwd_end = self.fwd_start + self.overlap_length - 1
        fused_seq = list(fwd[1][:fwd_end] + rev[1][rev_extra_start:])
        fused_quality = fwd_quality[:fwd_end] + rev_quality[rev_extra_start:]

        consecutive_mismatches = 0
        for i in xrange(self.overlap_length):
            a = self.fwd_start - 1 + i
            b = len(rev[1]) - self.rev_start - self.overlap_length + i + 1
            if fwd[1][a] == rev[1][b]:
                print(i, a, b, "match:", fwd[1][a], "==", rev[1][b])
                consecutive_mismatches = 0
                # take the highest quality value
                if rev_quality[b] > fwd_quality[a]:
                    fused_quality[a] = rev_quality[b]
            else: # mismatch
                if fwd_quality[a] == rev_quality[b]:
                    fused_seq[a] = 'X' # unresolvable
                elif rev_quality[b] > fwd_quality[a]:
                    fused_seq[a] = rev[1][b]
                    fused_quality[a] = rev_quality[b]
                else:
                    pass # overlap region already same as fwd
                consecutive_mismatches += 1
                if consecutive_mismatches > 1:
                    # use the aligner
                    fused_seq = list()
                    fused_quality = list()
                    traceback = self.aligner.align(fwd[1], rev[1])
                    if any(t[2] in ("insertion", "deletion") for t in 
                           traceback):
                        fused_seq = None
                    else:
                        for x, y, cat, _ in traceback:
                            if cat == "match":
                                if rev_quality[y] > rev_quality[x]:
                                    fused_quality[x] = rev_quality[y]
                            elif cat == "mismatch":
                                if fwd_quality[x] == rev_quality[y]:
                                    fused_seq[x] = 'X' # unresolvable
                                elif rev_quality[y] > fwd_quality[x]:
                                    fused_seq[x] = rev[1][y]
                                    fused_quality[x] = rev_quality[y]
                                else:
                                    pass # overlap region already same as fwd
                            else:
                                raise EnrichError("Alignment result error")
                    break

        if fused_seq is None: # fusing failed due to indels
            return None
        else:
            fused = (fwd[0], "".join(fused_seq), 
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
            if self.filters['chastity']:
                if not filter_fastq_chastity(fwd):
                    filter_flags['chastity'] = True
                    if self.verbose:
                        self.report_filtered_read(fwd, filter_flags)
                if not filter_fastq_chastity(rev):
                    filter_flags['chastity'] = True
                    if self.verbose:
                        self.report_filtered_read(rev, filter_flags)
                if filter_flags['chastity']:
                    self.filter_stats['chastity'] += 1
                    self.filter_stats['total'] += 1
                    continue
            fused = fuse_reads(fwd, rev)
            if fused is None: # fuser failed
                if self.filters['remove overlap indels']:
                    self.filter_stats['remove overlap indels'] += 1
                    self.filter_stats['total'] += 1
                    filter_flags['remove overlap indels'] = True
                    if self.verbose:
                        self.report_filtered_read(fwd, filter_flags)
                        self.report_filtered_read(rev, filter_flags)
                else:
                    raise NotImplementedError("Indel resolution for "
                            "overlapping reads is not supported")
            else:
                if self.filters['unresolvable']:
                    if 'X' in fused[1]:
                        self.filter_stats['unresolvable'] += 1
                        filter_flags['unresolvable'] = True
                if self.filters['min quality'] > 0:
                    if fastq_min_quality(fused) < self.filters['min quality']:
                        self.filter_stats['min quality'] += 1
                        filter_flags['min quality'] = True
                if self.filters['avg quality'] > 0:
                    if fastq_mean_quality(fq) < self.filters['avg quality']:
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
