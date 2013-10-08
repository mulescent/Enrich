from __future__ import print_function
from sys import stdout, stderr
import os.path
import re
import json
from fastq_util import *
from collections import Counter
from itertools import izip_longest
from copy import deepcopy
from enrich_error import EnrichError


__all__ = ["BasicEnrichLibrary", "BarcodeEnrichLibrary", 
           "OverlapEnrichLibrary"]


# Standard codon table for translating wild type and variant DNA sequences
codon_table = {
        'TTT':'F', 'TCT':'S', 'TAT':'Y', 'TGT':'C',
        'TTC':'F', 'TCC':'S', 'TAC':'Y', 'TGC':'C',
        'TTA':'L', 'TCA':'S', 'TAA':'*', 'TGA':'*',
        'TTG':'L', 'TCG':'S', 'TAG':'*', 'TGG':'W',
        'CTT':'L', 'CCT':'P', 'CAT':'H', 'CGT':'R',
        'CTC':'L', 'CCC':'P', 'CAC':'H', 'CGC':'R',
        'CTA':'L', 'CCA':'P', 'CAA':'Q', 'CGA':'R',
        'CTG':'L', 'CCG':'P', 'CAG':'Q', 'CGG':'R',
        'ATT':'I', 'ACT':'T', 'AAT':'N', 'AGT':'S',
        'ATC':'I', 'ACC':'T', 'AAC':'N', 'AGC':'S',
        'ATA':'I', 'ACA':'T', 'AAA':'K', 'AGA':'R',
        'ATG':'M', 'ACG':'T', 'AAG':'K', 'AGG':'R',
        'GTT':'V', 'GCT':'A', 'GAT':'D', 'GGT':'G',
        'GTC':'V', 'GCC':'A', 'GAC':'D', 'GGC':'G',
        'GTA':'V', 'GCA':'A', 'GAA':'E', 'GGA':'G',
        'GTG':'V', 'GCG':'A', 'GAG':'E', 'GGG':'G'
    }


# Convert between single- and three-letter amino acid codes
aa_codes = {
        'Ala' : 'A', 'A' : 'Ala',
        'Arg' : 'R', 'R' : 'Arg',
        'Asn' : 'N', 'N' : 'Asn',
        'Asp' : 'D', 'D' : 'Asp',
        'Cys' : 'C', 'C' : 'Cys',
        'Glu' : 'E', 'E' : 'Glu',
        'Gln' : 'Q', 'Q' : 'Gln',
        'Gly' : 'G', 'G' : 'Gly',
        'His' : 'H', 'H' : 'His',
        'Ile' : 'I', 'I' : 'Ile',
        'Leu' : 'L', 'L' : 'Leu',
        'Lys' : 'K', 'K' : 'Lys',
        'Met' : 'M', 'M' : 'Met',
        'Phe' : 'F', 'F' : 'Phe',
        'Pro' : 'P', 'P' : 'Pro',
        'Ser' : 'S', 'S' : 'Ser',
        'Thr' : 'T', 'T' : 'Thr',
        'Trp' : 'W', 'W' : 'Trp',
        'Tyr' : 'Y', 'Y' : 'Tyr',
        'Val' : 'V', 'V' : 'Val',
        'Ter' : '*', '*' : 'Ter'
}


# Information strings for the different kinds of read filtering
# Used for properly formatting error messages
filtered_read_messages = {
        'remove unresolvable' : "unresolvable mismatch",
        'min quality' : "single-base quality",
        'avg quality' : "average quality",
        'mutations' : "excess mutations"
}



def has_indel(variant):
    """
    Returns True if the variant contains an indel mutation.

    variant is a string containing mutations in HGVS format.
    """
    return any(x in variant for x in ("ins", "del", "dup"))


def needleman_wunsch(a, b):
    """
    Align the two sequences using the Needleman-Wunsch algorithm.

    Returns the alignment matrix.
    """
    return m



class EnrichLibrary(object):
    """
    Abstract class for data from a single Enrich sequencing library.

    Implements core functionality for assessing variants, handling 
    configuration files, and other shared processes.

    Subclasess should implement the required functionality to get the  
    variant DNA sequences that are being assessed.
    """
    def __init__(self, config):
        # set values from the config object
        try:
            self.set_wt(config['wild type']['sequence'], 
                        coding=config['wild type']['coding'])
            self.reference_offset = config['wild type']['reference offset']
        except KeyError as key:
            raise EnrichError("Missing required config value '%s'" % key)

        # initialize values to be set by EnrichExperiment
        self.verbose = False
        self.log = None

        # initialize data
        self.variants = Counter()
        self.mutations_nt = Counter()
        self.mutations_aa = Counter()
        self.filters = None


    def is_coding(self):
        return self.wt_protein is not None


    def count(self):
        raise NotImplementedError("must be implemented by subclass")


    def set_filters(self, config, class_default_filters):
        self.filters = class_default_filters

        for key in self.filters:
            if key in config['filters']:
                self.filters[key] = config['filters'][key]

        unused = list()
        for key in config['filters']:
            if key not in self.filters:
                unused.append(key)
        if len(unused) > 0:
            print("Warning: unused filter parameters (%s)" % \
                  ', '.join(unused), file=stderr)

        self.filter_stats = dict()
        for key in self.filters:
            self.filter_stats[key] = 0


    def report_filtered_read(self, fq, filter_flags):
        print("Filtered read (%s)" % (', '.join(filtered_read_messages(x)
                    for x in filter_flags if filter_flags[x])),
              file=self.log)
        print_fastq(fq, file=self.log)


    def set_wt(self, sequence, coding=True, codon_table=codon_table):
        """
        Set the wild type DNA sequence (and protein sequence if applicable).

        The DNA sequence must not contain ambiguity characters, but may 
        contain whitespace (which will be removed). If the sequence encodes a
        protein, it must be in-frame.
        """
        self.wt_dna = None
        self.wt_protein = None
        sequence = "".join(sequence.split()) # remove whitespace

        if not re.match("^[ACGTacgt]+$", sequence):
            raise EnrichError("WT DNA sequence contains unexpected characters")

        if len(sequence) % 3 != 0 and coding:
            raise EnrichError("WT DNA sequence contains incomplete codons")
        
        self.wt_dna = sequence.upper()
        if coding:
            self.wt_protein = ""
            for i in xrange(0, len(self.wt_dna), 3):
                self.wt_protein += codon_table[self.wt_dna[i:i + 3]]


    def count_variant(self, variant_dna, copies=1, include_indels=True, 
                      codon_table=codon_table):
        """
        Identify and count the variant.

        The algorithm attempts to call variants by comparing base-by-base.
        If the variant and wild type DNA are different lengths, or if there
        are an excess of mismatches (indicating a possible indel), local
        alignment is performed.

        Each variant is stored as a tab-delimited string of mutations in HGVS 
        format. Returns a list of HGSV variant strings. Returns an empty list 
        if the variant is wild type. Returns None if the variant was discarded
        due to excess mismatches.
        """
        if len(variant_dna) != len(self.wt_dna):
            mutations = self.align_variant(variant_dna)
        else:
            mutations = list()
            for i in xrange(len(variant_dna)):
                if variant_dna[i] != self.wt_dna[i] and 
                        variant_dna[i].upper() != 'N':  # ignore N's
                    mutations.append((i, "%s>%s" % \
                                       (self.wt_dna[i], variant_dna[i])))
                    if len(mutations) > self.filters['max mutations']:
                        matrix = needleman_wunsch(variant_dna, self.wt_dna)
                        mutations = list()
                        # PROCESS THE MATRIX
                        insertion = "_%dins%s" % (pos + 1, nuc)
                        deletion = "_%ddel" % (pos)
                        duplication = "dup%s" % (nuc)
                        break

        if len(mutations) > self.filters['max mutations']: # post-alignment
            # discard this variant
            return None

        mutation_strings = list()
        if self.is_coding():
            variant_protein = ""
            for i in xrange(0, len(variant_dna), 3):
                variant_protein += codon_table[variant_dna[i:i + 3]]

            for pos, change in mutations:
                ref_dna_pos = pos + self.reference_offset + 1
                ref_pro_pos = (pos + self.reference_offset) / 3 + 1
                mut = "c.%d%s" % (ref_dna_pos, change)
                if has_indel(change):
                    mut += " (p.%s%dfs)" % \
                            (aa_codes[self.wt_protein[pos / 3]], ref_pro_pos)
                elif variant_protein[pos / 3] == self.wt_protein[pos / 3]:
                    mut += " (p.=)"
                else:
                    mut += " (p.%s%d%s)" % \
                            (aa_codes[self.wt_protein[pos / 3]], ref_pro_pos,
                             aa_codes[variant_protein[pos / 3]])
                mutation_strings.append(mut)
        else:
            for pos, change in mutations:
                ref_dna_pos = pos + self.reference_offset + 1
                mut = "n.%d%s" % (ref_dna_pos, change)
                mutation_strings.append(mut)

        if len(mutation_strings) > 0:
            variant_string = '\t'.join(mutation_strings)
            if copies == 1:
                self.variants.update([variant_string])
            else:
                self.variants += Counter({variant_string : copies})

        return mutation_strings


    def count_mutations(self, include_indels=False):
        """
        Count the individual mutations in all variants.

        If include_indels is false, all mutations in a variant that contains
        an insertion/deletion/duplication will not be counted.

        For coding sequences, amino acid substitutions are counted
        independently of the corresponding nucleotide change.
        """
        self.mutations_nt.clear()
        self.mutations_aa.clear()

        for variant, count in self.variants.iteritems():
            if not has_indel(variant) or include_indels:
                m = variant.split('\t')
                self.mutations_nt += Counter(dict(
                        izip_longest(m, [], fillvalue=count)))
                if self.is_coding():
                    m = re.findall("p\.\w{3}\d+\w{3}", variant)
                    self.mutations_aa += Counter(dict(
                            izip_longest(m, [], fillvalue=count)))



class BarcodeEnrichLibrary(EnrichLibrary):
    def __init__(self):
        EnrichLibrary.__init__(self)


    def count(self):
        pass


    def read_barcode_map(self, fname):
        pass


    def filter_barcodes(self):
        pass



class OverlapEnrichLibrary(EnrichLibrary):
    def __init__(self, config):
        EnrichLibrary.__init__(self, config)
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



class BasicEnrichLibrary(EnrichLibrary):
    def __init__(self, config):
        EnrichLibrary.__init__(self, config)
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



