from __future__ import print_function
from sys import stdout, stderr
import os.path
import re
import json
from fastq_util import *
from collections import Counter
from itertools import izip_longest
from copy import deepcopy


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


def has_indel(variant):
    """
    Returns True if the variant contains an indel mutation.

    variant is a string containing mutations in HGVS format.
    """
    return any(x in variant for x in ("ins", "del", "dup"))


def read_config_file(config_file):
    ext = os.path.splitext(config_file)[-1].lower()
    if ext != ".json":
        if len(ext) > 0:
            print("Warning: unrecognized config file extension '%s'" \
                  % ext, file=stderr)
        else:
            print("Warning: config file has no extension", file=stderr)

    with open(config_file, "U") as handle:
        config = json.load(handle)
    return config


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
    def __init__(self):
        self.config = None
        self.wt_dna = None
        self.wt_protein = None
        self.variants = Counter()
        self.mutations_nt = Counter()
        self.mutations_aa = Counter()
        self.max_mutations = 2
        self.verbose_filter = None


    def is_coding(self):
        return self.wt_protein is not None


    def read_config(self, config):
        """
        Set instance variables based on the configuration object.

        This method should be extended by derived classes to set class-
        specific instance variables.
        """
        self.config = config
        # set instance variables


    def parse_fastq(self):
        raise NotImplementedError("Method must be implemented by subclass")


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
            print("Error: WT DNA sequence contains unexpected characters",
                  file=stderr)
            return

        if len(sequence) % 3 != 0 and coding:
            print("Error: WT DNA sequence contains incomplete codons",
                  file=stderr)
            return
        
        self.wt_dna = sequence.upper()
        if coding:
            self.wt_protein = ""
            for i in xrange(0, len(self.wt_dna), 3):
                self.wt_protein += codon_table[self.wt_dna[i:i + 3]]

        return


    def count_variant(self, variant_dna, copies=1, reference_offset=0,
                      include_indels=True, codon_table=codon_table):
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
                    if len(mutations) > self.max_variants: # abort and align
                        matrix = needleman_wunsch(variant_dna, self.wt_dna)
                        mutations = list()
                        # PROCESS THE MATRIX
                        insertion = "_%dins%s" % (pos + 1, nuc)
                        deletion = "_%ddel" % (pos)
                        duplication = "dup%s" % (nuc)
                        break

        if len(mutations) > self.max_variants: # post-alignment
            # discard this variant
            return None

        mutation_strings = list()
        if self.is_coding():
            variant_protein = ""
            for i in xrange(0, len(variant_dna), 3):
                variant_protein += codon_table[variant_dna[i:i + 3]]

            for pos, change in mutations:
                ref_dna_pos = pos + reference_offset + 1
                ref_pro_pos = (pos + reference_offset) / 3 + 1
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
                ref_dna_pos = pos + reference_offset + 1
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


    def read_config(self, config_file):
        EnrichLibrary.read_config(self, config_file)


    def build_filter_function(self):
        pass


    def read_barcode_map(self, fname):
        pass


    def filter_barcodes(self):
        pass



class OverlapEnrichLibrary(EnrichLibrary):
    def __init__(self):
        EnrichLibrary.__init__(self)
        self.overlap = self.config['overlap']
        self.filters = None
        self.forward_fname = None
        self.reverse_fname = None


    def read_config(self, config_file):
        EnrichLibrary.read_config(self, config_file)
        self.set_filters()


    def set_filters(self):
        self.filters = {'unresolvable' : False, 'min quality' : 0,
                        'avg quality' : 0}

        for key in self.filters:
            if key in self.config['filters']:
                self.filters[key] = self.config['filters'][key]


    def fuse_reads(self, fwd, rev):
        rev = reverse_fastq(rev)
        fwd_quality = fastq_quality_str2list(fwd)
        rev_quality = fastq_quality_str2list(rev)

        rev_extra_start = self.overlap['reverse start'] - 1 + \
                          self.overlap['length'] + 1
        fused_seq = list(fwd[1] + rev[1][rev_extra_start:])
        fused_quality = fwd_quality + rev_quality[rev_extra_start:]

        consecutive_mismatches = 0
        for i in xrange(self.overlap['length']):
            i_f = self.overlap['foward start'] - 1 + i
            i_r = self.overlap['reverse start'] - 1 + i
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
        if self.overlap['overlap only']:
            fused = trim_fastq_length(fused, self.overlap['forward start'], 
                                      self.overlap['length'])
        return fused


    def GOOD_FUNCTION_NAME(self, forward_file, reverse_file):
        # initialize filtered read counters
        self.filters['stats'] = {'unresolvable' : 0, 'min quality' : 0,
                                 'avg quality' : 0, 'total' : 0}

        excess_mismatches = 0
        if self.verbose_filter is not None:
            excess_mismatches_list = list()

        fq_iterator = read_fastq_multi([forward_file, reverse_file])
        for fwd, rev in fq_iterator:
            # initialize flags for verbose output of filtered reads
            unresolvable_flag = False
            below_min_quality_flag = False
            below_avg_quality_flag = False
            filtered_flag = False

            # filter the read based on specified quality settings
            fused = fuse_reads(fwd, rev)
            if self.filters['unresolvable']:
                if 'X' in fused[1]:
                    self.filters['stats']['unresolvable'] += 1
                    unresolvable_flag = True
            if self.filters['min quality'] > 0:
                if fastq_min_quality(fused) < self.filters['min quality']:
                    self.filters['stats']['min quality'] += 1
                    below_min_quality_flag = True
            if self.filters['avg quality'] > 0:
                if fastq_average_quality(fq) < self.filters['avg quality']:
                    self.filters['stats']['avg quality'] += 1
                    below_avg_quality_flag = True
            if any((unresolvable_flag, below_min_quality_flag, 
                    below_avg_quality_flag)):
                self.filters['stats']['total'] += 1
                if self.verbose_filter is not None:
                    filter_messages = list()
                    if unresolvable_flag:
                        filter_messages.append("unresolvable mismatch")
                    if below_min_quality_flag:
                        filter_messages.append("single-base quality")
                    if below_avg_quality_flag:
                        filter_messages.append("average quality")
                    print("Filtered read (%s)" % (', '.join(filter_messages)),
                          file=self.verbose_filter)
                    print_fastq(fused, file=self.verbose_filter)
            else:
                mutations = self.count_variant(fused[1])
                if mutations is None: # fused read was discarded
                    excess_mismatches += 1
                    if self.verbose_filter is not None:
                        excess_mismatches_list.append(fused[1])

        if self.verbose_filter is not None:
            print("Variants with excess mismatches", file=self.verbose_filter)
            for x in excess_mismatches_list:
                print(x, file=self.verbose_filter)



class BasicEnrichLibrary(EnrichLibrary):
    def __init__(self):
        EnrichLibrary.__init__(self)


    def read_config(self, config_file):
        EnrichLibrary.read_config(self, config_file)


    def build_filter_function(self):
        pass


