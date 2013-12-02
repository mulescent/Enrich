from __future__ import print_function
import re
import time
from sys import stdout, stderr
from collections import Counter
from itertools import izip_longest
from enrich_error import EnrichError
from aligner import Aligner

WILD_TYPE_VARIANT = "_wt"

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

def has_indel(variant):
    """
    Returns True if the variant contains an indel mutation.

    variant is a string containing mutations in HGVS format.
    """
    return any(x in variant for x in ("ins", "del", "dup"))


class SeqLib(object):
    """
    Abstract class for data from a single Enrich sequencing library.

    Implements core functionality for assessing variants, handling 
    configuration files, and other shared processes.

    Subclasess should implement the required functionality to get the  
    variant DNA sequences that are being assessed.
    """

    _filter_messages = {
            'remove unresolvable' : "unresolvable mismatch",
            'min quality' : "single-base quality",
            'avg quality' : "average quality",
            'max mutations' : "excess mutations",
			'chastity' : "not chaste",
            'remove overlap indels' : "indel in read overlap"
    }


    def __init__(self, config):
        self.name = "Unnamed" + self.__class__.__name__
        self.verbose = False
        self.log = None
        self.wt_dna = None
        self.wt_protein = None

        try:
            self.name = config['name']
            self.set_wt(config['wild type']['sequence'], 
                        coding=config['wild type']['coding'])
        except KeyError as key:
            raise EnrichError("Missing required config value '%s'" % key, 
                              self.name)

        if 'reference offset' in config['wild type']:
            try:
                self.reference_offset = int(config['wild type']
                                                  ['reference offset'])
            except ValueError:
                raise EnrichError("Invalid reference offset value", self.name)
        else:
            self.reference_offset = 0

        # initialize data
        self.counters = {'variants' : Counter(), 
                         'excluded_variants' : Counter()}
        self.filters = None
        self.aligner = Aligner()
        self.libtype = None
        self.frequencies = dict()
        self.ratios = dict()


    def enable_logging(self, log):
        self.verbose = True
        self.log = log
        try:
            print("# Logging started for '%s': %s" % 
                        (self.name, time.asctime()), file=self.log)
        except (IOError, AttributeError):
            raise EnrichException("Could not write to log file", self.name)


    def is_coding(self):
        return self.wt_protein is not None


    def count(self):
        raise NotImplementedError("must be implemented by subclass")


    def set_filters(self, config, class_default_filters):
        self.filters = class_default_filters

        for key in self.filters:
            if key in config['filters']:
                try:
                    self.filters[key] = int(config['filters'][key])
                except ValueError:
                    raise EnrichError("Invalid filter value for '%s'" % key, 
                                      self.name)

        unused = list()
        for key in config['filters']:
            if key not in self.filters:
                unused.append(key)
        if len(unused) > 0:
            raise EnrichError("Unused filter parameters (%s)" % 
                              ', '.join(unused), self.name)

        self.filter_stats = dict()
        for key in self.filters:
            self.filter_stats[key] = 0
        self.filter_stats['total'] = 0


    def report_filtered_read(self, fq, filter_flags):
        print("Filtered read (%s) [%s]" % \
                (', '.join(SeqLib._filter_messages[x] 
                 for x in filter_flags if filter_flags[x]), self.name), 
              file=self.log)
        print(fq, file=self.log)


    def set_wt(self, sequence, coding=True, codon_table=codon_table):
        """
        Set the wild type DNA sequence (and protein sequence if applicable).

        The DNA sequence must not contain ambiguity characters, but may 
        contain whitespace (which will be removed). If the sequence encodes a
        protein, it must be in-frame.
        """
        sequence = "".join(sequence.split()) # remove whitespace

        if not re.match("^[ACGTacgt]+$", sequence):
            raise EnrichError("WT DNA sequence contains unexpected "
                              "characters", self.name)

        if len(sequence) % 3 != 0 and coding:
            raise EnrichError("WT DNA sequence contains incomplete codons", 
                              self.name)
        
        self.wt_dna = sequence.upper()
        if coding:
            self.wt_protein = ""
            for i in xrange(0, len(self.wt_dna), 3):
                self.wt_protein += codon_table[self.wt_dna[i:i + 3]]
        else:
            self.wt_protein = None


    def align_variant(self, variant_dna):
        mutations = list()
        traceback = self.aligner.align(self.wt_dna, variant_dna)
        for x, y, cat, length in traceback:
            if cat == "match":
                continue
            elif cat == "mismatch":
                mut = "%s>%s" % (self.wt_dna[x], variant_dna[y])
            elif cat == "insertion":
                if y > length:
                    dup = variant_dna[y:y + length]
                    if dup == variant_dna[y - length:y]:
                        mut = "dup%s" % dup
                    else:
                        mut = "_%dins%s" % (x + 2, dup)
                else:                                    
                    mut = "_%dins%s" % (x + 2, 
                                        variant_dna[y:y + length])
            elif cat == "deletion":
                mut = "_%ddel" % (x + length)
            mutations.append((x, mut))
        return mutations


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
        if not re.match("^[ACGTNXacgtnx]+$", variant_dna):
            raise EnrichError("Variant DNA sequence contains unexpected "
                              "characters", self.name)

        variant_dna = variant_dna.upper()

        if len(variant_dna) != len(self.wt_dna):
            mutations = self.align_variant(variant_dna)
        else:
            mutations = list()
            for i in xrange(len(variant_dna)):
                if variant_dna[i] != self.wt_dna[i] and \
                        variant_dna[i] != 'N':  # ignore N's
                    mutations.append((i, "%s>%s" % \
                                       (self.wt_dna[i], variant_dna[i])))
                    if len(mutations) > self.filters['max mutations']:
                        mutations = self.align_variant(variant_dna)
                        break

        if len(mutations) > self.filters['max mutations']: # post-alignment
            # discard this variant
            return None

        mutation_strings = list()
        if self.is_coding():
            variant_protein = ""
            for i in xrange(0, len(variant_dna), 3):
                try:
                    variant_protein += codon_table[variant_dna[i:i + 3]]
                except KeyError: # garbage codon due to indel
                    variant_protein += '?'

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
            variant_string = ', '.join(mutation_strings)
        else:
            variant_string = WILD_TYPE_VARIANT
        if copies == 1:
            self.counters['variants'].update([variant_string])
        else:
            self.counters['variants'] += Counter({variant_string : copies})

        return mutation_strings


    def count_mutations(self, include_indels=False):
        """
        Count the individual mutations in all variants.

        If indels is False, all mutations in a variant that contains an 
        insertion/deletion/duplication will not be counted.

        For coding sequences, amino acid substitutions are counted
        independently of the corresponding nucleotide change.
        """
        self.counters['mutations_nt'] = Counter()
        if self.is_coding():
            self.counters['mutations_aa'] = Counter()

        for variant, count in self.counters['variants'].iteritems():
            if not has_indel(variant) or include_indels:
                muts = [x.strip() for x in variant.split(',')]
                self.counters['mutations_nt'] += \
                        Counter(dict(izip_longest(muts, [], fillvalue=count)))
                if self.is_coding():
                    muts = re.findall("p\.\w{3}\d+\w{3}", variant)
                    self.counters['mutations_aa'] += \
                            Counter(dict(izip_longest(muts, [], 
                                                      fillvalue=count)))


    def print_counter(self, key, min_count=0, indels=True, file=stdout):
        for variant, count in self.counters[key].most_common():
            if count >= min_count:
                if not has_indel(variant) or indels:
                    print("%d\t%s" % (count, variant), file=file)


    def calc_frequencies(self):
        for key in self.counters:
            self.frequencies[key] = dict()
            total = float(sum(self.counters[key].values()))
            for x in self.counters[key]:
                self.frequencies[key][x] = self.counters[key][x] / total


    def calc_ratios(self, inputlib):
        for key in self.frequencies:
            self.ratios[key] = dict()
            for x in self.frequencies[key]:
                if x in inputlib.frequencies[key]:
                    self.ratios[key][x] = self.frequencies[key][x] / \
                            inputlib.frequencies[key][x]


    def exclude_variant(self, variant):
        """
        Move a variant's data to the excluded list.

        Removes the variant's count and frequency information from the 
        'variants' data structure and adds it to the 'excluded_variants' data
        structure. Does nothing of the variant is not found in this SeqLib.
        """
        if variant in self.counters['variants']:
            self.counters['excluded_variants'] += \
                    Counter({variant : self.counters['variants'][variant]})
            del self.counters['variants'][variant]
        if variant in self.frequencies['variants']:
            self.frequencies['excluded_variants'][variant] = \
                    self.frequencies['variants'][variant]
            del self.frequencies['variants'][variant]
        if variant in self.ratios['variants']:
            self.ratios['excluded_variants'][variant] = \
                    self.ratios['variants'][variant]
            del self.ratios['variants'][variant]

