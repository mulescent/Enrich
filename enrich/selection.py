from __future__ import print_function
from enrich_error import EnrichError
from scipy import stats
from seqlib.basic import BasicSeqLib
from seqlib.barcode import BarcodeSeqLib, BarcodeMap
from seqlib.overlap import OverlapSeqLib
import os
import math
import itertools
import time
from collections import namedtuple


EnrichmentStats = namedtuple("EnrichmentStats", "score, r_sq")


class Selection(object):
    def __init__(self, config):
        self.name = "Unnamed" + self.__class__.__name__
        self.libraries = list()
        self.verbose = False
        self.log = None

        try:
            self.name = config['name']
            if 'barcodes' in config:
                if 'map file' in config['barcodes']:
                    self.barcode_map = BarcodeMap(config['barcodes']
                                                        ['map file'])
                else:
                    self.barcode_map = None

            for lib in config['libraries']:
                if 'barcodes' in lib:
                    new = BarcodeSeqLib(lib, barcode_map=self.barcode_map)
                elif 'overlap' in lib:
                    new = OverlapSeqLib(lib)
                else:
                    new = BasicSeqLib(lib)
                self.libraries.append(new)

        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, self.name)

        if len(self.libraries) == 0:
            raise EnrichError("Selection class has no libraries", self.name)
        self.check_wt()

        try:
            if 'correction' in config:
                if config['correction']['method'] == "stop":
                    if not self.libraries[0].is_coding():
                        raise EnrichError("Invalid correction method for "
                                          "noncoding sequences", self.name)
                    else:
                        config['correction']['length percentile'] # must exist
                        self.correction = config['correction']
            else:
                self.correction = None
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, self.name)

        self.enrichments = dict()


    def enable_logging(self, log):
        self.verbose = True
        self.log = log
        try:
            print("# Logging started for '%s': %s" % 
                        (self.name, time.asctime()), file=self.log)
        except (IOError, AttributeError):
            raise EnrichException("Could not write to log file", self.name)
        for lib in self.libraries:
            lib.enable_logging(log)


    def check_wt(self):
        try:
            coding = [lib.is_coding() for lib in self.libraries]
            sequence = [lib.wt_dna for lib in self.libraries]
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, self.name)
        
        if len(set(coding)) != 1 or len(set(sequence)) != 1:
            raise EnrichError("Inconsistent wild type sequences", self.name)


    def count(self):
        for lib in self.libraries:
            lib.count()
            lib.calc_frequencies()
            lib.calc_ratios(self.libraries[0])

        for key in self.libraries[0].frequencies:
            if all(key in lib.frequencies for lib in self.libraries):
                self.enrichments[key] = dict()

        self.correct_frequencies()
        for key in self.enrichments:
            self.calc_enrichments(key)
        if 'barcodes' in self.enrichments:
            self.filter_barcodes()


    def count_mutations(self):
        # needs to happen all the filtering/exclusion of variants
        for lib in self.libraries:
            lib.count_mutations()

        for key in ('mutations_nt', 'mutations_aa'):
            if all(key in lib.frequencies for lib in self.libraries):
                self.enrichments[key] = dict()
                self.calc_enrichments(key)


    def correct_frequencies(self):
        if self.correction is not None:
            if self.correction['method'] == "stop":
                maxpos = len(self.libraries[0].wt_protein) * \
                        self.correction['length percentile'] / 100.0
                for lib in self.libraries:
                    stops = 0
                    for variant in lib.counters['variants']:
                        m = variant.split('\t')
                        # is it a stop?
                        # is it early enough in the sequence?
                        # if so, count it
                    lib.correction = 1 - stops / float(sum(lib.counters['cariants'].values()))
                

            # correct things
            pass
        else:
            pass


    def write_enrichments(self, directory):
        enrichment_dir = os.path.join(directory, "enrichments")
        if not os.path.exists(enrichment_dir):
            os.makedirs(enrichment_dir)
        for key in (self.enrichments):
            fname = os.path.join(enrichment_dir, key + ".tsv")
            handle = open(fname, "w")
            header = ["sequence"]
            header += ["lib%d.count\tlib%d.freq\tlib%d.ratio" % (i + 1, i + 1, i + 1) for i in xrange(len(self.libraries))]
            header += EnrichmentStats._fields
            print("\t".join(header), file=handle)
            for (k, counts), (_, frequencies), (_, ratios) in \
                    zip(self.counter_data(key), self.frequency_data(key), self.ratio_data(key)):
                values = [str(x) for x in itertools.chain.from_iterable( \
                            zip(counts, frequencies, ratios))]
                stats = [str(x) for x in self.enrichments[key][k]]
                print(k, "\t".join(values), "\t".join(stats), sep="\t", 
                      file=handle)


    def frequency_data(self, key):
        """
        Generator function for getting frequency data from SeqLibs.

        Yields a list of values from the SeqLib frequency dicts. If a given 
        value doesn't have data in one of SeqLibs, nan is given. All lists 
        of values therefore have the same length (the number of SeqLibs).
        """
        frequencies = [lib.frequencies[key] for lib in self.libraries]
        allkeys = set()
        for f in frequencies:
            allkeys.update(f.keys())

        for k in sorted(allkeys):
            values = list()
            for f in frequencies:
                if k in f:
                    values.append(f[k])
                else:
                    values.append(float("NaN"))
            yield k, values


    def ratio_data(self, key):
        """
        Generator function for getting ratio data from SeqLibs.

        Yields a list of values from the SeqLib ratio dicts. If a given 
        value doesn't have data in one of SeqLibs, nan is given. All lists 
        of values therefore have the same length (the number of SeqLibs).
        """
        ratios = [lib.ratios[key] for lib in self.libraries]
        allkeys = set()
        for r in ratios:
            allkeys.update(r.keys())

        for k in sorted(allkeys):
            values = list()
            for r in ratios:
                if k in r:
                    values.append(r[k])
                else:
                    values.append(float("NaN"))
            yield k, values


    def counter_data(self, key):
        """
        Generator function for getting count data from SeqLibs.

        Yields a list of values from the SeqLib counters. If a given value 
        doesn't have data in one of SeqLibs, 0 is given. All lists of values 
        therefore have the same length (the number of SeqLibs).
        """
        counters = [lib.counters[key] for lib in self.libraries]
        allkeys = set()
        for c in counters:
            allkeys.update(c.keys())

        for k in sorted(allkeys):
            values = list()
            for c in counters:
                if k in c:
                    values.append(c[k])
                else:
                    values.append(0)
            yield k, values


    def calc_enrichments(self, key):
        self.enrichments[key] = dict()
        for k, values in self.ratio_data(key):
            if math.isnan(values[0]): # must be present in first library
                score = float("NaN")
            else:
                values = [x for x in values if not math.isnan(x)]
                if len(values) == 1: # can't calculate
                    score = float("NaN")
                    r = float("NaN")
                elif len(values) == 2: # simple slope
                    score = values[1] - values[0]
                    r = float("NaN")
                else:
                    xs = range(0, len(values))
                    # slope, intercept, r_value, p_value, std_err
                    score, _, r, _, _ = stats.linregress(xs, values)
            self.enrichments[key][k] = EnrichmentStats(score, r_sq=r ** 2)


    def filter_barcodes(self):
        for variant, bc_list in self.barcode_map.variants.iteritems():
            for lib in self.libraries:
                freqs = list()
                for bc in bc_list:
                    if bc in lib.frequencies['barcodes']:
                        freqs.append(lib.frequencies['barcodes'][bc])
                    else:
                        freqs.append(0.0)
            variant_cov = stats.variation(freqs)
#            if variant_cov > self.


    def exclude_variant(self, variant):
        """
        Move a variant's data to the excluded list.

        Removes the variant's enrichment information from the 'variants' data 
        structure and adds it to the 'excluded_variants' data structure. 
        Removes the variant's count and frequency information from all owned
        SeqLibs. Does nothing of the variant is not found in this Selection.
        """
        if variant in self.enrichments['variants']:
            self.enrichments['excluded_variants'][variant] = \
                    self.enrichments['variants'][variant]
            del self.enrichments['variants'][variant]
        for lib in self.libraries:
            lib.exclude_variant(variant)



