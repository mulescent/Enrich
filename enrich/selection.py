from __future__ import print_function
from enrich_error import EnrichError
from scipy import stats
from seqlib.basic import BasicSeqLib
from seqlib.barcode import BarcodeSeqLib, BarcodeMap
from seqlib.overlap import OverlapSeqLib
import os.path
import math
import itertools


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
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, self.name)

        self.enrichments = dict()
        for key in self.libraries[0].frequencies:
            if all(key in lib.frequencies for lib in self.libraries):
                self.enrichments[key] = dict()


    def enable_logging(self, log):
        self.verbose = True
        self.log = log
        try:
            print("# Logging started for '%s': %s" % 
                        (self.name, time.asctime()), file=self.log)
        except (IOError, AttributeError):
            raise EnrichException("Could not write to log file", self.name)


    def check_wt(self):
        try:
            coding = [lib['wild type']['coding'] for lib in self.libraries]
            sequence = [lib['wild type']['sequence'] for lib in self.libraries]
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, self.name)
        
        if len(set(coding)) != 1 or len(set(sequence)) != 1:
            raise EnrichError("Inconsistent wild type sequences", self.name)


    def count(self):
        for lib in self.libraries:
            lib.count()
            lib.count_mutations()
            lib.calc_frequencies()

        self.correct_frequencies()
        for key in self.enrichments:
            self.calc_enrichments[key]
        if 'barcodes' in self.enrichments:
            self.filter_barcodes()

    
    def correct_frequencies(self):
        if self.correction is not None:
            if self.correction['method'] == "stop":
                maxpos = len(self.libraries[0].wt_protein) * 
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
        for key in (self.enrichments):
            fname = os.path.join(directory, "enrichments", key + ".tsv")
            handle = open(fname, "w")
            header = ["sequence"]
            header += ["lib%d.count\tlib%d.freq" % (i + 1, i + 1) for i in 
                       xrange(len(self.libraries))]
            header += ["score"]
            print("\t".join(header), file=handle)
            for k, frequencies, _, counts in 
                    itertools.chain.from_iterable(zip(
                        self.frequency_data(key), self.counter_data(key))):
                values = itertools.chain(zip(counts, frequencies))
                print(k, "\t".join(values), self.enrichments[key][k])


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
        edict = dict()
        for k, values in self.frequency_data(key):
            if math.isnan(values[0]): # must be present in first library
                score = float("NaN")
            else:
                values = [x for x in values if not math.isnan(x)]
                if len(values) == 1: # can't calculate
                    score = float("NaN")
                elif len(values) == 2: # ratio
                    score = values[1] / float(values[0])
                else:
                    xs = range(0, len(values))
                    # slope, intercept, r_value, p_value, std_err
                    score, _, _, _, _ = stats.linregress(xs, values)
            edict[k] = score

        self.enrichments[key] = edict


    def filter_barcodes(self):

        # filter barcodes on consistency
        pass


