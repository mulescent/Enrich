from __future__ import print_function
from enrich_error import EnrichError
from scipy import stats
import seqlib.basic
import seqlib.barcode
import seqlib.overlap
import os.path
import math


def new_seqlib(config):
    if 'barcodes' in config:
        return seqlib.barcode.BarcodeSeqLib(config)
    elif 'overlap' in config:
        return seqlib.overlap.OverlapSeqLib(config)
    else:
        return seqlib.basic.BasicSeqLib(config)


class Selection(object):
    def __init__(self, config):
        self.libraries = list()
        try:
            for lib in config['libraries']:
                self.libraries.append(new_seqlib(lib))
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key)

        if len(self.libraries) == 0:
            raise EnrichError("Selection class has no libraries")
        self.check_wt()

        if 'carryover' in config:
            if config['carryover'] in ("stop"):
                self.carryover = config['carryover']
            else:
                self.carrover = None
                raise EnrichError("Invalid carryover correction option \"%s\"" % config['carryover'])
        else:
            self.carryover = None

        self.enrichments = dict()
        for key in self.libraries[0].counters:
            if all(key in lib.counters for lib in self.libraries):
                self.enrichments[key] = dict()


    def check_wt(self):
        try:
            coding = [lib['wild type']['coding'] for lib in self.libraries]
            sequence = [lib['wild type']['sequence'] for lib in self.libraries]
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key)
        
        if len(set(coding)) != 1 or len(set(sequence)) != 1:
            raise EnrichError("Inconsistent wild type sequences")


    def count(self):
        # NOT FINISHED
        for lib in self.libraries:
            lib.count()
            lib.count_mutations()

        for key in self.enrichments:
            self.calc_enrichments[key]
        if 'barcodes' in self.enrichments:
            self.filter_barcodes()

        
    def write_enrichments(self, directory):
        for key in (self.enrichments):
            fname = os.path.join(directory, "enrichments", key + ".tsv")
            handle = open(fname, "w")
            header = ["sequence"]
            header += ["lib%d.count" % (i + 1) for i in 
                       xrange(len(self.libraries))]
            header += ["score"]
            print("\t".join(header), file=handle)
            for k, values in self.counter_data(key):
                print(k, "\t".join(values), self.enrichments[key][k])


    def counter_data(self, key):
        """
        Generator function for getting count data from SeqLibs.

        Yields a list of values from the SeqLib counters. Extracts data from
        the counters of type key. If a given value doesn't have data in one of
        the counters, nan is given. All lists of values therefore have the 
        same length (the number of SeqLibs).
        """
        counters = [lib.counters[key] for lib in self.libraries]
        allkeys = set()
        for c in counters:
            allkeys.update(c.keys())

        for k in allkeys:
            values = list()
            for c in counters:
                if k in c:
                    values.append(c[k])
                else:
                    values.append(float("NaN"))
            yield k, values


    def calc_enrichments(self, key):
        edict = dict()
        for k, values in self.counter_data(key):
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


