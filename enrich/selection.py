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
import pandas as pd
import numpy as np
from collections import namedtuple


from sys import stdout, stderr


EnrichmentStats = namedtuple("EnrichmentStats", "score, r_sq")

def enrichment_apply_func(row, timepoints):
    if math.isnan(row[0]):
        # not present in input library
        score = float("NaN")
        r_sq = float("NaN")
    else:
        row = row.values
        ratios = row[np.invert(np.isnan(row))]
        times = timepoints[np.invert(np.isnan(row))]
        if len(ratios) == 1:
            # only present in input library
            score = float("NaN")
            r_sq = float("NaN")
        elif len(ratios) == 2:
            # rise over run
            score = (ratios[1] - ratios[0]) / (times[1] - times[0])
            r_sq = float("NaN")
        else:
            score, _, r, _, _ = stats.linregress(times, ratios)
            r_sq = r ** 2

    return pd.Series({'score' : score, 'r_sq' : r_sq})


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
            if len(set([lib.name for lib in self.libraries])) != \
                    len(self.libraries):
                raise EnrichError("Non-unique library names", self.name)
            self.libraries.sort(key=lambda x: x.timepoint)

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
            lib.calc_ratios(self.libraries[0])

#        self.correct_frequencies()
        for key in self.libraries[0].data:
            if all(key in lib.data for lib in self.libraries):
                self.calc_enrichments(key)
        if 'barcodes' in self.enrichments:
            self.filter_barcodes()


    def count_mutations(self):
        # needs to happen all the filtering/exclusion of variants
        for lib in self.libraries:
            lib.count_mutations()

        for key in ('mutations_nt', 'mutations_aa'):
            self.enrichments[key] = dict()
            self.calc_enrichments(key)


    def correct_frequencies(self):
        pass

    def write_enrichments(self, directory):
        enrichment_dir = os.path.join(directory, "enrichments")
        if not os.path.exists(enrichment_dir):
            os.makedirs(enrichment_dir)
        for key in (self.enrichments):
            fname = os.path.join(enrichment_dir, key + ".tsv")
            # create a new DataFrame for output
            # start with the first SeqLibs
            output_df = self.libraries[0].data[key].drop(['ratio', 'excluded'], axis=1)
            # append additional SeqLibs
            for lib in self.libraries[1:]:
                output_df = output_df.join(lib.data[key].drop('excluded', axis=1), how='outer',
                        rsuffix=".lib%d" % lib.timepoint)
            # append enrichment data
            output_df = output_df.join(self.enrichments[key], how='outer')
            output_df = output_df[np.invert(output_df['excluded'])]

            # write the file
            output_df.to_csv(fname, sep="\t", na_rep="NaN", 
                             index_label="sequence")


    def calc_enrichments(self, key):
        # make a DataFrame containing all the ratios
        ratio_df = self.libcolumns(key, 'ratio')
        timepoints = np.asarray([lib.timepoint for lib in self.libraries])
        self.enrichments[key] = ratio_df.apply(enrichment_apply_func, axis=1, 
                                               args=[timepoints])
        self.enrichments[key]['excluded'] = False


    def libcolumns(self, key, column):
        lib_df = pd.DataFrame(self.libraries[0].data[key][column])
        cnames = ['%s.lib%d' % (column, self.libraries[0].timepoint)]
        for lib in self.libraries[1:]:
            lib_df = lib_df.join(pd.DataFrame(lib.data[key][column] \
                    [np.invert(lib.data[key]['excluded'])]),
                    rsuffix="%d" % lib.timepoint, how='outer')
            cnames.append('%s.lib%d' % (column, lib.timepoint))
        lib_df.columns = cnames
        return lib_df


    def filter_barcodes(self):
        pass
