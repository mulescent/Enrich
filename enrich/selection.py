from __future__ import print_function
from enrich_error import EnrichError
from scipy import stats
from seqlib.basic import BasicSeqLib
from seqlib.barcodevariant import BarcodeVariantSeqLib, BarcodeMap
from seqlib.overlap import OverlapSeqLib
import os
import math
import itertools
import time
import pandas as pd
import numpy as np
from collections import Counter


from sys import stdout, stderr


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
        self.libraries = dict()
        self.data = dict()
        self.timepoints = list()
        self.verbose = False
        self.log = None

        # PARAMETERIZE
        self.pickle_dir = "."

        try:
            self.name = config['name']
            if 'barcodes' in config:
                if 'map file' in config['barcodes']:
                    self.barcode_map = BarcodeMap(config['barcodes']
                                                        ['map file'])
                else:
                    self.barcode_map = None

            libnames = list()
            for lib in config['libraries']:
                if 'barcodes' in lib:
                    if 'wild type' in lib:
                        new = BarcodeVariantSeqLib(lib, 
                                    barcode_map=self.barcode_map)
                    else:
                        new = BarcodeSeqLib(lib)
                elif 'overlap' in lib:
                    new = OverlapSeqLib(lib)
                else:
                    new = BasicSeqLib(lib)

                if new.timepoint not in self.libraries:
                    self.libraries[new.timepoint] = list()
                self.libraries[new.timepoint].append(new)
                libnames.append(new.name)
            self.timepoints = sorted(self.libraries.keys())

            if len(set(libnames)) != len(libnames):
                raise EnrichError("Non-unique library names", self.name)

        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, 
                              self.name)

        if len(self.libraries.keys()) < 2:
            raise EnrichError("Insufficient number of timepoints", 
                              self.name)

        if 0 not in self.timepoints:
            raise EnrichError("Missing timepoint 0", self.name)
        if self.timepoints[0] != 0:
            raise EnrichError("Invalid negative timepoint", self.name)

        # identify what kind of counts data is present in all timepoints
        dtype_counts = list()
        for tp in self.timepoints:
            for lib in self.libraries[tp]:
                dtype_counts.extend(lib.counts.keys())
        dtype_counts = Counter(dtype_counts)
        for dtype in dtype_counts:
            if dtype_counts[dtype] == len(config['libraries']):
                self.data[dtype] = dict()
        if 'barcodes_unmapped' in self.data.keys(): # special case for BarcodeVariantSeqLib
            del self.data['barcodes_unmapped']
        if len(self.data.keys()) == 0:
            raise EnrichError("No count data present across all timepoints", 
                              self.name)

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
        for timepoint in self.libraries:
            for lib in self.libraries[timepoint]:
                lib.enable_logging(log)


    def get_counts(self):
        # calculate counts for each timepoint
        for tp in self.timepoints:
            for lib in self.libraries[tp]:
                lib.count()
            for dtype in self.data:
                tp_data = self.libraries[tp][0].counts[dtype]
                for lib in self.libraries[tp][1:]:
                    tp_data = tp_data.add(lib.counts[dtype], fill_value=0)
                self.data[dtype][tp] = tp_data
            for lib in self.libraries[tp]:
                lib.save_counts(self.pickle_dir)

        # combine into a single dataframe
        for dtype in self.data:
            tp_data = self.data[dtype][0]
            cnames = ["%s.0" % x for x in self.data[dtype][0].columns]
            for tp in self.timepoints[1:]:
                tp_data = tp_data.join(self.data[dtype][tp], 
                                     how="outer",
                                     rsuffix="%s" % tp)
                cnames += ["%s.%d" % (x, tp) for x in \
                           self.data[dtype][tp].columns]
            tp_data.columns = cnames
            # remove data that are not in the initial timepoint
            tp_data = tp_data[np.invert(np.isnan(tp_data['count.0']))]
            self.data[dtype] = tp_data



    def calc_frequencies(self):
        for dtype in self.data:
            for tp in self.timepoints:
                self.data[dtype]['frequency.%d' % tp] =  \
                    self.data[dtype]['count.%d' % tp] / \
                    float(self.data[dtype]['count.%d' % tp].sum())


    def calc_ratios(self):
        for dtype in self.data:
            for tp in self.timepoints:
                if tp == 0: # input library
                    self.data[dtype]['ratio.%d' % tp] = 1.0
                else:
                    self.data[dtype]['ratio.%d' % tp] =  \
                            self.data[dtype]['frequency.%d' % tp] / \
                            self.data[dtype]['frequency.0']


    def calc_enrichments(self):
        for dtype in self.data:
            # apply the enrichment-calculating function to a DataFrame
            # containing only ratio data
            ratio_df = self.data[dtype][['ratio.%d' % x for x in self.timepoints]]
            enrichments = ratio_df.apply(enrichment_apply_func, axis=1, args=[np.asarray(self.timepoints)])
            self.data[dtype] = self.data[dtype].join(enrichments)


    def count_mutations(self):
        # needs to happen all the filtering/exclusion of variants
        for lib in self.libraries:
            lib.count_mutations()

        for key in ('mutations_nt', 'mutations_aa'):
            self.enrichments[key] = dict()
            self.calc_enrichments(key)


    def adjust_frequencies(self, frequencies, timepoint, dtype):
        """
        Replace frequencies with a corrected frequency for given timepoint.

        Corrected values overwrite the "raw" frequencies in the data
        DataFrame. To adjust all frequencies in the DataFrame, this method 
        must be invoked once for each timepoint.

        Frequencies must be a DataFrame with one column. Indices present in 
        freuqences that are not present in the Selection DataFrame will not 
        be added to the DataFrame.
        """
        if timepoint not in self.timepoints:
            raise EnrichError("Timepoint not found", self.name)
        else:
            self.data[dtype]['frequency.%d' % timepoint] = frequencies[0]


    def write_enrichments(self, directory):
        enrichment_dir = os.path.join(directory, "enrichments")
        if not os.path.exists(enrichment_dir):
            os.makedirs(enrichment_dir)
        for dtype in self.data:
            fname = os.path.join(enrichment_dir, dtype + ".tsv")
            self.data[dtype].to_csv(fname, sep="\t", na_rep="NaN", 
                             index_label="sequence")


    def filter_barcodes(self):
        pass
