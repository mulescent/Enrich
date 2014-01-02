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


def barcode_variation_apply_fn(row, barcode_data, mapping):
    """
    ``pandas`` apply function for calculating the coefficient of variation 
    for a variant's barcodes.
    """
    bc_scores = barcode_data.ix[mapping.variants[row.name]]['score']
    bc_scores = bc_scores[np.invert(np.isnan(bc_scores))]
    cv = stats.variation(bc_scores)
    return pd.Series({'scored.unique.barcodes' : len(bc_scores), \
                      'barcode.cv' : cv})


def barcode_count_apply_fn(row, mapping):
    """
    ``pandas`` apply function for counting the number of unique barcodes for a
    variant.
    """
    return len(mapping.variants[row.name])


def barcode_varation_filter(row, cutoff):
    """
    Filtering function for barcode coefficient of variation.
    """
    if row['barcode.cv'] > cutoff:
        return False
    else:
        return True


def min_count_filter(row, cutoff):
    """
    Filtering function for minimum counts across all timepoints.
    """
    counts = row[[x for x in row.index if x.startswith("count")]].values
    if counts.min() < cutoff:
        return False
    else:
        return True


def min_input_count_filter(row, cutoff):
    """
    Filtering function for minimum count in input timepoint.
    """
    if row['count.0'] < cutoff:
        return False
    else:
        return True


def min_rsq_filter(row, cutoff):
    """
    Filtering function for minimum r-squared value. Entries with no 
    r-squared value are retained.
    """
    if np.isnan(row['r_sq']):
        return True
    elif row['r_sq'] < cutoff:
        return False
    else:
        return True


def enrichment_apply_fn(row, timepoints):
    """
    ``pandas`` apply function for calculating enrichment scores and 
    r-squared values.
    """
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
    """
    Class for a single selection experiment, consisting of multiple timepoints. 
    This class coordinates :py:class:`SeqLib` objects. Creating a :py:class:`Selection` requires a valid *config* object, usually from a  
    ``.json`` configuration file.
    """
    def __init__(self, config):
        self.name = "Unnamed" + self.__class__.__name__
        self.libraries = dict()
        self.data = dict()
        self.data_file = dict()
        self.timepoints = list()
        self.verbose = False
        self.log = None
        self.filters = None
        self.filter_stats = None

        # PARAMETERIZE
        self.hdf_dir = "."

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

            self.set_filters(config, {'min count' : 0,
                                      'min input count' : 0,
                                      'min rsquared' : 0.0,
                                      'max barcode variation' : None})

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
                self.data[dtype] = True
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
        """
        Turns on log output for this object. Messages will be sent to the 
        open file handle *log*. 

        .. note:: One log file is usually shared by all objects in the \
        analysis. This method is invoked by :py:class:`Experiment` logging functions.
        """
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


    def set_filters(self, config_filters, default_filters):
        """
        Sets the filtering options using the values from the 
        *config_filters* dictionary and *default_filters* dictionary. 
        This method is used by the ``__init__`` method of *SeqLib* subclasses.

        .. note:: To help prevent user error, *config_filters* must be a \
        subset of *default_filters*.
        """
        self.filters = default_filters

        for key in self.filters:
            if key in config_filters:
                self.filters[key] = config_filters[key]

        unused = list()
        for key in config_filters:
            if key not in self.filters:
                unused.append(key)
        if len(unused) > 0:
            raise EnrichError("Unused filter parameters (%s)" % 
                              ', '.join(unused), self.name)

        self.filter_stats = dict()
        for key in self.filters:
            self.filter_stats[key] = 0
        self.filter_stats['total'] = 0


    def count_timepoints(self):
        """
        Combine :py:class:`SeqLib` objects into individual timepoints and 
        tabulate counts for each timepoint. Counts are stored in the local
        ``pandas`` DataFrame. To tabulate counts for individual mutations 
        (not variants), see :py:meth:`count_mutations`.
        """
        # calculate counts for each SeqLib
        for tp in self.timepoints:
            for lib in self.libraries[tp]:
                lib.count()
        for dtype in self.data:
            self.calc_counts(dtype)
        for tp in self.timepoints:
            for lib in self.libraries[tp]:
                lib.save_counts(self.hdf_dir, clear=True)


    def calc_counts(self, dtype):
        """
        Tabulate counts for each timepoint and create the DataFrame indicated by 
        *dtype* ('variant' or 'barcode'). All :py:class:`SeqLib` objects need to 
        be counted before calling this method.
        """
        # combine all libraries for a given timepoint
        self.data[dtype] = dict()
        for tp in self.timepoints:
            tp_data = self.libraries[tp][0].counts[dtype]
            for lib in self.libraries[tp][1:]:
                tp_data = tp_data.add(lib.counts[dtype], fill_value=0)
            self.data[dtype][tp] = tp_data

        tp_frame = self.data[dtype][0]
        cnames = ["%s.0" % x for x in self.data[dtype][0].columns]
        for tp in self.timepoints[1:]:
            tp_frame = tp_frame.join(self.data[dtype][tp], how="outer", 
                                     rsuffix="%s" % tp)
            cnames += ["%s.%d" % (x, tp) for x in \
                       self.data[dtype][tp].columns]
        tp_frame.columns = cnames
        # remove data that are not in the initial timepoint
        tp_frame = tp_frame[np.invert(np.isnan(tp_frame['count.0']))]
        self.data[dtype] = tp_frame


    def count_mutations(self):
        # needs to happen all the filtering/exclusion of variants
        for lib in self.libraries:
            lib.count_mutations()

        for dtype in ('mutations_nt', 'mutations_aa'):
            if dtype in self.data.keys():
                self.calc_counts(dtype)
                self.calc_frequencies(dtype)
                self.calc_ratios(dtype)
                self.calc_enrichments(dtype)


    def calc_all(self):
        """
        Wrapper method to calculate counts, frequencies, ratios, and enrichment scores 
        for all data in the :py:class:`Selection`.
        """
        self.get_counts()
        for dtype in self.data:
            self.calc_frequencies(dtype)
            self.calc_ratios(dtype)
            self.calc_enrichments(dtype)
        if 'variants' in self.data and 'barcodes' in self.data:
            self.calc_barcode_variation()


    def calc_frequencies(self, dtype):
        """
        Calculate frequencies for each element in the DataFrame indicated by 
        *dtype* ('variant' or 'barcode').
        """
        for tp in self.timepoints:
            self.data[dtype]['frequency.%d' % tp] =  \
                self.data[dtype]['count.%d' % tp] / \
                float(self.data[dtype]['count.%d' % tp].sum())


    def calc_ratios(self, dtype):
        """
        Calculate ratios for each element in the DataFrame indicated by 
        *dtype* ('variant' or 'barcode'). Assumes frequencies have been 
        calculated by :py:meth:`calc_frequencies`.
        """
        for tp in self.timepoints:
            if tp == 0: # input library
                self.data[dtype]['ratio.%d' % tp] = 1.0
            else:
                self.data[dtype]['ratio.%d' % tp] =  \
                        self.data[dtype]['frequency.%d' % tp] / \
                        self.data[dtype]['frequency.0']


    def calc_enrichments(self, dtype):
        """
        Calculate enrichment scores and r-squared values for each element in the DataFrame indicated by 
        *dtype* ('variant' or 'barcode'). Assumes ratios have been 
        calculated by :py:meth:`calc_ratios`. Calculations performed using 
        :py:func:`enrichment_apply_fn`.
        """
        # apply the enrichment-calculating function to a DataFrame
        # containing only ratio data
        ratio_df = self.data[dtype][['ratio.%d' % x for x in self.timepoints]]
        enrichments = ratio_df.apply(enrichment_apply_fn, axis=1, args=[np.asarray(self.timepoints)])
        self.data[dtype]['score'] = enrichments['score']
        self.data[dtype]['r_sq'] = enrichments['r_sq']


    def calc_barcode_variation(self):
        """
        Calculate the `coefficient of variation <http://en.wikipedia.org/wiki/Coefficient_of_variation>`_ 
        for each variant's barcode enrichment scores. Requires both variant and barcode 
        data for all timepoints.
        """
        self.data['variants']['barcode.count'] = \
                self.data['variants'].apply(barcode_count_apply_fn, 
                axis=1, args=[self.barcode_map]).astype("int32")
        barcode_cv = self.data['variants'].apply(\
                barcode_variation_apply_fn, axis=1,
                args=[self.data['barcodes'], self.barcode_map])
        self.data['variants']['scored.unique.barcodes'] = \
                barcode_cv['scored.unique.barcodes'].astype("int32")
        self.data['variants']['barcode.cv'] = barcode_cv['barcode.cv']


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


    def save_data(self, directory, keys=None, clear=False):
        if keys is None:
            keys = self.data.keys()
        for key in keys:
            output_dir = os.path.join(directory, key)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            fname = "".join(c for c in self.name if c.isalnum() or c in (' ._~'))
            fname = fname.replace(' ', '_')
            self.data_file[key] = os.path.join(output_dir, fname + ".h5")
            self.data[key].to_csv(self.data_file[key], 
                    sep="\t", na_rep="NaN", float_format="%.4g", 
                    index_label="sequence")
            if clear:
                self.data[key] = None


    def load_data(self):
        for key in self.data_file:
            self.data[key] = pd.from_csv(self.data_file[key], sep="\t")


    def filter_data(self):
        self.save_data(os.path.join(self.hdf_dir, "selection_prefilter"), 
                       clear=False)
        # for each filter that's specified
        # apply the filter
        if self.filters['max barcode variation']:
            nrows = len(self.data['variants'])
            self.data['variants'] = \
                    self.data['variants'][self.data['variants'].apply(\
                        barcode_varation_filter, axis=1, 
                        args=[self.filters['max barcode variation']])]
            self.filter_stats['max barcode variation'] = \
                    nrows - len(self.data['variants'])
        if self.filters['min count'] > 0:
            nrows = len(self.data['variants'])
            self.data['variants'] = \
                    self.data['variants'][self.data['variants'].apply(\
                        min_count_filter, axis=1, 
                        args=[self.filters['min count']])]
            self.filter_stats['min count'] = \
                    nrows - len(self.data['variants'])
        if self.filters['min input count'] > 0:
            nrows = len(self.data['variants'])
            self.data['variants'] = \
                    self.data['variants'][self.data['variants'].apply(\
                        min_input_count_filter, axis=1, 
                        args=[self.filters['min count']])]
            self.filter_stats['min input count'] = \
                    nrows - len(self.data['variants'])
        if self.filters['min rsquared'] > 0.0:
            nrows = len(self.data['variants'])
            self.data['variants'] = \
                    self.data['variants'][self.data['variants'].apply(\
                        min_rsq_filter, axis=1, 
                        args=[self.filters['min rsquared']])]
            self.filter_stats['min rsquared'] = \
                    nrows - len(self.data['variants'])

        self.filter_stats['total'] = sum(self.filter_stats.values())

        # recalculate with the updated (post-filter) frequencies
        self.calc_frequencies('variants')
        self.calc_ratios('variants')
        self.calc_enrichments('variants')

