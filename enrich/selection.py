from __future__ import print_function
from enrich_error import EnrichError
from scipy import stats
from seqlib.basic import BasicSeqLib
from seqlib.barcodevariant import BarcodeVariantSeqLib, BarcodeMap
from seqlib.barcode import BarcodeSeqLib
from seqlib.overlap import OverlapSeqLib
from config_check import seqlib_type
import os
import re
import math
import itertools
import time
import pandas as pd
import numpy as np
from collections import Counter


def nonsense_ns_carryover_apply_fn(row, position):
    """
    :py:meth:`pandas.DataFrame.apply` function for determining which rows 
    contribute counts to nonspecific carryover calculations. Returns ``True`` 
    if the variant has a change to stop at or before amino acid number 
    *position*.
    """
    m = re.search("p\.[A-Z][a-z][a-z](\d+)Ter", row.name)
    if m is not None:
        if int(m.group(1)) <= position:
            return True
        else:
            return False
    else:
        return False


def barcode_variation_apply_fn(row, barcode_data, mapping):
    """
    :py:meth:`pandas.DataFrame.apply` function for calculating the coefficient 
    of variation for a variant's barcodes.
    """
    bc_scores = barcode_data.ix[mapping.variants[row.name]]['score']
    bc_scores = bc_scores[np.invert(np.isnan(bc_scores))]
    cv = stats.variation(bc_scores)
    return pd.Series({'scored.unique.barcodes' : len(bc_scores), \
                      'barcode.cv' : cv})


def barcode_count_apply_fn(row, mapping):
    """
    :py:meth:`pandas.DataFrame.apply` function for counting the number of 
    unique barcodes for a variant.
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
    :py:meth:`pandas.DataFrame.apply` apply function for calculating 
    enrichment scores and r-squared values.
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


class Selection(DataContainer):
    """
    Class for a single selection experiment, consisting of multiple 
    timepoints. This class coordinates :py:class:`~seqlib.seqlib.SeqLib` 
    objects. Creating a :py:class:`~selection.Selection` requires a valid 
    *config* object, usually from a ``.json`` configuration file.

    Example config file for a :py:class:`~selection.Selection`:

    .. literalinclude:: config_examples/selection.json

    :download:`Download this JSON file <config_examples/selection.json>`

    The ``"libraries"`` config entry is a list of all sequencing library 
    configuration entries. Counts for :py:class:`~seqlib.seqlib.SeqLib` 
    objects from the same ``timepoint`` are combined. Scores are calculated 
    based on the resulting time series. 

    All filters listed in the example config are optional. Unlike 
    :py:class:`~seqlib.seqlib.SeqLib` objects, filtered data are not written 
    to a log file. Instead, the data are written to a file before filtering 
    and stored for future comparison. 

    .. note:: ``"max barcode variation"`` filtering is only applicable if \
    all sequencing data has both barcode and variant data (*i.e.* are \
    :py:class:`~seqlib.barcodevariant.BarcodeVariantSeqLib` objects).

    If the optional ``"carryover correction"`` is specified, scores will be 
    corrected based on an estimate of nonspecific carryover (retention of 
    non-functional variants during the selection). Currently, this is  
    implemented for the ``"nonsense"`` option, as described by 
    `Araya and Fowler`_. However, additional methods can be added using this 
    as an example.
    """
    def __init__(self, config):
        DataContainer.__init__(self, config)
        self.libraries = dict()
        self.timepoints = list()

        try:
            if 'barcodes' in config:
                if 'map file' in config['barcodes']:
                    self.barcode_map = BarcodeMap(config['barcodes']
                                                        ['map file'])
                else:
                    self.barcode_map = None
            else:
                self.barcode_map = None

            libnames = list()
            for lib in config['libraries']:
                libtype = seqlib_type(lib)
                if libtype is None:
                    raise EnrichError("Unrecognized SeqLib config", self.name)
                elif libtype == "BarcodeVariantSeqLib":
                    new = BarcodeVariantSeqLib(lib, barcode_map=self.barcode_map)
                else:
                    new = globals()[libtype](lib)

                if new.timepoint not in self.libraries:
                    self.libraries[new.timepoint] = list()
                self.libraries[new.timepoint].append(new)
                libnames.append(new.name)
            self.timepoints = sorted(self.libraries.keys())

            if len(set(libnames)) != len(libnames):
                raise EnrichError("Non-unique library names", self.name)

            self.set_filters(config['filters'], {'min count' : 0,
                                      'min input count' : 0,
                                      'min rsquared' : 0.0,
                                      'max barcode variation' : None})

            if 'carryover correction' in config:
                if config['carrover correction']['method'] == "nonsense":
                    self.ns_carryover_fn = nonsense_ns_carryover_apply_fn
                    self.ns_carryover_kwargs = {'position' : int(config['carryover correction']['position'])}
                # add additional methods here using "elif" blocks
                else:
                    raise EnrichError("Unrecognized nonspecific carryover correction", self.name)
            else:
                self.ns_carryover_fn = None
                self.ns_carryover_kwargs = None

        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, 
                              self.name)
        except ValueError as value:
            raise EnrichError("Invalid parameter value %s" % value, self.name)

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
                self.df_dict[dtype] = True
        if 'barcodes_unmapped' in self.df_dict.keys(): # special case for BarcodeVariantSeqLib
            del self.df_dict['barcodes_unmapped']
        if len(self.df_dict.keys()) == 0:
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


    def count_timepoints(self):
        """
        Combine :py:class:`~seqlib.seqlib.SeqLib` objects into individual timepoints and 
        tabulate counts for each timepoint. Counts are stored in the local
        :py:class:`pandas.DataFrame`. To tabulate counts for individual mutations 
        (not variants), see :py:meth:`count_mutations`.
        """
        # calculate counts for each SeqLib
        for tp in self.timepoints:
            for lib in self.libraries[tp]:
                lib.calculate()
        for dtype in self.df_dict:
            self.calc_counts(dtype)
        for tp in self.timepoints:
            for lib in self.libraries[tp]:
                lib.save_data(self.save_dir, clear=True)


    def calc_counts(self, dtype):
        """
        Tabulate counts for each timepoint and create the :py:class:`pandas.DataFrame` indicated by 
        *dtype* ('variant' or 'barcode'). All :py:class:`~seqlib.seqlib.SeqLib` objects need to 
        be counted before calling this method.
        """
        # combine all libraries for a given timepoint
        self.df_dict[dtype] = dict()
        for tp in self.timepoints:
            tp_data = self.libraries[tp][0].counts[dtype]
            for lib in self.libraries[tp][1:]:
                tp_data = tp_data.add(lib.counts[dtype], fill_value=0)
            self.df_dict[dtype][tp] = tp_data

        tp_frame = self.df_dict[dtype][0]
        cnames = ["%s.0" % x for x in self.df_dict[dtype][0].columns]
        for tp in self.timepoints[1:]:
            tp_frame = tp_frame.join(self.df_dict[dtype][tp], how="outer", 
                                     rsuffix="%s" % tp)
            cnames += ["%s.%d" % (x, tp) for x in \
                       self.df_dict[dtype][tp].columns]
        tp_frame.columns = cnames
        # remove data that are not in the initial timepoint
        tp_frame = tp_frame[np.invert(np.isnan(tp_frame['count.0']))]
        self.df_dict[dtype] = tp_frame


    def count_mutations(self):
        """
        Creates and populates :py:class:`pandas.DataFrame` objects for individual mutations. This method 
        should be called after all filtering has been completed. The new 
        :py:class:`pandas.DataFrame` objects have dtype 'mutations_nt' and 'mutations_aa' (only if the 
        data set is coding).
        """
        # needs to happen all the filtering/exclusion of variants
        for lib in self.libraries:
            lib.count_mutations()

        for dtype in ('mutations_nt', 'mutations_aa'):
            if dtype in self.df_dict.keys():
                self.calc_counts(dtype)
                self.calc_frequencies(dtype)
                self.calc_ratios(dtype)
                self.calc_enrichments(dtype)


    def calculate(self):
        """
        Wrapper method to calculate counts, frequencies, ratios, and enrichment scores 
        for all data in the :py:class:`Selection`.
        """
        self.count_timepoints()
        if self.ns_carryover_fn is not None:
            self.nonspecific_carryover(self.ns_carryover_fn, **self.ns_carryover_kwargs)
        for dtype in self.df_dict:
            self.calc_frequencies(dtype)
            self.calc_ratios(dtype)
            self.calc_enrichments(dtype)
        if 'variants' in self.df_dict and 'barcodes' in self.df_dict:
            self.calc_barcode_variation()


    def calc_frequencies(self, dtype):
        """
        Calculate frequencies for each element in the :py:class:`pandas.DataFrame` indicated by 
        *dtype* ('variant' or 'barcode').
        """
        for tp in self.timepoints:
            self.df_dict[dtype]['frequency.%d' % tp] =  \
                self.df_dict[dtype]['count.%d' % tp] / \
                float(self.df_dict[dtype]['count.%d' % tp].sum())


    def calc_ratios(self, dtype):
        """
        Calculate ratios for each element in the :py:class:`pandas.DataFrame` indicated by 
        *dtype* ('variant' or 'barcode'). Assumes frequencies have been 
        calculated by :py:meth:`calc_frequencies`.
        """
        for tp in self.timepoints:
            if tp == 0: # input library
                self.df_dict[dtype]['ratio.%d' % tp] = 1.0
            else:
                self.df_dict[dtype]['ratio.%d' % tp] =  \
                        self.df_dict[dtype]['frequency.%d' % tp] / \
                        self.df_dict[dtype]['frequency.0']


    def calc_enrichments(self, dtype):
        """
        Calculate enrichment scores and r-squared values for each element in the :py:class:`pandas.DataFrame` indicated by 
        *dtype* ('variant' or 'barcode'). Assumes ratios have been 
        calculated by :py:meth:`calc_ratios`. Calculations performed using 
        :py:func:`enrichment_apply_fn`.
        """
        # apply the enrichment-calculating function to a DataFrame
        # containing only ratio data
        ratio_df = self.df_dict[dtype][['ratio.%d' % x for x in self.timepoints]]
        enrichments = ratio_df.apply(enrichment_apply_fn, axis=1, args=[np.asarray(self.timepoints)])
        self.df_dict[dtype]['score'] = enrichments['score']
        self.df_dict[dtype]['r_sq'] = enrichments['r_sq']


    def calc_barcode_variation(self):
        """
        Calculate the `coefficient of variation <http://en.wikipedia.org/wiki/Coefficient_of_variation>`_ 
        for each variant's barcode enrichment scores. Requires both variant and barcode 
        data for all timepoints.
        """
        self.df_dict['variants']['barcode.count'] = \
                self.df_dict['variants'].apply(barcode_count_apply_fn, 
                axis=1, args=[self.barcode_map]).astype("int32")
        barcode_cv = self.df_dict['variants'].apply(\
                barcode_variation_apply_fn, axis=1,
                args=[self.df_dict['barcodes'], self.barcode_map])
        self.df_dict['variants']['scored.unique.barcodes'] = \
                barcode_cv['scored.unique.barcodes'].astype("int32")
        self.df_dict['variants']['barcode.cv'] = barcode_cv['barcode.cv']


    def nonspecific_carryover(self, ns_apply_fn, **kwargs):
        """
        Correct the counts in the 'variants' :py:class:`pandas.DataFrame` for nonspecific carryover. 
        Nonspecific counts are defined by *ns_apply_fn* and its *kwargs*, which 
        takes a row as an argument and returns ``True`` if the row's counts 
        are nonspecific.
        """
        dtype = 'variants'
        ns_data = self.df_dict[dtype][self.df_dict[dtype].apply(ns_apply_fn, axis=1, **kwargs)]
        read_totals = [self.df_dict[dtype]['count.%d' % tp].sum() \
                       for tp in self.timepoints]
        read_totals = dict(zip(self.timepoints, read_totals))
        ns_frequencies = [ns_data['count.%d' % tp].sum() / \
                             read_totals[tp] for tp in self.timepoints]
        ns_frequencies = dict(zip(self.timepoints, ns_frequencies))
        for tp in self.timepoints[1:]: # don't modify time 0
            ns_mod = (ns_frequencies[tp] / ns_frequencies[0])
            self.df_dict[dtype]['count.%d' % tp] = self.df_dict[dtype]['count.%d' % tp] - self.df_dict[dtype]['count.%d' % tp] * ns_mod
            self.df_dict[dtype]['count.%d' % tp] = self.df_dict[dtype]['count.%d' % tp].astype("int32")


    def save_data(self, directory, keys=None, clear=False):
        """
        Save the :py:class:`pandas.DataFrame` objects as tab-separated files in *directory*. The 
        file names are stored in the ``self.df_file`` dictionary.

        The optional *keys* parameter is a list of types of counts to be 
        saved. By default, all counts are saved.

        If *clear* is ``True``, saved data will be set to ``None`` after 
        writing. This is used to save memory. If needed later, the data 
        can be restored using :py:meth:`load_data`.
        """
        if keys is None:
            keys = self.df_dict.keys()
        for key in keys:
            output_dir = os.path.join(directory, key)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            fname = "".join(c for c in self.name if c.isalnum() or c in (' ._~'))
            fname = fname.replace(' ', '_')
            self.df_file[key] = os.path.join(output_dir, fname + ".tsv")
            self.df_dict[key].to_csv(self.df_file[key], 
                    sep="\t", na_rep="NaN", float_format="%.4g", 
                    index_label="sequence")
            if clear:
                self.df_dict[key] = None


    def load_data(self):
        """
        Load the data from the ``.tsv`` files in the ``self.df_file`` 
        dictionary.

        The optional *keys* parameter is a list of types of counts to be 
        loaded. By default, all counts are loaded.
        """
        if keys is None:
            keys = self.df_file.keys()
        else:
            if not all(key in self.df_file.keys() for key in keys):
                raise EnrichError("Cannot load unsaved data", self.name)
        for key in keys:
            self.counts[key] = pd.from_csv(self.df_file[key], sep="\t")


    def filter_data(self):
        """
        Apply the filtering functions to the data, based on the filter 
        options present in the configuration object. Filtering is performed 
        using the appropriate apply function. Frequencies, ratios, and 
        enrichments must be recalculated after filtering.
        """
        self.save_data(os.path.join(self.save_dir, "selection_prefilter"), 
                       clear=False)
        # for each filter that's specified
        # apply the filter
        if self.filters['max barcode variation']:
            nrows = len(self.df_dict['variants'])
            self.df_dict['variants'] = \
                    self.df_dict['variants'][self.df_dict['variants'].apply(\
                        barcode_varation_filter, axis=1, 
                        args=[self.filters['max barcode variation']])]
            self.filter_stats['max barcode variation'] = \
                    nrows - len(self.df_dict['variants'])
        if self.filters['min count'] > 0:
            nrows = len(self.df_dict['variants'])
            self.df_dict['variants'] = \
                    self.df_dict['variants'][self.df_dict['variants'].apply(\
                        min_count_filter, axis=1, 
                        args=[self.filters['min count']])]
            self.filter_stats['min count'] = \
                    nrows - len(self.df_dict['variants'])
        if self.filters['min input count'] > 0:
            nrows = len(self.df_dict['variants'])
            self.df_dict['variants'] = \
                    self.df_dict['variants'][self.df_dict['variants'].apply(\
                        min_input_count_filter, axis=1, 
                        args=[self.filters['min count']])]
            self.filter_stats['min input count'] = \
                    nrows - len(self.df_dict['variants'])
        if self.filters['min rsquared'] > 0.0:
            nrows = len(self.df_dict['variants'])
            self.df_dict['variants'] = \
                    self.df_dict['variants'][self.df_dict['variants'].apply(\
                        min_rsq_filter, axis=1, 
                        args=[self.filters['min rsquared']])]
            self.filter_stats['min rsquared'] = \
                    nrows - len(self.df_dict['variants'])

        self.filter_stats['total'] = sum(self.filter_stats.values())


