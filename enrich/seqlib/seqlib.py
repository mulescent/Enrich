from __future__ import print_function
import time
import pandas as pd
import base64
from sys import stdout, stderr
from enrich_error import EnrichError
import os.path


class SeqLib(object):
    """
    Abstract class for handling count data from a single sequencing library.
    Creating a *SeqLib* requires a valid *config* object, usually from a  
    ``.json`` configuration file.
    """

    _filter_messages = {
            'remove unresolvable' : "unresolvable mismatch",
            'min quality' : "single-base quality",
            'avg quality' : "average quality",
            'max mutations' : "excess mutations",
			'chastity' : "not chaste",
            'remove overlap indels' : "indel in read overlap",
            'fuser failure' : "unable to fuse reads"
    }


    def __init__(self, config):
        self.name = "Unnamed" + self.__class__.__name__
        self.verbose = False
        self.log = None

        try:
            self.name = config['name']
            self.timepoint = int(config['timepoint'])
            if 'align variants' in config:
                if config['align variants']:
                    self.aligner = Aligner()
                else:
                    self.aligner = None
            else:
                self.aligner = None

        except KeyError as key:
            raise EnrichError("Missing required config value '%s'" % key, 
                              self.name)
        except ValueError as value:
            raise EnrichError("Invalid parameter value %s" % value, self.name)

        # initialize data
        self.counts = dict()        # pandas dataframes
        self.counts_file = dict()   # paths to pickled counts
        self.filters = None         # dictionary
        self.filter_stats = None    # dictionary


    def enable_logging(self, log):
        """
        Turns on log output for this object. Messages will be sent to the 
        open file handle *log*. 

        .. note:: One log file is usually shared by all objects in the \
        analysis. This method is invoked by **Selection** logging functions.
        """
        self.verbose = True
        self.log = log
        try:
            print("# Logging started for '%s': %s" % 
                        (self.name, time.asctime()), file=self.log)
        except (IOError, AttributeError):
            raise EnrichError("Could not write to log file", self.name)


    def count(self):
        """
        This method defines how the data is counted. It must be implemented 
        by each subclass.
        """
        raise NotImplementedError("must be implemented by subclass")


    def set_filters(self, config_filters, class_default_filters):
        """
        Sets the filtering options using the values from the 
        *config_filters* dictionary and *class_default_filters* dictionary. 
        This method is used by the ``__init__`` method of *SeqLib* subclasses.

        .. note:: To help prevent user error, *config_filters* must be a \
        subset of *class_default_filters*.
        """
        self.filters = class_default_filters

        for key in self.filters:
            if key in config_filters:
                self.filters[key] = int(config['filters'][key])

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


    def report_filtered_read(self, handle, fq, filter_flags):
        """
        Outputs the **FQRead** object *fq* to *handle* (usually the object's
        log file). The dictionary *filter_flags* contains ``True`` values for 
        each filtering option that applies to *fq*. Keys in *filter_flags* 
        are converted to messages using the ``SeqLib._filter_messages`` 
        dictionary.
        """
        print("Filtered read (%s) [%s]" % \
                (', '.join(SeqLib._filter_messages[x] 
                 for x in filter_flags if filter_flags[x]), self.name), 
              file=self.log)
        print(fq, file=self.log)


    def report_filtered(self, handle):
        """
        Outputs a summary of filtered reads to *handle*. Internal filter 
        names are converted to messages using the ``SeqLib._filter_messages`` 
        dictionary.
        """
        print('Number of filtered reads for "%s"' % self.name, file=handle)
        for key in sorted(self.filter_stats):
            if key == 'total':
                continue # print this last
            if self.filter_stats[key] > 0:
                print("", SeqLib._filter_messages[key], self.filter_stats[key], 
                      sep="\t", file=handle)

        print("", 'total', self.filter_stats['total'], sep="\t", file=handle)


    def save_counts(self, directory, keys=None, clear=False):
        """
        Save the counts DataFrame as a tab-separated file in *directory*. The 
        counts are saved as ``barcode.tsv`` and/or ``variant.tsv`` as 
        appropriate for the analysis. The file names are stored in the 
        ``self.counts_file`` dictionary.

        The optional *keys* parameter is a list of types of counts to be 
        saved. By default, all counts are saved.

        If *clear* is ``True``, saved counts will be set to ``None`` after 
        writing. This is used to save memory. If needed later, the counts 
        can be restored using ``load_counts``.
        """
        if keys is None:
            keys = self.counts.keys()
        for key in keys:
            output_dir = os.path.join(directory, key)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            fname = "".join(c for c in self.name if c.isalnum() or c in (' ._~'))
            fname = fname.replace(' ', '_')
            self.counts_file[key] = os.path.join(output_dir, fname + ".txt")
            self.counts[key].to_csv(self.counts_file[key], 
                    sep="\t", na_rep="NaN", float_format="%.4g", 
                    index_label="sequence")
            if clear:
                self.counts[key] = None


    def load_counts(self, keys=None):
        """
        Load the counts from the ``.tsv`` files in the ``self.counts_file`` 
        dictionary.

        The optional *keys* parameter is a list of types of counts to be 
        loaded. By default, all counts are loaded.
        """
        if keys is None:
            keys = self.counts_file.keys()
        else:
            if not all(key in self.counts_file.keys() for key in keys):
                raise EnrichError("Cannot load unsaved counts", self.name)
        for key in keys:
            self.counts[key] = pd.from_csv(self.counts_file[key], sep="\t")

