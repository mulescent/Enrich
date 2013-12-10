from __future__ import print_function
import time
import pandas as pd
import base64
from sys import stdout, stderr
from enrich_error import EnrichError
import os.path


class SeqLib(object):
    """
    Abstract class for data from a single Enrich sequencing library.

    Implements core functionality for handling count data, parsing 
    configuration objects, and other shared processes.
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
        self.verbose = True
        self.log = log
        try:
            print("# Logging started for '%s': %s" % 
                        (self.name, time.asctime()), file=self.log)
        except (IOError, AttributeError):
            raise EnrichException("Could not write to log file", self.name)


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


    def report_filtered(self, handle):
        print('Number of filtered reads for "%s"' % self.name, file=handle)
        for key in sorted(self.filter_stats):
            if key == 'total':
                continue # print this last
            if self.filter_stats[key] > 0:
                print("", SeqLib._filter_messages[key], self.filter_stats[key], 
                      sep="\t", file=handle)

        print("", 'total', self.filter_stats['total'], sep="\t", file=handle)


    def save_counts(self, directory, keys=None):
        """
        Save the counts DataFrame as a pickled object and remove the counts.

        The counts filename is generated procedurally using base64 encoding.
        This sacrifices human readability in favor of ensuring a valid
        filename is generated regardless of the SeqLib name.
        """
        if keys is None:
            keys = self.counts.keys()
        for key in keys:
            counts_dir = os.path.join(directory, "pickled_counts", key)
            if not os.path.exists(counts_dir):
                os.makedirs(counts_dir)
            fname = "".join(c for c in self.name if c.isalnum() or c in (' ._~'))
            self.counts_file[key] = os.path.join(counts_dir, fname + ".pickle")
            self.counts[key].to_pickle(self.counts_file[key])
            self.counts[key] = None


    def load_counts(self):
        """
        Load the pickled counts.
        """
        for key in self.counts_file:
            self.counts[key] = pd.read_pickle(self.counts_file[key], table)

