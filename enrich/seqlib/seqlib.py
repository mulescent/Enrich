from __future__ import print_function
import time
import pandas as pd
from enrich_error import EnrichError
from datacontainer import DataContainer
import os.path


class SeqLib(DataContainer):
    """
    Abstract class for handling count data from a single sequencing library.
    Creating a :py:class:`seqlib.seqlib.SeqLib` requires a valid *config* 
    object, usually from a ``.json`` configuration file.

    .. note:: Example configuration files can be found in the documentation \
    for derived classes.
    """
    def __init__(self, config):
        DataContainer.__init__(self, config)

        try:
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
        self.counts_file = dict()   # paths to saved counts
        self.filters = None         # dictionary
        self.filter_stats = None    # dictionary


    def calculate(self):
        """
        This method defines how the data is counted. It must be implemented 
        by each subclass.
        """
        raise NotImplementedError("must be implemented by subclass")


    def report_filtered_read(self, handle, fq, filter_flags):
        """
        Outputs the :py:class:`~fqread.FQRead` object *fq* to *handle* 
        (usually a log file). The dictionary *filter_flags* contains ``True`` 
        values for each filtering option that applies to *fq*. Keys in 
        *filter_flags* are converted to messages using the 
        ``DataContainer._filter_messages`` dictionary.
        """
        print("Filtered read (%s) [%s]" % \
                (', '.join(DataContainer._filter_messages[x] 
                 for x in filter_flags if filter_flags[x]), self.name), 
              file=self.log)
        print(fq, file=self.log)
