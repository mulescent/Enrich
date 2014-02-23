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


    def save_data(self, directory, keys=None, clear=False):
        """
        Save the counts DataFrame as a tab-separated file in *directory*. 
        The file names are stored in the ``self.counts_file`` dictionary.

        The optional *keys* parameter is a list of types of counts to be 
        saved. By default, all counts are saved.

        If *clear* is ``True``, saved counts will be set to ``None`` after 
        writing. This is used to save memory. If needed later, the counts 
        can be restored using :py:meth:`load_data`.
        """
        if keys is None:
            keys = self.counts.keys()
        for key in keys:
            output_dir = os.path.join(directory, key)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            fname = "".join(c for c in self.name if c.isalnum() or c in (' ._~'))
            fname = fname.replace(' ', '_')
            self.counts_file[key] = os.path.join(output_dir, fname + ".tsv")
            self.counts[key].to_csv(self.counts_file[key], 
                    sep="\t", na_rep="NaN", float_format="%.4g", 
                    index_label="sequence")
            if clear:
                self.counts[key] = None


    def load_data(self, keys=None):
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

