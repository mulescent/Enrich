from __future__ import print_function
from enrich_error import EnrichError
import sys
import time
import os

class DataContainer(object):
    """
    Abstract class for all data-containing classes 
    (:py:class:`seqlib.seqlib.SeqLib`, :py:class:`selection.Selection`, and 
    :py:class:`experiment.Experiment`). Creating a 
    :py:class:`datacontainer.DataContainer`  requires a valid *config* object, 
    usually from a ``.json`` configuration file.

    .. note:: Example configuration files can be found in the documentation \
    for derived classes.

    Log file messages for filtered reads are defined here in the 
    ``_filter_messages`` dictionary. New read filtering options must have an 
    associated log file output message added to the dictionary.

    .. literalinclude:: ../datacontainer.py
        :lines: 26-41
    """

    # Note: the following block is referenced by line number above
    # When adding new messages, update the documentation line numbers also!
    _filter_messages = {
            # SeqLib messages
            'remove unresolvable' : "unresolvable mismatch",
            'min quality' : "single-base quality",
            'avg quality' : "average quality",
            'max mutations' : "excess mutations",
            'chastity' : "not chaste",
            'remove overlap indels' : "indel in read overlap",
            'merge failure' : "unable to merge reads",
            # Selection messages
            'min count' : "not enough variant reads",
            'min input count' : "not enough variant reads in input",
            'min rsquared' : "low r-squared",
            'max barcode variation' : "high barcode CV"
            # Experiment messages
        }


    def __init__(self, config):
        self.name = "Unnamed" + self.__class__.__name__
        self.verbose = False
        self.log = None
        self.df_dict = dict()
        self.df_file = dict()
        self.filters = None
        self.filter_stats = None
        self.save_dir = "."
        
        try:
            self.name = config['name']
            
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, 
                              self.name)
                              
        #self.save_dir = config['output directory']

    def enable_logging(self, log):
        """
        Turns on log output for this object. Messages will be sent to the 
        open file handle *log*. 

        .. note:: One log file is usually shared by all objects in the \
        analysis.
        """
        self.verbose = True
        self.log = log
        try:
            print("# Logging started for '%s': %s" % 
                        (self.name, time.asctime()), file=self.log)
        except (IOError, AttributeError):
            raise EnrichException("Could not write to log file", self.name)


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


    def set_filters(self, config_filters, default_filters):
        """
        Sets the filtering options using the values from the *config_filters* 
        dictionary and *default_filters* dictionary. 

        .. note:: To help prevent user error, *config_filters* must be a \
        subset of *default_filters*.
        """
        self.filters = default_filters

        for key in self.filters:
            if key in config_filters:
                self.filters[key] = int(config_filters[key])

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


    def write_filter_stats(self, handle=None):
        if handle is None:
            if self.log is not None:
                handle = self.log
            else:
                handle = sys.stdout

        for key in sorted(self.filter_stats, key=self.filter_stats.__getitem__, reverse=True):
            if key != 'total':
                print(DataContainer._filter_messages[key], self.filter_stats[key], file=handle)
        print('total', self.filters_stats['total'], file=handle)


    def calculate(self):
        raise NotImplementedError("must be implemented by subclass")


    def filter_data(self):
        raise NotImplementedError("must be implemented by subclass")


    def make_plots(self):
        raise NotImplementedError("must be implemented by subclass")
