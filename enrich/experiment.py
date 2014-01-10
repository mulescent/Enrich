from __future__ import print_function
from enrich_error import EnrichError
import selection
import time


class Experiment(object):
    """
    Class for a coordinating multiple :py:class:`Selection` objects. The :py:class:`Selection` objects 
    are assigned as replicates of experimental conditions. 
    Creating an :py:class:`Experiment` requires a valid *config* object, usually from a  
    ``.json`` configuration file.
    """
    def __init__(self, config):
        self.name = "Unnamed" + self.__class__.__name__
        self.verbose = False
        self.log = None
        self.conditions = dict()
        self.control = None

        try:
            self.name = config['name']
            for cnd in config['conditions']:
                if not cnd['label'].isalnum():
                    raise EnrichError("Alphanumeric label required for condition '%s'" % cnd['label'], self.name)
                self.conditions[cnd['label']] = [selection.Selection(x) for x in cnd['selections']]
                if cnd['control']:
                    if self.control is None:
                        self.control = self.conditions[cnd['label']]
                    else:
                        raise EnrichError("Multiple control conditions", self.name)
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, 
                              self.name)


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


    def set_filters(self, config_filters, default_filters):
        """
        Sets the filtering options using the values from the 
        *config_filters* dictionary and *default_filters* dictionary. Filtering 
        options include consistency of scores across replicates.

        .. note:: To help prevent user error, *config_filters* must be a \
        subset of *default_filters*.
        """
        pass


    def calc_selection_scores(self):
        """
        Calculate scores for all :py:class:`Selection` objects.
        """
        for c in self.conditions:
            for sel in self.conditions[c]:
                sel.calc_all()


    def filter_data(self):
        """
        Filter everything.
        """
        pass