from __future__ import print_function
from enrich_error import EnrichError
import selection


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
                self.conditions[cnd['label']] = cnd['selections']
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
