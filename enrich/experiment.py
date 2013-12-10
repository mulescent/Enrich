from __future__ import print_function
from enrich_error import EnrichError
import selection


class Experiment(object):
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
        self.verbose = True
        self.log = log
        try:
            print("# Logging started for '%s': %s" % 
                        (self.name, time.asctime()), file=self.log)
        except (IOError, AttributeError):
            raise EnrichException("Could not write to log file", self.name)
