from __future__ import print_function
from enrich_error import EnrichError
import selection


class Experiment(object):
    def __init__(self, config):
        self.name = "Unnamed" + self.__class__.__name__
        self.verbose = False
        self.log = None

        try:
            self.name = config['name']
            pass
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
