from __future__ import print_function
from enrich_error import EnrichError
import selection


class Experiment(object):
    def __init__(self, config):
        try:
            pass
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key)


    