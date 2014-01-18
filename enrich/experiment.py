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
        self.df_dict = dict()
        self.use_scores = True

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

        for dtype in self.conditions.values()[0].df_dict:
            if all(dtype in x.df_dict for x in self.conditions.values()):
                self.df_dict[dtype] = True
        if len(self.df_dict.keys()) == 0:
            raise EnrichError("No enrichment data present across all selections", 
                              self.name)

        for cnd in self.conditions:
            if len(cnd.timepoints) == 2:
                self.use_scores = False


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
        first = True
        for c in self.conditions:
            s_id = 1
            for s in self.conditions[c]:
                s_label = "%s.%d" % (c, s_id)
                s_id += 1
                s.calc_all()
                if self.use_scores: # keep the score and r_sq columns
                    if first:
                        for dtype in self.df_dict:
                            self.df_dict[dtype] = s.df_dict[dtype][['score', 'r_sq']]
                            cnames = ["%s.%s" % (x, s_label) for x in ['score', 'r_sq']]
                        first = False
                    else:
                        for dtype in self.df_dict:
                            self.df_dict[dtype] = self.df_dict[dtype].join(s.df_dict[dtype][['score', 'r_sq']],
                                how="outer", rsuffix="%s" % s_label)
                            cnames += ["%s.%s" % (x, s_label) for x in ['score', 'r_sq']]
                else:               # only two timepoints, so keep the ratio
                    for dtype in self.df_dict:
                            self.df_dict[dtype] = s.df_dict[dtype][['ratio.%d' % s.timepoints[1]]]
                            cnames = ["ratio.%s" % s_label]
                        first = False
                    else:
                        for dtype in self.df_dict:
                            self.df_dict[dtype] = self.df_dict[dtype].join(s.df_dict[dtype][['ratio.%d' % s.timepoints[1]]],
                                how="outer", rsuffix="%s" % s_label)
                            cnames.append("ratio.%s" % s_label)
                s.save_data(clear=True)
        for dtype in self.df_dict:
            self.df_dict[dtype].columns = cnames


    def filter_data(self):
        """
        Filter everything.
        """
        pass

