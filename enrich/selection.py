from __future__ import print_function
from enrich_error import EnrichError
from scipy import stats
import seqlib.basic
import seqlib.barcode
import seqlib.overlap
import os.path
import math


def new_seqlib(config):
    if 'barcodes' in config:
        return seqlib.barcode.BarcodeSeqLib(config)
    elif 'overlap' in config:
        return seqlib.overlap.OverlapSeqLib(config)
    else:
        return seqlib.basic.BasicSeqLib(config)


class Selection(object):
    def __init__(self, config):
        self.libraries = list()
        try:
            for lib in config['libraries']:
                self.libraries.append(new_seqlib(lib))
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key)

        if len(self.libraries) == 0:
            raise EnrichError("Selection class has no libraries")
        self.check_wt()

        if 'carryover' in config:
            if config['carryover'] in ("stop"):
                self.carryover = config['carryover']
            else:
                self.carrover = None
                raise EnrichError("Invalid carryover correction option \"%s\"" % config['carryover'])
        else:
            self.carryover = None

        self.enrichments = None
        self.enrichments_nt = None
        self.enrichments_aa = None
        self.barcode_enrichments = None


    def check_wt(self):
        try:
            coding = [lib['wild type']['coding'] for lib in self.libraries]
            sequence = [lib['wild type']['sequence'] for lib in self.libraries]
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key)
        
        if len(set(coding)) != 1 or len(set(sequence)) != 1:
            raise EnrichError("Inconsistent wild type sequences")


    def count(self):
        # NOT FINISHED
        for lib in self.libraries:
            lib.count()
            lib.count_mutations()


    def count_enrichments(self):
        for counter_type, enrichment_type in 
                (('variants', 'enrichments'), 
                 ('mutations_nt', 'enrichments_nt'), 
                 ('mutations_aa', 'enrichments_aa'), 
                 ('barcode_counts', 'barcode_enrichments')):
            seqlib_counters = [getattr(lib, counter_type, None) \
                               for lib in self.libraries]
            if all(c is None for c in seqlib_counters):
                continue
            elif None in seqlib_counters: # absent in only some SeqLibs
                raise EnrichError("Inconsistent SeqLib count data")

            counted_set = set()
            for counts in seqlib_counters:
                counted_set.update(counts.keys())

            enrichment_dict = dict()
            for k in counted_set: # all encountered variants
                enrichment_dict[k] = {'counts' : list(), 
                                      'score' : float("NaN")}
                for counts in seqlib_counters:
                    if k in counts:
                        enrichment_dict[k]['counts'].append(counts[k])
                    else: # not found in this SeqLib
                        enrichment_dict[k]['counts'].append(float("NaN"))

            setattr(self, enrichment_type, enrichment_dict)


    def score_enrichments(self):
        for enrichment_type in ('enrichments', 'enrichments_nt', 
                                'enrichments_aa', 'barcode_enrichments'):
            enrichment_dict = getattr(self, enrichment_type)
            if enrichment_dict is not None:
                for k in enrichment_dict:
                    if math.isnan(enrichment_dict[k]['counts'][0]):
                        # must be present in first library
                        score = float("NaN")
                    else:
                        values = [x for x in enrichment_dict[k]['counts']
                                  if not math.isnan(x)]
                        if len(values) == 1:
                            score = float("NaN")
                        elif len(values) == 2:
                            score = values[1] / float(values[0])
                        else:
                            xs = range(0, len(values))
                            # slope, intercept, r_value, p_value, std_err
                            score, _, _, _, _ = stats.linregress(xs, values)
                    enrichment_dict[k]['score'] = score

        
    def write_enrichments(self, directory):
        for enrichment_type in ('enrichments', 'enrichments_nt', 
                                'enrichments_aa', 'barcode_enrichments'):
            enrichment_dict = getattr(self, enrichment_type)
            if enrichment_dict is not None:
                fname = os.path.join(directory, enrichment_type + ".tsv")
                handle = open(fname, "w")
                header = ["sequence"]
                header += ["count%d" % i for i in xrange(len(self.libraries))]
                header += ["score"]
                print("\t".join(header), file=handle)
                for k in enrichment_dict:
                    print(k, "\t".join(enrichment_dict[k]['counts']), 
                          enrichment_dict[k]['score'])


    def calc_enrichments(self, counter_name, enrichment_name):
        counters = [getattr(x, counter_name, None) for x in self.libraries]
        if all(c is None for c in counters):
            return None
        elif None in counters: # absent in only some SeqLibs
            raise EnrichError("Inconsistent SeqLib count data")

        allkeys = set()
        for c in counters:
            allkeys.update(c.keys())

        edict = dict()
        for k in allkeys:
            values = list()
            for c in counters:
                if k in c:
                    values.append(c[k])
                else:
                    values.append(float("NaN"))
            if math.isnan(values[0]): # must be present in first library
                score = float("NaN")
            else:
                values = [x for x in values if not math.isnan(x)]
                if len(values) == 1: # can't calculate
                    score = float("NaN")
                elif len(values) == 2: # ratio
                    score = values[1] / float(values[0])
                else:
                    xs = range(0, len(values))
                    # slope, intercept, r_value, p_value, std_err
                    score, _, _, _, _ = stats.linregress(xs, values)
            edict[k] = score

        setattr(self, enrichment_name, edict)


    def filter_barcodes(self):
        # filter barcodes on consistency
        pass


