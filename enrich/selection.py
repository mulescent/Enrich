from enrich_error import EnrichError
import seqlib.basic
import seqlib.barcode
import seqlib.overlap

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

        if all(lib.libtype == "barcode" for lib in self.libraries):
            self.barcode_enrichments = dict()
        self.enrichments = dict()
        self.enrichments_nt = dict()
        self.enrichments_aa = dict()


    def count(self):
        for lib in self.libraries:
            lib.count()
            lib.count_mutations()


    def calc_enrichments(self):
        if len(self.libraries) > 2:
            # regression
        elif len(self.libraries) == 2:
            # ratio
        else:
            pass


