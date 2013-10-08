from seqlib import SeqLib
from enrich_error import EnrichError
from fastq_util import *


class BarcodeSeqLib(SeqLib):
    def __init__(self, config):
        SeqLib.__init__(self, config)
        try:
            self.barcode_map_file = config['barcodes']['map file']
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key)

        self.barcode_map = dict()
        self.variant_barcodes = dict()
        self.read_barcode_map(self.barcode_map_file)

        self.set_filters(config, {'max mutations' : len(self.wt_dna)})

        self.barcode_counts = dict()


    def count(self):
        # flags for verbose output of filtered reads
        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        for fq in read_fastq(self.reads):
            if self.reverse_reads:
                fq = reverse_fastq(fq)

            for key in filter_flags:
                filter_flags[key] = False

            # filter the read based on specified quality settings
            if self.filters['min quality'] > 0:
                if fastq_min_quality(fq) < self.filters['min quality']:
                    self.filter_stats['min quality'] += 1
                    filter_flags['min quality'] = True
            if self.filters['avg quality'] > 0:
                if fastq_average_quality(fq) < self.filters['avg quality']:
                    self.filter_stats['avg quality'] += 1
                    filter_flags['avg quality'] = True
            if not any(filter_flags.values()): # passed quality filtering
                mutations = self.count_variant(fq[1])
                if mutations is None: # fused read has too many mutations
                    self.filter_stats['max mutations'] += 1
                    filter_flags['max mutations'] = True
            if any(filter_flags.values()):
                self.filter_stats['total'] += 1
                if self.verbose:
                    self.report_filtered_read(fq, filter_flags)


    def read_barcode_map(self):
        try:
            handle = open(fname, "U")
        except IOError:
            raise EnrichError("Could not open barcode map file '%s'" % fname)

        for line in handle:
            # skip comments and whitespace-only lines
            if len(line.strip()) == 0 or line[0] == '#':
                continue

            try:
                barcode, variant = line.strip().split()
            except ValueError:
                raise EnrichError("Unexpected barcode-variant line")

            if not re.match("^[ACGTacgt]+$", barcode):
                raise EnrichError("Barcode DNA sequence contains unexpected "
                                  "characters")
            if not re.match("^[ACGTacgt]+$", variant):
                raise EnrichError("Variant DNA sequence contains unexpected "
                                  "characters")

            barcode = barcode.upper()
            variant = variant.upper()
            if barcode in self.barcode_map:
                if self.barcode_map[barcode] != variant:
                    raise EnrichError("Barcode '%s' assigned to multiple "
                                      "unique variants" % barcode)
            else:
                self.barcode_map[barcode] = variant

            if variant not in self.variant_barcodes:
                self.variant_barcodes[variant] = set()
            self.variant_barcodes[variant].add(barcode)
        handle.close()


    def read_barcodes(self, fname):
        pass


    def count(self):
        pass


    def filter_barcodes(self):
        pass

