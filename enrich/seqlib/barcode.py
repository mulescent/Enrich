from __future__ import print_function
import re
from seqlib import SeqLib
from enrich_error import EnrichError
from fastq_util import read_fastq, check_fastq

# debugging
from sys import stdout, stderr


class BarcodeMap(dict):
    def __init__(self, mapfile):
        try:
            handle = open(mapfile, "U")
        except IOError:
            raise EnrichError("Could not open barcode map file '%s'" % mapfile, self.name)

        self.filename = mapfile
        for line in handle:
            # skip comments and whitespace-only lines
            if len(line.strip()) == 0 or line[0] == '#':
                continue

            try:
                barcode, variant = line.strip().split()
            except ValueError:
                raise EnrichError("Unexpected barcode-variant line format", self.name)

            if not re.match("^[ACGTacgt]+$", barcode):
                raise EnrichError("BarcoBde DNA sequence contains unexpected "
                                  "characters", self.name)
            if not re.match("^[ACGTNacgtn]+$", variant):
                raise EnrichError("Variant DNA sequence contains unexpected "
                                  "characters", self.name)

            barcode = barcode.upper()
            variant = variant.upper()
            if barcode in self:
                if self[barcode] != variant:
                    raise EnrichError("Barcode '%s' assigned to multiple "
                                      "unique variants" % barcode, self.name)
            else:
                self[barcode] = variant
        handle.close()

        self.set_variants()


    def set_variants(self):
        """
        Create a reverse-dictionary with variants as keys.

        Variants are keys and values are a list of barcodes associated with 
        that variant. Creation of this dictionary is optional for memory 
        reasons.
        """
        self.variants = dict()
        for bc, v in self.iteritems():
            if v not in self.variants:
                self.variants[v] = list()
            self.variants[v].append(bc)


    def del_variants(self):
        """
        Remove reference to the variants dictionary to save memory.
        """
        self.variants = None



class BarcodeSeqLib(SeqLib):
    def __init__(self, config, barcode_map=None):
        SeqLib.__init__(self, config)
        try:
            if 'map file' in config['barcodes']:
                self.barcode_map = BarcodeMap(config['barcodes']['map file'])
            else:
                self.barcode_map = None

            if 'forward' in config['fastq'] and 'reverse' in config['fastq']:
                raise EnrichError("Multiple FASTQ files specified", self.name)
            elif 'forward' in config['fastq']:
                self.reads = config['fastq']['forward']
                self.reverse_reads = False
            elif 'reverse' in config['fastq']:
                self.reads = config['fastq']['reverse']
                self.reverse_reads = True
            else:
                raise KeyError("'forward' or 'reverse'")

            if 'start' in config['fastq']:
                self.bc_start = config['fastq']['start']
            else:
                self.bc_start = 1
            if 'length' in config['fastq']:
                self.bc_length = config['fastq']['length']
            else:
                self.bc_length = 2147483647 # longer than any read... for now

            self.set_filters(config, {'min quality' : 0,
                                      'avg quality' : 0,
                                      'chastity' : False,
                                      'max mutations' : len(self.wt_dna)})
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, self.name)

        try:
            check_fastq(self.reads)
        except IOError as fqerr:
            raise EnrichError("FASTQ file error: %s" % fqerr, self.name)

        if self.barcode_map is None: # not in local config
            if barcode_map is None:  # not provided on object creation
                raise EnrichError("Barcode map not specified", self.name)
            else:
                self.barcode_map = barcode_map

        self.counters['barcodes'] = dict()


    def count(self):
        # flags for verbose output of filtered reads
        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        # count all the barcodes
        for fq in read_fastq(self.reads):
            fq.trim_length(self.bc_length, start=self.bc_start)
            if self.reverse_reads:
                fq.reverse()

            for key in filter_flags:
                filter_flags[key] = False

            # filter the barcode based on specified quality settings
            if self.filters['chastity']:
                if not fq.is_chaste():
                    self.filter_stats['chastity'] += 1
                    filter_flags['chastity'] = True
            if self.filters['min quality'] > 0:
                if fq.min_quality() < self.filters['min quality']:
                    self.filter_stats['min quality'] += 1
                    filter_flags['min quality'] = True
            if self.filters['avg quality'] > 0:
                if fq.mean_quality() < self.filters['avg quality']:
                    self.filter_stats['avg quality'] += 1
                    filter_flags['avg quality'] = True
            if any(filter_flags.values()): # failed quality filtering
                self.filter_stats['total'] += 1
                if self.verbose:
                    self.report_filtered_read(fq, filter_flags)
            else: # passed quality filtering
                try:
                    self.counters['barcodes'][fq.sequence.upper()] += 1
                except KeyError:
                    self.counters['barcodes'][fq.sequence.upper()] = 1

        # count variants associated with the barcodes
        for bc, count in self.counters['barcodes'].iteritems():
            if bc in self.barcode_map:
                variant = self.barcode_map[bc]
                mutations = self.count_variant(variant, copies=count)
                if mutations is None: # variant has too many mutations
                    self.filter_stats['max mutations'] += 1
                    self.filter_stats['total'] += count # total counts reads
                    if self.verbose:
                        self.report_filtered_variant(variant, count)

        if self.verbose:
            self.report_filter_stats(self.log)
        self.report_filter_stats(stderr)
        self.initialize_df()


    def orphan_barcodes(self, mincount=0):
        orphans = dict()
        for bc, count in self.counters['barcodes'].iteritems():
            if count > mincount:
                if bc not in self.barcode_map:
                    orphans[bc] = count
        return orphans


    def report_filtered_variant(self, variant, count):
        print("Filtered variant (%s)" % \
                    (SeqLib._filter_messages['max mutations']), file=self.log)
        print(variant, file=self.log)
        print("quantity=", count, sep="", file=self.log)

