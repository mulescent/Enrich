from __future__ import print_function
from sys import stderr
import os.path
import re
import string
import itertools
from array import array


__all__ = ["FQRead", "check_fastq", "read_fastq", "read_fastq_multi"]


# Matches FASTQ headers based on the following pattern (modify as needed):
# @<MachineName>:<Lane>:<Tile>:<X>:<Y>:<Chastity>#<IndexRead>/<ReadNumber>
header_pattern = re.compile("@(?P<MachineName>.+)"
                            ":(?P<Lane>\d+)"
                            ":(?P<Tile>\d+)"
                            ":(?P<X>\d+)"
                            ":(?P<Y>\d+)"
                            ":(?P<Chastity>[01])"
                            "#(?P<IndexRead>\d)"
                            "/(?P<ReadNumber>\d)")


BUFFER_SIZE = 100000 # empirically optimized for reading FASTQ files


dna_trans = string.maketrans("actgACTG", "tgacTGAC")


class FQRead(object):
    # no instance dictionary, more memory efficient
    __slots__ = ('header', 'sequence', 'header2', 'quality', 'qbase')


    def __init__(self, header, sequence, header2, quality, qbase=33):
        if len(sequence) != len(quality):
            raise ValueError('different lengths for sequence and quality')
        else:
            self.header = header
            self.sequence = sequence
            self.header2 = header2
            # quality is a list of integers
            self.quality = [x - qbase for x in array('b', quality).tolist()]
            self.qbase = qbase


    def __str__(self):
        """
        Reformat as a FASTQ read.
        """
        return '\n'.join([self.header, self.sequence, self.header2, 
                array('b', [x + self.qbase for x in self.quality]).tostring()])


    def __len__(self):
        """
        Length of the sequence.
        """
        return len(self.sequence)


    def trim(self, start=1, end=-1):
        """
        Trim this read to contain bases between start and end.

        Bases are numbered starting at 1. Start and end are inclusive. Header 
        information is unchanged.
        """
        self.sequence = self.sequence[start - 1:end]
        self.quality = self.quality[start - 1:end]


    def trim_length(fq, length, start=1):
        """
        Trim this read to contain length bases including start.

        Bases are numbered starting at 1. Header information is unchanged.
        """
        self.trim(start=start, end=start + length - 1)


    def reverse(self):
        """
        Reverse-complement the sequence in place.

        The sequence is reverse-complemented and the quality values are reversed.
        Header information is unchanged.
        """
        self.sequence = self.sequence.translate(dna_trans)[::-1]
        self.quality = self.quality[::-1]


    def header_information(self, pattern=header_pattern):
        """
        Parses the first FASTQ header and returns a dictionary.

        Dictionary keys are based on named groups in the regular expression 
        pattern. Unnamed matches are ignored. Integer dictionary values are 
        converted.
        """
        match = pattern.match(self.header)
        if match is None:
            return None
        else:
            header_dict = match.groupdict()
            for key in header_dict:
                if header_dict[key].isdigit():
                    header_dict[key] = int(header_dict[key])
            return header_dict


    def min_quality(self):
        """
        Return the minimum base quality in the read.
        """
        return min(self.quality)


    def mean_quality(self):
        """
        Return the average base quality in the read.
        """
        return float(sum(self.quality)) / len(self)


    def is_chaste(self):
        """
        Returns True if the chastity bit is set in the FASTQ header.
        """
        try:
            if self.header_information()['Chastity'] == "1":
                return True
            else:
                return False
        except KeyError, TypeError: # no 'Chastity' in pattern or no match
            return False



def check_fastq(fname):
    """
    Check that the file exists and has a valid FASTQ file extension.

    Returns if the file exists and the extension is recognized (.fastq 
    or .fq), otherwise raise an IOError.
    """
    if os.path.isfile(fname):
        ext = os.path.splitext(fname)[-1].lower()
        if ext in (".fq", ".fastq"):
            return None
        else:
            raise IOError("improper file extension for '%s'" % fname)
    else:
        raise IOError("file '%s' doesn't exist" % fname)



def read_fastq(fname, filter_function=None, buffer_size=BUFFER_SIZE, qbase=33):
    """
    Generator function for reading from a FASTQ file.

    Yields the FASTQ record's name, sequence, and quality values as a tuple.

    The filter_function must operate on a FASTQ tuple and return True (pass)
    or False (fail). FASTQ records that fail filtering will be skipped, so
    this feature should not be used when reading files in parallel. Use 
    read_fastq_multi filtering instead.
    """
    check_fastq(fname)
    handle = open(fname, "U")

    eof = False
    leftover = ''

    while not eof:
        buf = handle.read(buffer_size)
        if len(buf) < buffer_size:
            eof = True

        buf = leftover + buf # prepend partial record from previous buffer
        lines = buf.split('\n')
        fastq_count = len(lines) / 4

        if not eof: # handle lines from the trailing partial FASTQ record
            dangling = len(lines) % 4
            if dangling == 0: # quality line (probably) incomplete
                dangling = 4
                fastq_count = fastq_count - 1
            # join the leftover lines back into a string
            leftover = '\n'.join(lines[len(lines) - dangling:])

        # index into the list of lines to pull out the FASTQ records
        for i in xrange(fastq_count):
            # (header, sequence, header2, quality)
            fq = FQRead(*lines[i * 4:(i + 1) * 4], qbase=qbase)
            if filter_function is None: # no filtering
                yield fq
            elif filter_function(fq):   # passes filtering
                yield fq
            else:                       # fails filtering
                continue

    handle.close()



def read_fastq_multi(fnames, filter_function=None, buffer_size=BUFFER_SIZE,
                     match_lengths=True, qbase=33):
    """
    Generator function for reading from multiple FASTQ files in parallel.

    Yields the a tuple of tuples containing the FASTQ record's name, 
    sequence, and quality values.

    fnames must be an iterable of FASTQ file names.

    The filter_function must operate on a FASTQ tuple and return True (pass)
    or False (fail). If any of the FASTQ records read in a given set fails
    filtering, the set will be skipped (not returned). Note that the filtering
    is performed on each FASTQ record independently.

    If match_lengths is True, the generator will yield None if the files
    do not contain the same number of FASTQ records. Otherwise, it will 
    silently ignore partial records.
    """
    fq_generators = list()
    for f in fnames:
        fq_generators.append(read_fastq(f, filter_function=None,
                             buffer_size=BUFFER_SIZE, qbase=qbase))

    for records in itertools.izip_longest(*fq_generators, fillvalue=None):
        if None in records: # mismatched file lengths
            if match_lengths:
                yield None
            else:
                break # shortest FASTQ file is empty, so we're done
        if filter_function is None:                     # no filtering
            yield records
        elif all(filter_function(x) for x in records):  # pass filtering
            yield records
        else:                                           # fail filtering
            continue



def fastq_filter_chastity(fq):
    """
    Filtering function for read_fastq and read_fastq_multi.

    Returns True if the chastity bit is set in the FASTQ header.
    """
    return fq.is_chaste()


