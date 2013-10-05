from __future__ import print_function
from sys import stdout, stderr
import os.path
import string
import re
from itertools import izip_longest
from array import array


__all__ = ["check_fastq_extension", "fastq_quality_str2list", "fastq_quality_list2str", "filter_fastq_chastity", 
           "parse_fastq_header", "print_fastq", "read_fastq", 
           "read_fastq_multi", "reverse_fastq", "trim_fastq", "trim_fastq_length"]


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


def check_fastq_extension(fname, quiet=False, error_handle=stderr):
    """
    Check for FASTQ file extension.

    Returns True if the extension is recognized (.fastq or .fq), else False.
    Prints an error message to error_handle, suppressed if quiet is True.
    """
    ext = os.path.splitext(fname)[-1].lower()
    if ext in (".fq", ".fastq"):
        return True
    else:
        if len(ext) > 0:
            if not quiet:
                print("Warning: unrecognized FASTQ file extension '%s'" % ext,
                      file=stderr)
            return False
        else:
            if not quiet:
                print("Warning: FASTQ file has no file extension",
                      file=stderr)
            return False


def read_fastq(fname, filter_function=None, buffer_size=BUFFER_SIZE):
    """
    Generator function for reading from a FASTQ file.

    Yields the FASTQ record's name, sequence, and quality values as a tuple.

    The filter_function must operate on a FASTQ tuple and return True (pass)
    or False (fail). FASTQ records that fail filtering will be skipped, so
    this feature should not be used when reading files in parallel. Use 
    read_fastq_multi instead.
    """
    try:
        handle = open(fname, "U")
    except IOError:
        print("Error: could not open FASTQ file '%s'" % fname, file=stderr)
        return

    check_fastq_extension(fname, quiet=False, error_handle=stderr)

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
            # (header, sequence, quality)
            fq = (lines[i * 4], lines[i * 4 + 1], lines[i * 4 + 3])
            if filter_function is None: # no filtering
                yield fq
            elif filter_function(fq):   # passes filtering
                yield fq
            else:                       # fails filtering
                continue

    handle.close()


def read_fastq_multi(fnames, filter_function=None, buffer_size=BUFFER_SIZE,
                     match_lengths=True):
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
                             buffer_size=BUFFER_SIZE))

    for records in izip_longest(*fq_generators, fillvalue=None):
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


def print_fastq(fq, file=stdout):
    """
    Prints a FASTQ record to the specified file handle.

    Information in the '+' header is lost.
    """
    print(fq[0], fq[1], "+", fq[2], sep="\n", file=file)


def trim_fastq(fq, start, end):
    """
    Returns a new FASTQ tuple containing bases between start and end.

    Bases are numbered starting at 1. Start and end are inclusive. Header 
    information is unchanged.
    """
    new_sequence = fq[1][start - 1:end]
    new_quality = fq[2][start - 1:end]
    return fq[0], new_sequence, new_quality


def trim_fastq_length(fq, start, length):
    """
    Returns a new FASTQ tuple containing length bases including start.

    Bases are numbered starting at 1. Header information is unchanged.
    """
    end = start + length - 1
    return trim_fastq(fq, start, end)


def reverse_fastq(fq):
    """
    Returns a new FASTQ tuple with reverse-complemented sequence.

    The sequence is reverse-complemented and the quality values are reversed.
    Header information is unchanged.
    """
    new_sequence = fq[1].translate(dna_trans)[::-1]
    new_quality = fq[2][::-1]
    return fq[0], new_sequence, new_quality


def parse_fastq_header(header, pattern=header_pattern):
    """
    Parses a FASTQ header and returns a dictionary.

    Dictionary keys are based on named groups in the regular expression 
    pattern. Unnamed matches are ignored. Integer dictionary values are 
    converted.
    """
    match = pattern.match(header)
    if match is None:
        return None
    else:
        header_dict = match.groupdict()
        for key in header_dict:
            if header_dict[key].isdigit():
                header_dict[key] = int(header_dict[key])
        return header_dict


def filter_fastq_chastity(fq):
    """
    Filtering function for read_fastq and read_fastq_multi.

    Returns True if the chastity bit is set in the FASTQ header.
    """
    match = header_pattern.match(fq[0])
    if match is None:
        return False
    else:
        chastity = match.groupdict()['Chastity']
        if chastity == "1":
            return True
        else:
            return False


def fastq_quality(fq, base=33):
    """
    Convert the quality string to a list of integers.

    base is 33 for Sanger and Illumina 1.8, or 64 for Illumina 1.3 and 1.5.
    """
    quality = [x - base for x in array('b', fq[2]).tolist()]
    return quality


def fastq_quality_reconvert(quality, base=33):
    """
    Convert a list of integers into a quality string.

    base is 33 for Sanger and Illumina 1.8, or 64 for Illumina 1.3 and 1.5.
    """
    fq_quality = [x + base for x in quality]
    return array('b', fq_quality).tostring()


def fastq_min_quality(fq, base=33):
    """
    Return the minimum base quality in the read.
    """
    return min(fastq_quality(fq, base))


def fastq_avg_quality(fq, base=33):
    """
    Return the average base quality in the read.
    """
    return float(sum(fastq_quality(fq))) / len(fq[1])


