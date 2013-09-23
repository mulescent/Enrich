from __future__ import print_function
from sys import stdout, stderr
import os.path
import string
import re


__BUFFER_SIZE = 100000
__dna_trans = string.maketrans('actgACTG', 'tgacTGAC')


# @<MachineName>:<Lane>:<Tile>:<X>:<Y>:<Chastity>#<IndexRead>/<ReadNumber>
header_pattern = re.compile("@(?P<MachineName>.+)"
                            ":(?P<Lane>\d+)"
                            ":(?P<Tile>\d+)"
                            ":(?P<X>\d+)"
                            ":(?P<Y>\d+)"
                            ":(?P<Chastity>[01])"
                            "#(?P<IndexRead>\d)"
                            "/(?P<ReadNumber>\d)")


def read_fastq(fname, filter_function=None, buffer_size=__BUFFER_SIZE):
    """Generator function for reading records from a FASTQ file.

    Yields the FASTQ record's name, sequence, and quality values.
    """
    try:
        handle = open(fname, 'U')
    except IOError:
        print("Error: could not open FASTQ file '%s'" % fname, file=stderr)
        return

    ext = os.path.splitext(fname)[-1].lower()
    if ext not in ('.fq', '.fastq'):
        if len(ext) > 0:
            print("Warning: unrecognized FASTQ file extension '%s'" % ext, 
                  file=stderr)
        else:
            print("Warning: FASTQ file has no file extension", file=stderr)

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


def print_fastq(fq, handle=stdout):
    print(fq[0], fq[1], "+", fq[2], sep="\n", file=handle)


def trim_fastq(fq, start, end):
    new_sequence = fq[1][start - 1:end]
    new_quality = fq[2][start - 1:end]
    return fq[0], new_sequence, new_quality


def reverse_fastq(fq):
    new_sequence = fq[1].translate(__dna_trans)[::-1]
    new_quality = fq[2][::-1]
    return fq[0], new_sequence, new_quality


def parse_header(header, pattern=header_pattern):
    match = pattern.match(header)
    if match is None:
        return None
    else:
        header_dict = match.groupdict()
        for key in header_dict:
            if header_dict[key].isdigit():
                header_dict[key] = int(header_dict[key])
        return header_dict
