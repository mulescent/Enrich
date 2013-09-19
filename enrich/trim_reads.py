from __future__ import print_function
from sys import stderr
import os.path
import argparse
from fastq_util import read_fastq, print_fastq, trim_fastq, reverse_fastq


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq", help="sequence data file in FASTQ format")
    parser.add_argument("-s", "--start", default=1, metavar="N",
                        help="starting base")
    parser.add_argument("-e", "--end", default=None, metavar="M",
                        help="ending base")
    parser.add_argument("-o", "--output", default=".", metavar="DIR",
                        help="output directory")
    parser.add_argument("-r", "--reverse", action="store_true", default=False,
                        help="reverse-complement the output")
    args = parser.parse_args()

    name, ext = os.path.splitext(os.path.basename(args.fastq))
    output_file_name = name + ".trim"
    if args.reverse:
        output_file_name += ".rev"
    output_file_name += ext
    output_file_name = os.path.join(args.output, output_file_name)
    with open(output_file_name, "w") as output_handle:
        for fq in read_fastq(args.fastq):
            fq = trim_fastq(fq, args.start, args.end)
            if args.reverse:
                fq = reverse_fastq(fq)
            print_fastq(fq, file=output_handle)
