from __future__ import print_function
from sys import stderr
import os.path
import argparse
import json
import itertools
from fastq_util import read_fastq, print_fastq


def assign_library_ids(config):
    """Assign unique identifiers to each libary in the config object.

    Checks user-specified libary identifiers for uniqueness. If nonunique
    identifiers are encountered, returns None.

    Returns config.
    """
    id_list = []

    # build list of all user-specified library identifiers
    duplicate_error = False
    for lib in config['libraries']:
        if 'id' in lib:
            lib['id'] = str(lib['id'])
            if lib['id'] in id_list:
                print("Error: duplicate library identifier '%s'" % lib['id'],
                      file=stderr)
                duplicate_error = True
            else:
                id_list.append(lib['id'])
    if duplicate_error:
        return None

    # assign identifiers as needed
    current_id = 1
    for lib in config['libraries']:
        if 'id' not in lib:
            while str(current_id) in id_list:
                current_id += 1
            lib['id'] = str(current_id)
            current_id += 1

    return config


def assign_library_mismatch_thresholds(config, mismatch_threshold=None):
    unassigned_error = False
    noninteger_error = False
    for lib in config['libraries']:
        if 'mismatch_threshold' not in lib['index']:
            if mismatch_threshold is None:
                print("Error: no mismatch_threshold specified for '%s'" %
                    lib['name'], file=stderr)
                unassigned_error = True
            else:
                lib['index']['mismatch_threshold'] = mismatch_threshold
        else:
            try:
                lib['index']['mismatch_threshold'] = \
                    int(lib['index']['mismatch_threshold'])
            except ValueError:
                print("Error: non-integer mismatch_threshold for '%s'" %
                    lib['name'], file=stderr)
                noninteger_error = True

    if unassigned_error or noninteger_error:
        return None
    else:
        return config


def split_fastq(config, outdir, index, forward, reverse):
    """
    """
    if index is None:
        print("Error: no index file specified for split_fastq", file=stderr)
        return

    # build an iterator to process the files in parallel
    fq_handles = {} # output file handles and index read sequences
    if forward is not None and reverse is not None:
        fq_iterator = itertools.izip_longest(read_fastq(index), 
                                             read_fastq(forward),
                                             read_fastq(reverse),
                                             fillvalue=None)

        for library in config['libraries']:
            name, ext = os.path.splitext(os.path.basename(index))
            index_name = name + ".id%s" % library['id'] + ext
            index_name = os.path.join(outdir, index_name)

            name, ext = os.path.splitext(os.path.basename(forward))
            forward_name = name + ".id%s" % library['id'] + ext
            forward_name = os.path.join(outdir, forward_name)

            name, ext = os.path.splitext(os.path.basename(reverse))
            reverse_name = name + ".id%s" % library['id'] + ext
            reverse_name = os.path.join(outdir, reverse_name)

            fq_handles[library['id']] = \
                (open(index_name, "w"), open(forward_name, "w"),
                 open(reverse_name, "w"))

    elif forward is not None:
        fq_iterator = itertools.izip_longest(read_fastq(index), 
                                             read_fastq(forward),
                                             fillvalue=None)

        for library in config['libraries']:
            name, ext = os.path.splitext(os.path.basename(index))
            index_name = name + ".id%s" % library['id'] + ext
            index_name = os.path.join(outdir, index_name)

            name, ext = os.path.splitext(os.path.basename(forward))
            forward_name = name + ".id%s" % library['id'] + ext
            forward_name = os.path.join(outdir, forward_name)

            fq_handles[library['id']] = \
                (open(index_name, "w"), open(forward_name, "w"))

    elif reverse is not None:
        fq_iterator = itertools.izip_longest(read_fastq(index),
                                             read_fastq(reverse),
                                             fillvalue=None)

        for library in config['libraries']:
            name, ext = os.path.splitext(os.path.basename(index))
            index_name = name + ".id%s" % library['id'] + ext
            index_name = os.path.join(outdir, index_name)

            name, ext = os.path.splitext(os.path.basename(reverse))
            reverse_name = name + ".id%s" % library['id'] + ext
            reverse_name = os.path.join(outdir, reverse_name)

            fq_handles[library['id']] = \
                (open(index_name, "w"), open(reverse_name, "w"))

    else:
        print("Warning: no forward or reverse files specfied for split_fastq",
              file=stderr)
        return

    for t in fq_iterator:
        if None in t:
            print("Error: FASTQ files are not the same length", file=stderr)
            break
        index_read = t[0][1]
        library_match = None
        for library in config['libraries']:
            mismatches = 0
            for i in xrange(len(library['index']['sequence'])):
                if index_read[i] != library['index']['sequence'][i]:
                    mismatches += 1
                    if mismatches > library['index']['mismatch_threshold']:
                        break
            if mismatches <= library['index']['mismatch_threshold']:
                library_match = library['id']
                break

        if library_match:
            for i in xrange(len(t)):
                print_fastq(t[i], file=fq_handles[library_match][i])

    # close all the files
    for handle_tuple in fq_handles.values():
        for h in handle_tuple:
            h.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="configuration file in JSON format")
    parser.add_argument("-f", "--forward", metavar="FQ", 
                        help="forward read FASTQ file")
    parser.add_argument("-r", "--reverse", metavar="FQ", 
                        help="reverse read FASTQ file")
    parser.add_argument("-i", "--index", metavar="FQ", 
                        help="index read FASTQ file")
    parser.add_argument("-o", "--output", default=".", metavar="DIR",
                        help="output directory")
    parser.add_argument("-m", "--mismatches", metavar="N",
                        help="index read mismatch threshold")
    args = parser.parse_args()

    config = assign_library_ids(json.load(open(args.config)))
    config = assign_library_mismatch_thresholds(config, 
                mismatch_threshold=args.mismatches)

    split_fastq(config, args.output, args.index, args.forward, args.reverse)
    
    name, ext = os.path.splitext(os.path.basename(args.config))
    postrun_config_name = name + ".postrun" + ext
    postrun_config_name = os.path.join(args.output, postrun_config_name)
    with open(postrun_config_name, "w") as config_handle:
        json.dump(config, config_handle, indent=2, sort_keys=True)
