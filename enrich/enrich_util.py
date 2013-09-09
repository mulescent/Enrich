#!/usr/bin/env python
'''
enrich_util: this module contains a number of small utility functions, and is intended to be extensible
'''

__author__ = "Douglas M. Fowler"
__copyright__ = "Copyright 2011"
__credits__ = ["Douglas M Fowler", "Carlos L. Araya"]
__license__ = "FreeBSD"
__version__ = "0.2"
__maintainer__ = "Douglas M. Fowler"
__email__ = "dfowler@uw.edu"

def valuesort(input_dict):
    '''Returns the keys of dictionary sorted by their values'''

    items=input_dict.items()
    backitems=[[v[1],v[0]] for v in items]
    backitems.sort()
    return [ backitems[i][1] for i in range(0,len(backitems)) ]
    
def build_tally_dict(input_file):
    '''define a function to build a tally of sequences observed from filtered reads'''
    properties_dict = {}
    tally_dict = {}
    input_key = open(input_file)
    header = input_key.readline()
    input_line = input_key.readline()
    while input_line:
        readID, sequence, match_count, mutation_count, mutation_location, mutation_identity, max_mutation_run = input_line.rstrip('\n').split('\t')
        if not mutation_location in tally_dict:
            tally_dict[mutation_location] = {}
            properties_dict[mutation_location] = {}
        if mutation_identity in tally_dict[mutation_location]:
            tally_dict[mutation_location][mutation_identity] += 1
        else:
            tally_dict[mutation_location][mutation_identity] = 1
            properties_dict[mutation_location][mutation_identity] = [sequence, match_count, mutation_count, mutation_location, mutation_identity, max_mutation_run]
        input_line = input_key.readline()   
    return tally_dict, properties_dict

def norm_count_dict(count_dict):
    '''define a function that normalizes the counts in a counts dictionary'''

    norm_dict = {}
    total = 0
    for seqID in count_dict:
        total += count_dict[seqID]
    for seqID in count_dict:
        norm_dict[seqID] = float(count_dict[seqID])/total
    return norm_dict
    

def build(input_file):
    '''define a function to build a dictionary of input lines'''

    output_dict = {}
    input_key = open(input_file)
    header = input_key.readline()
    input_line = input_key.readline()
    while input_line:
        input_items = input_line.rstrip('\n').split('\t')
        lineID = input_items[0]
        if not lineID in output_dict:
            output_dict[lineID] = input_items[1:]
        else:
            print "ERROR: ID is not unique to line in file"
        input_line = input_key.readline()
    return output_dict

def build_value_dict(input_file, id_index, value_index, mode):
    '''define a function that builds a dictionary representation of a column of data'''

    value_dict = {}
    input_key = open(input_file)
    header = input_key.readline()
    input_line = input_key.readline()
    
    if mode == "str":
        def format(x):
            return x
    elif mode == "float":
        def format(x):
            return float(x)     
    elif mode == "int":
        def format(x):
            return int(x)   
            
    while input_line:
        input_items = input_line.rstrip('\n').split('\t')
        value_dict[input_items[id_index]] = format(input_items[value_index])
        input_line = input_key.readline()
    return value_dict
    
def gen_ratio_dict(A_dict, B_dict):
    '''define a function that returns the ratio between two frequency dictionaries'''

    ratio_dict = {}
    for seqID in set(A_dict).intersection(set(B_dict)):
        ratio_dict[seqID] = float(A_dict[seqID])/float(B_dict[seqID])
    return ratio_dict
