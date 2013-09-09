#!/usr/bin/env python
'''
map_counts: The map_counts module identifies unique variants in read_aligner output files and produces an output file that contains a list of those unique sequences as well as the number of times they appear normalized by the total number of counts in the input file (i.e. the frequency).
'''

import sys, time, optparse # import general libraries
import enrich_util # import project-specific libraries

__author__ = "Douglas M. Fowler"
__copyright__ = "Copyright 2011"
__credits__ = ["Douglas M Fowler", "Carlos L. Araya"]
__license__ = "FreeBSD"
__version__ = "0.2"
__maintainer__ = "Douglas M. Fowler"
__email__ = "dfowler@uw.edu"

def main(project_directory, input_file, grid = 'L'):

    if grid != 'L': #print logging information, if being called from the cluster
        print time.asctime(time.localtime())
        print project_directory, input_file
    
    if "_DNA_" in input_file:
        filter_chars = "N"
    elif "_PRO_" in input_file:
        filter_chars = "X"
    
    try:
        #build a dictionary of the sequences observed and counts:
        tally_dict, properties_dict = enrich_util.build_tally_dict(project_directory + 'data/tmp/' + input_file)
                
        #build a dictionary of the sequence codes and counts:
        count_dict = {}
        for mutation_location in sorted(tally_dict):
            for mutation_identity in sorted(tally_dict[mutation_location]):
                seqID = mutation_location + '-' + mutation_identity
                count_dict[seqID] = tally_dict[mutation_location][mutation_identity]
        seqIDs = enrich_util.valuesort(count_dict)
        
        #build a dictionary of the normalized tally:
        norm_dict = enrich_util.norm_count_dict(count_dict)

    except:
        print 'Error: could not build dictionary of counts from input file'
        return(1)
        
    try:
        #print sequences and counts:
        f = open(project_directory + 'data/output/' + 'counts_' + input_file,'w')
        f_1 = open(project_directory + 'data/output/' + 'counts_' + input_file + '.m1', 'w')
        f_2 = open(project_directory + 'data/output/' + 'counts_' + input_file + '.m2', 'w')

        print >>f, '\t'.join(["seqID","sequence","match_count","mutation_count","mutation_location","mutation_identity","max_mutation_run","sequence_frequency","sequence_count"] )
        print >>f_1, '\t'.join(["seqID","sequence","match_count","mutation_count","mutation_location","mutation_identity","max_mutation_run","sequence_frequency","sequence_count"] )
        print >>f_2, '\t'.join(["seqID","sequence","match_count","mutation_count","mutation_location","mutation_identity","max_mutation_run","sequence_frequency","sequence_count"] )
        
        seqIDs.reverse()
        for seqID in seqIDs:
            mutation_location, mutation_identity = seqID.split('-')
            norm = norm_dict[seqID]
            tally = tally_dict[mutation_location][mutation_identity]
            sequence = properties_dict[mutation_location][mutation_identity][0]
            if not filter_chars in sequence:
                print >>f, '\t'.join([seqID] + properties_dict[mutation_location][mutation_identity] + [str(norm), str(tally)])
             
                if len(mutation_location.split(',')) == 1 and mutation_location != 'NA':
                    print >>f_1, '\t'.join([seqID] + properties_dict[mutation_location][mutation_identity] + [str(norm), str(tally)])     
                
                if len(mutation_location.split(',')) == 2:
                    print >>f_2, '\t'.join([seqID] + properties_dict[mutation_location][mutation_identity] + [str(norm), str(tally)])  
                          
        f.close()
        f_1.close()
        f_2.close()
        
    except:
        print 'Error: write to output file'
        return(1)

    return(0)
            
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('--path', action = 'store', type = 'string', dest = 'path', help = 'exact location from / project directory')
    parser.add_option('--infile', action = 'store', type = 'string', dest = 'infile', help = 'input filename')
    (option, args) = parser.parse_args()
    parser.add_option('--local', action = 'store', type = 'string', dest = 'local', help = 'Is this a local (L) run or should an SGE (SGE) job be scheduled?')
    main(option.path, option.infile, option.local)
