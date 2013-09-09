#!/usr/bin/env python
'''
read_aligner: the read_aligner module translates each finished sequence and then aligns both the DNA and protein sequences to the wtDNA and wtPRO configuration elements.  

For the sake of speed, no alignment algorithm is used. Instead a simple base-by-base comparison is made. Read_aligner also filters the output sequences on a number of parameters. The unresolvable_max configuration element specifies the maximum number of unresolvable bases to permit, and should generally be set to 0. The maximum_mutation_run configuration element specifies the maximum number of consecutive DNA mutations to permit.  In a lightly mutagenized library, reads with a large number (>3) of consecutive mutations are generally sequencing errors. The avg_quality configuration element specifies the minimum read-average QPhred quality score to permit, and is generally set to 20. The chaste configuration element specifies whether reads that have failed the Illumina quality chastity filter should be included, and is generally set to 'n'.  The Ncount_max configuration element specifies the maximum number of 'N' bases a read may contain, and is generally set to 0.
'''

import optparse, sys, time, numpy #import standard modules

__author__ = "Douglas M. Fowler"
__copyright__ = "Copyright 2011"
__credits__ = ["Douglas M Fowler", "Carlos L. Araya"]
__license__ = "FreeBSD"
__version__ = "0.2"
__maintainer__ = "Douglas M. Fowler"
__email__ = "dfowler@uw.edu"

def main(path, infile, referenceDNA, referenceAA, gap_max, unresolvable_max, maxmutrun, avg_quality, chaste, Ncount_max, mode, grid = 'L'):
    codon_table = {'TTT':'F', 'TCT':'S', 'TAT':'Y', 'TGT':'C',
    'TTC':'F', 'TCC':'S', 'TAC':'Y', 'TGC':'C',
    'TTA':'L', 'TCA':'S', 'TAA':'*', 'TGA':'*',
    'TTG':'L', 'TCG':'S', 'TAG':'*', 'TGG':'W',
    'CTT':'L', 'CCT':'P', 'CAT':'H', 'CGT':'R',
    'CTC':'L', 'CCC':'P', 'CAC':'H', 'CGC':'R',
    'CTA':'L', 'CCA':'P', 'CAA':'Q', 'CGA':'R',
    'CTG':'L', 'CCG':'P', 'CAG':'Q', 'CGG':'R',
    'ATT':'I', 'ACT':'T', 'AAT':'N', 'AGT':'S',
    'ATC':'I', 'ACC':'T', 'AAC':'N', 'AGC':'S',
    'ATA':'I', 'ACA':'T', 'AAA':'K', 'AGA':'R',
    'ATG':'M', 'ACG':'T', 'AAG':'K', 'AGG':'R',
    'GTT':'V', 'GCT':'A', 'GAT':'D', 'GGT':'G',
    'GTC':'V', 'GCC':'A', 'GAC':'D', 'GGC':'G',
    'GTA':'V', 'GCA':'A', 'GAA':'E', 'GGA':'G',
    'GTG':'V', 'GCG':'A', 'GAG':'E', 'GGG':'G'}

    if grid != 'L': #print log information if running on a cluster
        print time.asctime(time.localtime())
        print path, infile, referenceDNA, referenceAA, gap_max, unresolvable_max, maxmutrun, avg_quality, chaste, Ncount_max, mode, grid

    try:
        gap_max, unresolvable_max, maxmutrun, avg_quality, Ncount_max = map(int, [gap_max, unresolvable_max, maxmutrun, avg_quality, Ncount_max])
    
    except:
        print 'Error: integer-only parameters were not integers'
        return 1

    #open the files for input and output
    try:
       f_infile = open((path + 'data/tmp/' + infile), 'U')
    
    except:
        print 'Error: input files could not be opened'
        return 1
            
    f_DNA_output = open((path + 'data/tmp/' + infile + '_DNA_qc'), 'w')
    f_protein_output = open((path + 'data/tmp/' + infile + '_PRO_qc'), 'w')
    f_aligner_removed_sequences = open((path + 'data/tmp/' + infile + '_qc_removed'), 'w')
    
    print >> f_DNA_output, '\t'.join(['readID', 'sequence', 'match_count', 'mutation_count', 'mutation_location', 'mutation_identity', 'max_mutation_run'])
    print >> f_protein_output, '\t'.join(['readID', 'sequence', 'match_count', 'mutation_count', 'mutation_location', 'mutation_identity', 'max_mutation_run'])
    print >> f_aligner_removed_sequences, 'readID'
    
    line = f_infile.readline() #read the first line and discard it, so that the header is not read
    
    while True:
        line = f_infile.readline().rstrip().split('\t')
        
        #check to see if EOF has arrived
        if len(line[0]) == 0:
            if __name__ == '__main__':
                print 'GRACEFUL EXIT: EOF'
                print time.asctime(time.localtime())
            
            if grid != 'L':
                print 'Job completed successfully at ' + time.asctime(time.localtime())
            
            return 0
                
        if mode == 'B':
            
            read = FuserData(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13])  
            
            #perform quality filtration
            if read.gap_count > gap_max or read.unresolvable_count > unresolvable_max or read.maxmutrun > maxmutrun or read.read1_avgquality < avg_quality or read.read2_avgquality < avg_quality or read.read1_Ncount > Ncount_max or read.read2_Ncount > Ncount_max:
                print >> f_aligner_removed_sequences, str(read.ID)

            elif chaste == 'y' and (read.read1_chastity == '0' or read.read2_chastity == '0' or read.read1_chastity == 'N' or read.read2_chastity == 'N'):
                print >> f_aligner_removed_sequences, str(read.ID)
    
            #ennumerate DNA matches/mismatches
            else:
                for i in xrange(0,len(referenceDNA)):
                    if read.sequence[i] == referenceDNA[i]:
                        read.DNAmatch_count = read.DNAmatch_count + 1
                    else: 
                        read.DNAmutation_count = read.DNAmutation_count + 1
                        read.DNAmutation_location.append(i)
                        read.DNAmutation_identity.append(read.sequence[i])
            
                #translate, and enumerate protein matches/mismatches
                read.AAsequence = ''
                
                for i in xrange(0, len(referenceDNA), 3):
                    read.AAsequence += codon_table[read.sequence[i:i+3]] 
                    
                for i in xrange(0,len(referenceAA)):
                    if read.AAsequence[i] == referenceAA[i]:
                        read.AAmatch_count = read.AAmatch_count + 1
                    else: 
                        read.AAmutation_count = read.AAmutation_count + 1
                        read.AAmutation_location.append(i)
                        read.AAmutation_identity.append(read.AAsequence[i]) 
                
                if read.DNAmutation_count == 0:
                    read.DNAmutation_location = 'NA'
                    read.DNAmutation_identity = 'NA'
                    print >> f_DNA_output, '\t'.join(map(str, [read.ID, read.sequence, read.DNAmatch_count, read.DNAmutation_count, read.DNAmutation_location, read.DNAmutation_identity, read.maxmutrun]))
                
                else:
                    print >> f_DNA_output, '\t'.join(map(str, [read.ID, read.sequence, read.DNAmatch_count, read.DNAmutation_count, ','.join(map(str, read.DNAmutation_location)), ','.join(map(str, read.DNAmutation_identity)), read.maxmutrun]))
                    
                if read.AAmutation_count == 0:
                    read.AAmutation_location = 'NA'
                    read.AAmutation_identity = 'NA'
                    print >> f_protein_output, '\t'.join(map(str, [read.ID, read.AAsequence, read.AAmatch_count, read.AAmutation_count, read.AAmutation_location, read.AAmutation_identity, read.maxmutrun]))
                
                else:
                    print >> f_protein_output, '\t'.join(map(str, [read.ID, read.AAsequence, read.AAmatch_count, read.AAmutation_count, ','.join(map(str, read.AAmutation_location)), ','.join(map(str, read.AAmutation_identity)), read.maxmutrun]))
        
        elif mode == 'R1' or mode == 'R2':
        
            read = FuserData(line[0], line[1], 0, 0, 0, 0, 0, line[7], line[8], line[9], line[10], 0, 0, 0)
            
            if read.maxmutrun > maxmutrun or read.read1_avgquality < avg_quality or read.read1_Ncount > Ncount_max:
                print >> f_aligner_removed_sequences, str(read.ID)
                
            elif chaste == 1 and (read.read1_chastity == '0' or read.read1_chastity == 'N'):
                print >> f_aligner_removed_sequences, str(read.ID)
                
            #ennumerate DNA matches/mismatches
            else:
                for i in range(0,len(referenceDNA)):
                    if read.sequence[i] == referenceDNA[i]:
                        read.DNAmatch_count = read.DNAmatch_count + 1
                    else: 
                        read.DNAmutation_count = read.DNAmutation_count + 1
                        read.DNAmutation_location.append(i)
                        read.DNAmutation_identity.append(read.sequence[i])
            
                #translate, and enumerate protein matches/mismatches
                read.AAsequence = ''
                for i in xrange(0, len(referenceDNA), 3):
                    read.AAsequence += codon_table[read.sequence[i:i+3]]
            
                for i in range(0,len(referenceAA)):
                    if read.AAsequence[i] == referenceAA[i]:
                        read.AAmatch_count = read.AAmatch_count + 1
                    else: 
                        read.AAmutation_count = read.AAmutation_count + 1
                        read.AAmutation_location.append(i)
                        read.AAmutation_identity.append(read.AAsequence[i]) 
    
                if read.DNAmutation_count == 0:
                    read.DNAmutation_location = 'NA'
                    read.DNAmutation_identity = 'NA'
                    print >> f_DNA_output, '\t'.join(map(str, [read.ID, read.sequence, read.DNAmatch_count, read.DNAmutation_count, read.DNAmutation_location, read.DNAmutation_identity, read.maxmutrun]))
                
                else:
                    print >> f_DNA_output, '\t'.join(map(str, [read.ID, read.sequence, read.DNAmatch_count, read.DNAmutation_count, ','.join(map(str, read.DNAmutation_location)), ','.join(map(str, read.DNAmutation_identity)), read.maxmutrun]))
                    
                if read.AAmutation_count == 0:
                    read.AAmutation_location = 'NA'
                    read.AAmutation_identity = 'NA'
                    print >> f_protein_output, '\t'.join(map(str, [read.ID, read.AAsequence, read.AAmatch_count, read.AAmutation_count, read.AAmutation_location, read.AAmutation_identity, read.maxmutrun]))
                
                else:
                    print >> f_protein_output, '\t'.join(map(str, [read.ID, read.AAsequence, read.AAmatch_count, read.AAmutation_count, ','.join(map(str, read.AAmutation_location)), ','.join(map(str, read.AAmutation_identity)), read.maxmutrun]))

    f_infile.close()
    f_protein_output.close()
    f_DNA_output.close()
    return(0)
    
class FuserData:
    '''FuserData: A class to hold data from the read_fuser.py script'''
    def __init__(self, ID, sequence, total_count, mismatch_count, match_count, gap_count, unresolvable_count, maxmutrun, read1_avgquality, read1_chastity, read1_Ncount, read2_avgquality, read2_chastity, read2_Ncount, DNAmatch_count = 0, DNAmutation_count = 0, DNAmutation_location = [], DNAmutation_identity = [], AAsequence = 0, AAmatch_count = 0, AAmutation_count = 0, AAmutation_location = [], AAmutation_identity = []):
        self.ID = str(ID) #holds the entire Solexa ID tag
        self.sequence = str(sequence) #this holds the input sequences, the reverse read is reverse complemented and both are trimmed
        self.gap_count = int(gap_count) #this holds the number of gaps in the alignment
        self.mismatch_count = int(mismatch_count) #holds count of mismatches
        self.match_count = int(match_count) #holds count of matches
        self.total_count = int(total_count) #holds count of total bases 
        self.unresolvable_count = int(unresolvable_count) #holds count of unresolvable base pairs
        self.maxmutrun = int(maxmutrun) #holds length of maximum run of mutations 
        self.read1_avgquality = float(read1_avgquality) #holds average quality score from the read
        self.read1_chastity = str(read1_chastity) #did the read pass the chastity filter?
        self.read1_Ncount = int(read1_Ncount) #number of "N" bases in the read
        self.read2_avgquality = float(read2_avgquality) #holds average quality score from the read
        self.read2_chastity = str(read2_chastity) #did the read pass the chastity filter?
        self.read2_Ncount = int(read2_Ncount) #number of "N" bases in the read
        self.DNAmatch_count = int(DNAmatch_count) #holds match count for comparison to WT
        self.DNAmutation_count = int(DNAmutation_count) #holds mutation count for comparison to WT
        self.DNAmutation_location = list(DNAmutation_location) #holds mutation locations for comparison to WT
        self.DNAmutation_identity = list(DNAmutation_identity) #holds the mutant bases
        self.AAsequence = str(AAsequence) #holds AA sequence, translated from sequence
        self.AAmatch_count = int(AAmatch_count) #holds match count for comparison to WT
        self.AAmutation_count = int(AAmutation_count) #holds mutation count for comparison to WT
        self.AAmutation_location = list(AAmutation_location) #holds mutation locations for comparison to WT
        self.AAmutation_identity = list(AAmutation_identity) #holds the mutant bases
            
if __name__ == '__main__':
    print time.asctime(time.localtime())
    
    parser = optparse.OptionParser()
    parser.add_option('--path', action = 'store', type = 'string', dest = 'path', help = 'path from script to files')
    parser.add_option('--infile', action = 'store', type = 'string', dest = 'infile', help = 'input file, is the output of Paired_read_fuser.py')
    parser.add_option('--referenceDNA', action = 'store', type = 'string', dest = 'referenceDNA', help = 'reference DNA sequence')
    parser.add_option('--referenceAA', action = 'store', type = 'string', dest = 'referenceAA', help = 'reference amino acid sequence')
    parser.add_option('--gap_max', action = 'store', type = 'int', dest = 'gap_max', help = 'maximum number of paired read alignment gaps')
    parser.add_option('--unresolvable_max', action = 'store', type = 'int', dest = 'unresolvable_max', help = 'maximum number of unresolved bases in the pair')
    parser.add_option('--maxmutrun', action = 'store', type = 'int', dest = 'maxmutrun', help = 'maximum number of consecutive mutated bases, relative to WT')
    parser.add_option('--avg_quality', action = 'store', type = 'float', dest = 'avg_quality', help = 'minimum average quality score for each read')
    parser.add_option('--chaste', action = 'store', type = 'int', dest = 'chaste', help = '1 = use only reads that passed the chastity filter, 0 = use all reads')
    parser.add_option('--Ncount_max', action = 'store', type = 'string', dest = 'Ncount_max', help = 'maximum number of "N" bases in the read')
    parser.add_option('--mode', action = 'store', type = 'string', dest = 'mode', help = 'B = both reads, R1 = forward read only, R2 = reverse read only')
    (option, args) = parser.parse_args()

    main(option.path, option.infile, option.referenceDNA, option.referenceAA, option.gap_max, option.unresolvable_max, option.maxmutrun, option.avg_quality, option.chaste, option.Ncount_max, option.mode, option.local)
