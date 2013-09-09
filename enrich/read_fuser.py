#!/usr/bin/env python
'''
read_fuser: The read_fuser module performs several important tasks centered on forming a finished sequence for comparison to the wtDNA configuration element. 

If overlapping read pairs have been acquired, read_fuser stitches these together using the read1_overlap_start, read1_overlap_end, read2_overlap_start, read2_overlap_end, and include_nonoverlap_region parameters. Counting is done from 0, rather than 1. Read2 is considered to be reverse-complemented for the purposes of calculating the overlap start and overlap end.  Thus, in the example above, R2_overlap_start would be 0, whereas read1 start would be at position 22. If the include_nonoverlap_region configuration element were set to y, then the nonoverlap regions above would be appended to the overlapped sequence. Before a full-scale Enrich run is initiated, it is advisable to check that the overlap parameters are set correctly on a small subset of data.
    
For each forward and reverse read pair, read_fuser checks within the overlap region to make sure that the read pairs agree on the identity of each base. If they do not, then read_fuser examines the quality score of each base and picks the base with the higher score.  If both reads have an identical quality score, the base is considered unresolvable and an 'X' is inserted in the finished sequence at the unresolvable position.  Read_fuser also calculates a compound quality score estimate at each position, equal to the product of the two quality scores at that position if the two reads agree on the identity of the position.  If the reads do not agree, the higher quality score is used.  Positions in nonoverlapped regions are assigned their single-read quality score.

In the event that only one read is present, the overlap start and end configuration elements can be used to specify the beginning and end of the region of interest. Finsihed sequences within that region are generated. If only a reverse read is specified, then each read will be reverse-complemented before production of a finished sequence.

If the run_aligner configuration element is set to y, then a Needleman-Wunsch alignment is performed for each read pair which disagrees at more positions than specified in the paired_mismatch_threshold configuration element. Generally, read pairs that contain a large number of mismatches represent deletions within the region of interest and the Needleman-Wunch alignment produces a gap.
'''

import optparse, time, sys, os, array #import standard modules
    
__author__ = "Douglas M. Fowler"
__copyright__ = "Copyright 2011"
__credits__ = ["Douglas M Fowler", "Carlos L. Araya"]
__license__ = "FreeBSD"
__version__ = "0.2"
__maintainer__ = "Douglas M. Fowler"
__email__ = "dfowler@uw.edu"

def main(path, read1, read2, paired_mismatch_threshold, read1_overlap_start, read2_overlap_start, read1_overlap_end, read2_overlap_end, include_nonoverlap_region, wtseq, mode, run_aligner, chaste, grid = 'L'):

    if run_aligner == 'y':
        try:
            from Bio import pairwise2
        
        except:
            print 'Error: Could not import BioPython, which is needed for Needleman-Wunsch alignment.  Either install this dependency or set run_aligner = n'
            return(1)
        
    if grid != 'L':
        print time.asctime(time.localtime())
        print path, read1, read2, paired_mismatch_threshold, read1_overlap_start, read2_overlap_start, read1_overlap_end, read2_overlap_end, include_nonoverlap_region, wtseq, mode, run_aligner, chaste, grid

    try:
        paired_mismatch_threshold, read1_overlap_start, read2_overlap_start, read1_overlap_end, read2_overlap_end = map(int, [paired_mismatch_threshold, read1_overlap_start, read2_overlap_start, read1_overlap_end, read2_overlap_end])
    
    except:
        print 'Error: integer-only parameters were not integers'
        return 1
        
    # specify output file handle:
    if mode == 'B' or mode == 'R1':
        handle = read1.rstrip('.fq')
    
    elif mode == 'R2':
        handle = read2.rstrip('.fq')
        
    else:
        print 'Error: read_fuser was passed an invalid mode'
        return 1

    if grid == 'L':
        #open the files
        try:
            f_read1 = open((path + read1), 'U')
        
        except:
            if mode == 'B' or mode == 'R1':
                print 'Error: Read 1 input file not found'
                return 1
        
        try:
            f_read2 = open((path + read2), 'U') 
        
        except:
            if mode == 'B' or mode == 'R2':
                print 'Error: Read 2 input file not found'
                return 1    
        
        f_qc1 = open((path + handle + '_' + mode), 'w')
        print >> f_qc1, '\t'.join(['read1.ID', 'fused_sequence', 'total_bases', 'paired_mismatch_count', 'paired_match_count', 'gap_count', 'paired_unresolvable_count', 'maxmutrun', 'read1_avgquality', 'read1_chastity', 'read1_Ncount', 'read2_avgquality', 'read2_chastity', 'read2_Ncount'])
        
        if run_aligner == 'y' and mode == 'B':
            f_gap = open((path + handle + '_' + mode + '_gap'), 'w')
            print >> f_gap, '\t'.join(['read1.ID', 'fused_sequence', 'total_bases', 'paired_mismatch_count', 'paired_match_count', 'gap_count', 'paired_unresolvable_count', 'maxmutrun', 'read1_avgquality', 'read1_chastity', 'read1_Ncount', 'read2_avgquality', 'read2_chastity', 'read2_Ncount', 'read1_align', 'read2_align'])
    
    elif grid != 'L':
        #open the files
        try:
            f_read1 = open((path + read1 + '.' + os.getenv('SGE_TASK_ID')), 'U')
        
        except:
            if mode == 'B' or mode == 'R1':
                print 'Error: Read 1 not found'
                print path, read1
                return 1
        
        try:
            f_read2 = open((path + read2 + '.' + os.getenv('SGE_TASK_ID')), 'U')    
        
        except:
            if mode == 'B' or mode == 'R2':
                print 'Error: Read 2 not found'
                return 1
        
        #All data is written into f_qc1.  In addition, misaligned sequences are also written to f_gap for later inspection of the run_aligners
        f_qc1 = open((path + handle + '_' + mode + '.' + os.getenv('SGE_TASK_ID')), 'w')
        print >> f_qc1, '\t'.join(['read1.ID', 'fused_sequence', 'total_bases', 'paired_mismatch_count', 'paired_match_count', 'gap_count', 'paired_unresolvable_count', 'maxmutrun', 'read1_avgquality', 'read1_chastity', 'read1_Ncount', 'read2_avgquality', 'read2_chastity', 'read2_Ncount'])
        
        if run_aligner == 'y' and mode == 'B':
            f_gap = open((path + handle + '_' + mode + '_gap.' + os.getenv('SGE_TASK_ID')), 'w')
            print >> f_gap, '\t'.join(['read1.ID', 'fused_sequence', 'total_bases', 'paired_mismatch_count', 'paired_match_count', 'gap_count', 'paired_unresolvable_count', 'maxmutrun', 'read1_avgquality', 'read1_chastity', 'read1_Ncount', 'read2_avgquality', 'read2_chastity', 'read2_Ncount', 'read1_align', 'read2_align'])

    else:
        print 'Error: choose a valid local mode (e.g. L or SGE)'
        return 1
    
    #declare some variables important for both modes
    complement_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'} #dictionary for complementing a nucleotide sequence
    firstrun = 1
    len_wtseq = len(wtseq)
    score_adjustment = 33 #this is the numeric offset for the quality scores   
    lenread = 100000 #lenread specifies the number of bytes read in from each input file 
    lenbit = 1 #lenbit is a flag used to indicate when the last chunk of input from the file has been read
    remainder_read1 = ''
    remainder_read2 = ''

    if mode == 'B':     
        
        while lenbit: 
            #read in chunks of each input file and prepend the remainder, which is the incomplete set of lines left over from the previous loop
            read1_bytes = remainder_read1 + f_read1.read(lenread)
            read2_bytes = remainder_read2 + f_read2.read(lenread)

            #check to see if the end of file has arrived (i.e. if the length of the read from the input file is less than the specified length
            if len(read1_bytes) < lenread:
                lenbit = 0
                
                if len(read1_bytes) != len(read2_bytes):
                    print 'Error: length of input files not equal'
                    return(1)
                    
            #split the read into lines 
            read1_lines = read1_bytes.split('\n')
            read2_lines = read2_bytes.split('\n')
            
            #iterate over the split reads in four line chunks (since each illumina read consists of four lines)
            for i in range(0, ((len(read1_lines)/4)-lenbit)):
                #declare variables and whatnot
                read1 = SolexaData()
                read2 = SolexaData()
                paired_mismatch_count = 0
                paired_match_count = 0
                paired_unresolvable_count = 0
                gap_count = 0
                fused_sequence = ''
                wt_mismatch_count =0 
                wt_mismatch_location = []
                max_wt_mismatch = 0

                #take the first four lines from each read and then remove them from the list
                read1_entry = read1_lines[0:4]
                read1_lines[0:4] = []
                
                read2_entry = read2_lines[0:4]
                read2_lines[0:4] = []
                                
                read1.ID = read1_entry[0]
                read2.ID = read1_entry[0]

                if chaste == 'y': #check to see if the chastity bit is going to be used for filtration
                    read1.chastity = read1_entry[0][read1_entry[0].index('#')-1]
                    read2.chastity = read2_entry[0][read2_entry[0].index('#')-1]
                                        
                else: #if it's not, just set chastity to 1 (e.g. all reads are chaste)
                    read1.chastity = 1
                    read2.chastity = 1

                read1.sequence = read1_entry[1]
                temp_seq2 = read2_entry[1][::-1]
                read2.sequence = ''.join([complement_dict[base] for base in temp_seq2])
                
                if firstrun == 1:
                    len_read1 = len(read1.sequence)
                    len_read2 = len(read2.sequence)
                    firstrun = 0
                    
                read1.quality = [score - score_adjustment for score in array.array('b', read1_entry[3]).tolist()]
                #neato "extended slice" method for reversing a string... syntax is [start:end:step] so [::-1] is simple reversal
                read2.quality = [score - score_adjustment for score in array.array('b', read2_entry[3][::-1]).tolist()]
                
                read1.avgquality = sum(read1.quality)/float(len_read1)
                read2.avgquality = sum(read2.quality)/float(len_read2)
            
                if read1.ID[0:23] != read2.ID[0:23]:
                    print 'ERROR: read IDs not equal'
                    return 1    
                
                #this case is when there is at least some overlap between read pairs that needs resolving 
                if read1_overlap_end-read1_overlap_start > 0:
                    #loop over the sequences, resolving read pair mismatches in favor of the base with the highest quality score and inserting an X if the difference is unresolvable        

                    for i in xrange(0, (read1_overlap_end - read1_overlap_start + 1)):
                        r1index = i + read1_overlap_start
                        r2index = i + read2_overlap_start
                        
                        if read1.sequence[r1index] == read2.sequence[r2index]:
                            fused_sequence = fused_sequence + read1.sequence[r1index]
                            paired_match_count += 1
                        
                        elif read1.sequence[r1index] != read2.sequence[r2index]:
                            
                            if read1.quality[r1index] > read2.quality[r2index]:
                                fused_sequence = fused_sequence + read1.sequence[r1index]
                            
                            elif read1.quality[r1index] < read2.quality[r2index]:
                                fused_sequence = fused_sequence + read2.sequence[r2index]
                            
                            elif read1.quality[r1index] == read2.quality[r2index]:
                                fused_sequence = fused_sequence + 'X'
                                paired_unresolvable_count += 1
                            paired_mismatch_count += 1
                            
                        elif read1.sequence[r1index] == read2.sequence[r2index]:
                            fused_sequence = fused_sequence + read1.sequence[r1index]
                            paired_match_count += 1
                        
                    #this case deals with read pairs that are probably gapped, either because of a real indel causing extra sequence on each end or some other reason, by performing a Needleman-Wunsch alignment
                    if paired_mismatch_count >= paired_mismatch_threshold:
                        
                        if run_aligner == 'y':
                            
                            #reset all counters and fused_sequence since we will regenerate them from the alignment
                            paired_mismatch_count = 0
                            paired_match_count = 0
                            paired_unresolvable_count = 0
                            gap_count = 0
                            fused_sequence = ''
                            
                            #set the MAX_ALIGNMENTS variable in the local pairwise2 environment, otherwise it will return lots of sequences
                            pairwise2.MAX_ALIGNMENTS=1
                    
                            #do the alignment, match +2, mismatch -2, gap start = -3, gap extension = -1
                            alignment = pairwise2.align.globalms(read1.sequence[read1_overlap_start:read1_overlap_end+1], read2.sequence[read2_overlap_start:read2_overlap_end+1], 2, -1, -3, -1)
                            
                            #pull the aligned sequences out of the alignment                    
                            read1.align = str(alignment[0][0])
                            read2.align = str(alignment[0][1])
        
                            #fix up the aligned sequence and generate gap_count 
                            for j in xrange(0, read1_overlap_end - read1_overlap_start):
                                i = j + read1_overlap_start
                                
                                if read1.align[j] == '-':
                                    fused_sequence = fused_sequence + read2.align[j]
                                    gap_count += 1
                                
                                elif read2.align[j] == '-':
                                    fused_sequence = fused_sequence + read1.align[j]
                                    gap_count += 1
                                
                                elif read1.align[j] != read2.align[j]:
                                    
                                    if read1.quality[i] > read2.quality[i+read2_overlap_start]:
                                        fused_sequence = fused_sequence + read1.align[j]
                                    
                                    elif read1.quality[i] < read2.quality[i+read2_overlap_start]:
                                        fused_sequence = fused_sequence + read2.align[j]
                                    
                                    elif read1.quality[i] == read2.quality[i+read2_overlap_start]:
                                        fused_sequence = fused_sequence + 'X'
                                        paired_unresolvable_count += 1
                                    paired_mismatch_count += 1
                                
                                elif read1.align[j] == read2.align[j]:
                                    fused_sequence = fused_sequence + read1.align[j]
                                    paired_match_count += 1
                                                    
                    #append the 5' overhangs if that option was specified
                    if include_nonoverlap_region == 'y':
                            fused_sequence = read1.sequence[0:read1_overlap_start] + fused_sequence + read2.sequence[(read2_overlap_end + 1)::]
                
                else:
                    fused_sequence = read1.sequence + read2.sequence
        
                #align to WT here if the fused sequence is of the appropriate length
                len_fused_sequence = len(fused_sequence) 
                if len_fused_sequence == len_wtseq:
    
                    for i in xrange(0, len_wtseq):
                        
                        if fused_sequence[i] != wtseq[i]: 
                            wt_mismatch_count += 1
                            wt_mismatch_location.append(i)
    
                #calculate maxmutrun
                len_wt_mismatch_location = len(wt_mismatch_location)
                if len_wt_mismatch_location == 1:
                    max_wt_mismatch = 1
                
                else: 
                    
                    for j in xrange(0, len_wt_mismatch_location - 1): #j iterates over the list index
                        
                        for i in xrange(0, len_wt_mismatch_location - j): #i iterates from the jth list element to the end
                            
                            if wt_mismatch_location[j] == (wt_mismatch_location[j+i] - i): 
                                
                                if max_wt_mismatch < i + 1:
                                    max_wt_mismatch = i + 1
                                    
                #print to files here.  All sequences go into f_qc1.  Additionally, sequences with gaps go into f_gap with the alignment appended as extra two columns
                if run_aligner == 'y':
                
                    if gap_count == 0 and len_fused_sequence == len_wtseq:
                        print >> f_qc1, '\t'.join(map(str, [read1.ID, fused_sequence, len(fused_sequence), paired_mismatch_count, paired_match_count, gap_count, paired_unresolvable_count, max_wt_mismatch, read1.avgquality, read1.chastity, read1.sequence.count("N"), read2.avgquality, read2.chastity, read2.sequence.count("N")])) 
                    
                    else:    
                        print >> f_gap, '\t'.join(map(str, [read1.ID, fused_sequence, len(fused_sequence), paired_mismatch_count, paired_match_count, gap_count, paired_unresolvable_count, max_wt_mismatch, read1.avgquality, read1.chastity, read1.sequence.count("N"), read2.avgquality, read2.chastity, read2.sequence.count("N"), read1.align, read2.align]))             
                
                if run_aligner == 'n':
                     print >> f_qc1, '\t'.join(map(str, [read1.ID, fused_sequence, len_fused_sequence, paired_mismatch_count, paired_match_count, gap_count, paired_unresolvable_count, max_wt_mismatch, read1.avgquality, read1.chastity, read1.sequence.count("N"), read2.avgquality, read2.chastity, read2.sequence.count("N")]))   
                    
            remainder_read1 = '\n'.join(read1_lines) #since the input files are being read in chunks of 100000 bytes, the end of each chunk contains an incomplete illumina read.  these remainder line(s) are prepended to the next chunk at the beginning of the while loop
            remainder_read2 = '\n'.join(read2_lines)
           
        #close file handles and end
        f_read1.close()
        f_read2.close() 
        f_qc1.close()
        
        if run_aligner == 'y':
            f_gap.close()   
                        
        if __name__ == '__main__':
            print "GRACEFUL EXIT: EOF"
            print time.asctime(time.localtime()) #time to run                                                   
        
        if grid != 'L':
            print 'Job completed successfully at ' + time.asctime(time.localtime())
        
        return 0

    #similar loop as above, but for R1 or R2 mode
    elif mode == 'R1' or mode =='R2':
        remainder_readX = ''
        mode = int(mode.lstrip('R'))
        
        print 'mode = RX, using only readX'.replace('X', str(mode))
        
        if mode == 1:
            f_readX = f_read1
            readX_start = read1_overlap_start
            readX_end = read1_overlap_end
     
        if mode == 2:
            f_readX = f_read2
            readX_start = read2_overlap_start
            readX_end = read2_overlap_end
    
        while lenbit: 
            #read in chunks of each input file and prepend the remainder, which is the incomplete set of lines left over from the previous loop
            readX_bytes = remainder_readX + f_readX.read(lenread)

            #check to see if the end of file has arrived (i.e. if the length of the read from the input file is less than the specified length
            if len(readX_bytes) < lenread:
                lenbit = 0
                    
            #split the read into lines 
            readX_lines = readX_bytes.split('\n')
            
            #iterate over the split reads in four line chunks (since each illumina read consists of four lines)
            for i in range(0, ((len(readX_lines)/4)-lenbit)):
                #declare variables and whatnot
                readX = SolexaData()
                paired_mismatch_count = 'NA'
                paired_match_count = 'NA'
                paired_unresolvable_count = 'NA'
                gap_count = 'NA'
                fused_sequence = ''
                wt_mismatch_count =0 
                wt_mismatch_location = []
                max_wt_mismatch = 0

                #take the first four lines from each read and then remove them from the list
                readX_entry = readX_lines[0:4]
                readX_lines[0:4] = []
                                                
                readX.ID = readX_entry[0]

                if chaste == 'y': #check to see if the chastity bit is going to be used for filtration
                    readX.chastity = readX_entry[0][readX_entry[0].index('#')-1]
                                        
                else: #if it's not, just set chastity to 1 (e.g. all reads are chaste)
                    readX.chastity = 1

                
                if mode == 1:
                    readX.sequence = readX_entry[1]
                    readX.quality = [score - score_adjustment for score in array.array('b', readX_entry[3]).tolist()]
                    
                if mode == 2:
                    temp_seq2 = readX_entry[1][::-1]      
                    readX.sequence = ''.join([complement_dict[base] for base in temp_seq2])
                    readX.quality = [score - score_adjustment for score in array.array('b', readX_entry[3][::-1]).tolist()]                    
                
                if firstrun == 1:
                    len_readX = len(readX.sequence)
                    firstrun = 0
                             
                readX.avgquality = sum(readX.quality)/float(len_readX)

                #loop over the sequences pulling out relevant residues          
                for i in xrange(readX_start, readX_end + 1):
                    fused_sequence = fused_sequence + readX.sequence[i]
                                                                    
                #align to WT here if the fused sequence is of the appropriate length
                len_fused_sequence = len(fused_sequence)
                if len_fused_sequence == len_wtseq:
                    
                    for i in xrange(0, len_wtseq):
                        if fused_sequence[i] != wtseq[i]: 
                            wt_mismatch_count += 1
                            wt_mismatch_location.append(i)
    
                #calculate maxmutrun
                len_wt_mismatch_location = len(wt_mismatch_location)
                if len_wt_mismatch_location == 1:
                    max_wt_mismatch = 1
                    
                else: 
                    
                    for j in xrange(0, len_wt_mismatch_location - 1): #j iterates over the list index
                        
                        for i in xrange(0, len_wt_mismatch_location - j): #i iterates from the jth list element to the end
                            
                            if wt_mismatch_location[j] == (wt_mismatch_location[j+i] - i): 
                                
                                if max_wt_mismatch < i + 1:
                                    max_wt_mismatch = i + 1
                                    
                #print to files here.  All sequences go into f_qc1.
                print >> f_qc1, '\t'.join(map(str, [readX.ID, fused_sequence, len_fused_sequence, paired_mismatch_count, paired_match_count, gap_count, paired_unresolvable_count, max_wt_mismatch, readX.avgquality, readX.chastity, readX.sequence.count("N"), 'NA', 'NA', 'NA']))
            
            remainder_readX = '\n'.join(readX_lines) #since the input files are being read in chunks of 100000 bytes, the end of each chunk contains an incomplete illumina read.  these remainder line(s) are prepended to the next chunk at the beginning of the while loop
           
        #close file handles and end
        f_readX.close()
        f_qc1.close() 
                        
        if __name__ == '__main__':
            print "GRACEFUL EXIT: EOF"
            print time.asctime(time.localtime()) #time to run                                                   
        
        if grid != 'L':
            print 'Job completed successfully at ' + time.asctime(time.localtime())
        
        return 0
        
    else:
        print 'No data processed: specify a valid mode using --mode'
        
class SolexaData:
    "A class to hold a Solexa sequencer output entry"
    def __init__(self, ID = 0, sequence = 0, quality = 0, align = 0, avgquality = 0, chastity = ''):
        self.ID = str(ID) #holds the entire Solexa ID tag
        self.sequence = str(sequence) #this holds the input sequences, the reverse read is reverse complemented and both are trimmed
        self.quality = str(quality) #this holds the quality scores, the reverse read quality is reversed and both are trimmed
        self.align = str(align) #this holds the aligned sequences, the reverse read is reverse complemented and both are trimmed
        self.avgquality = float()
        self.chastity = str() 

if __name__ == '__main__':
                    
    print time.asctime(time.localtime())
    
    parser = optparse.OptionParser()
    parser.add_option('--path', action = 'store', type = 'string', dest = 'path', help = 'path from script to files')
    parser.add_option('--read1', action = 'store', type = 'string', dest = 'read1', help = 'Solexa formatted first set of reads')
    parser.add_option('--read2', action = 'store', type = 'string', dest = 'read2', help = 'Solexa formatted second set of reads')
    parser.add_option('--paired_mismatch_threshold', action = 'store', type = 'int', dest = 'paired_mismatch_threshold', help = 'maximum number of read pair mismatches that can occur before a read gets flagged for Needleman-Wunsch alignment')
    parser.add_option('--read1_overlap_start', action = 'store', type = 'int', dest = 'read1_overlap_start', help = 'beginning of overlapped region from read1')
    parser.add_option('--read2_overlap_start', action = 'store', type = 'int', dest = 'read2_overlap_start', help = 'beginning of overlapped region from read2, in the reverse complement')
    parser.add_option('--read1_overlap_end', action = 'store', type = 'int', dest = 'read1_overlap_end', help = 'end of overlapped region from read1')
    parser.add_option('--read2_overlap_end', action = 'store', type = 'int', dest = 'read2_overlap_end', help = 'end of overlapped region from read2, in the reverse complement')
    parser.add_option('--include_nonoverlap_region', action = 'store', type = 'string', dest = 'include_nonoverlap_region', help = 'if y, 5prime nonoverlapping regions will be appended.')
    parser.add_option('--wtseq', action = 'store', type = 'string', dest = 'wtseq', help = 'WT DNA sequence of assembled sequence, ex: CAGTACGAAACCCTGCTGATCGAAACCGCTTCTTCTCTGGTTAAAAACGCT')
    parser.add_option('--mode', action = 'store', type = 'string', dest = 'mode', help = 'B = use both reads, R1 = use read1 only, R2 = use read2 only.  Note that R1 and R2 options will use the appropriate readN_overlap_start and readN_overlap_end to extract the sequence in question')
    parser.add_option('--run_aligner', action = 'store', type = 'string', dest = 'run_aligner', help = 'If TRUE, a Needleman-Wunsch alignment will be performed for all read pairs that exceed the mismatch threshold.  Read pairs that produce a gapped alignment will be removed from the filtered output file and a _gap file will be created containing them.')
    parser.add_option('--local', action = 'store', type = 'string', dest = 'local', help = 'Is this a local run (L) or should an SGE (SGE) job be scheduled?')
    (option, args) = parser.parse_args()
                
    main(option.path, option.read1, option.read2, option.paired_mismatch_threshold, option.read1_overlap_start, option.read2_overlap_start, option.read1_overlap_end, option.read2_overlap_end, option.include_nonoverlap_region, option.wtseq, option.mode, option.run_aligner, option.local)
 
