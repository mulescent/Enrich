#!/usr/bin/env python 
'''
fastq_filter: The fastq_filter module filters a pair of FASTQ files to identify and selected reads with a particular sequence (or a particular index read sequence). 

The fastq_filter module exists to deal with multiplexed sequencing data sets, in which multiple libraries have been mixed and sequenced together.  Fastq_filter can run in one of three modes, identified by the index_mode element of the enrich run configuration file.  Index mode indicates that a separate index read is present.  For each read, fastq_filter examines the index sequence.  If the index sequence matches the *index_sequence* configuration element with no more than the number of mismatches specified in the index_mismatch_threshold configuration element, the read is considered a match and written into the output file.  If not, the read is discarded.  N_include and N_exclude modes indicate that an index read is not present.  These modes both compare the first 20 bases of each read to the supplied wtDNA configuration element and consider the read a match if it has no more than the number of mismatches specified in the index_mismatch_threshold configuration element.  N_include mode counts 'N' bases as mismatches, whereas N_exclude mode does not.
'''

import optparse, time, sys

__author__ = "Douglas M. Fowler"
__copyright__ = "Copyright 2011"
__credits__ = ["Douglas M Fowler", "Carlos L. Araya"]
__license__ = "FreeBSD"
__version__ = "0.2"
__maintainer__ = "Douglas M. Fowler"
__email__ = "dfowler@uw.edu"

def main(path, infile, mismatch_threshold, wtseq, index_file, index_sequence, mode, grid = 'L'):
    try: #make sure that integer parameters are really integers
        mismatch_threshold = int(mismatch_threshold)
    
    except:
        print 'Error: mismatch_threshold not an integer'
        return 1
        
    if grid != 'L': #print log info if called during a cluster run
        print time.asctime(time.localtime())
        print path, infile, mismatch_threshold, wtseq, index_file, index_sequence, mode, grid
        
    #lenread specifies the number of bytes read in from each input file
    lenread = 100000
    
    #below is the case for when an index read is present
    if 'index' in mode:
        
        #open files, see if one or two input files are specified
        try:
            infile_F, infile_R = infile.split(',')
            rev_read = 1 #indicate that a reverse read is present
            
        except:
            infile_F = infile
            rev_read = 0
            
        handle_F = infile_F.rstrip('.fq')
        
        try:
            f_infile_F = open(path + infile_F, 'U')
            f_index = open(path + index_file, 'U')
        
        except:
            print 'Error: input files could not be opened'
            return 1
            
        f_output_F = open(path + '../tmp/' + handle_F + '_' + mode + '_filtered.fq', 'w')
        
        if rev_read == 1:
            try:
                f_infile_R = open(path + infile_R, 'U') 
            
            except:
                print 'Error: input files could not be opened'
                return 1
            
            handle_R = infile_R.rstrip('.fq')   
            f_output_R = open(path + '../tmp/' + handle_R + '_' + mode + '_filtered.fq', 'w')
            
        #compute the length of the index sequence
        lenindex = len(index_sequence)
        remainder_index = ''
        
        #a set is a hashable disordered set of objects that is supports much faster x in set searches
        IDs = set()
        
        #lenbit is a flag used to indicate when the last chunk of input from the file has been read
        lenbit = 1
        
        #this loop builds a dictionary of illumina sequence IDs that have the correct index read sequence
        while lenbit:
                
            #read in chunks of each input file and prepend the remainder, which is the incomplete set of lines left over from the previous loop
            index_read = remainder_index + f_index.read(lenread)

            #check to see if the end of file has arrived (i.e. if the length of the read from the input file is less than the specified length
            if len(index_read) < lenread:
                lenbit = 0
                
            #split the read into lines 
            lines_index = index_read.split('\n')
            
            #iterate over the split reads in four line chunks (since each illumina read consists of four lines)
            for i in range(0, ((len(lines_index)/4)-lenbit)):
                mismatch_count_index = 0
                
                #take the first four lines from each read and then remove them from the list
                read_index = lines_index[0:4]
                lines_index[0:4] = []
                
                #iterate over the length of the index read checking to see if each position matches the specified index (xrange is used here for speed)
                for i in xrange(0, lenindex):
                    
                    if read_index[1][i] != index_sequence[i]:
                        mismatch_count_index += 1
                    
                if mismatch_count_index <= mismatch_threshold:
                    IDs.add(read_index[0].rsplit('/')[0])
                    
            #since the input files are being read in chunks of 100000 bytes, the end of each chunk contains an incomplete illumina read.  these remainder line(s) are prepended to the next chunk at the beginning of the while loop
            remainder_index = '\n'.join(lines_index)
        
        #close the index file
        f_index.close()
        
        #reset lenbit
        lenbit = 1
        remainder_input_F = ''
        remainder_input_R = ''
        
        #this loop builds a dictionary of illumina sequence IDs that have the correct index read sequence
        while lenbit:
            
            #define loop-specific variables
            outputlines_F = ''
            outputlines_R = ''
            
            #read in chunks of each input file and prepend the remainder, which is the incomplete set of lines left over from the previous loop
            input_read_F = remainder_input_F + f_infile_F.read(lenread)
            
            if rev_read == 1:
                input_read_R = remainder_input_R + f_infile_R.read(lenread)

            #check to see if the end of file has arrived (i.e. if the length of the read from the input file is less than the specified length
            if len(input_read_F) < lenread:
                lenbit = 0
            
            #make sure that the reads from the file are of equal length 
            if rev_read == 1:
            
                if len(input_read_F) != len(input_read_R):
                    print 'Error, file read lengths unequal'
                    return 1
                        
            #split the read into lines 
            lines_input_F = input_read_F.split('\n')
            
            if rev_read == 1:
            
                lines_input_R = input_read_R.split('\n')
            
            #iterate over the split reads in four line chunks (since each illumina read consists of four lines).  there are two cases, depending on whether a reverse read is present or noth
            if rev_read == 1:
            
                for i in range(0, ((len(lines_input_F)/4)-lenbit)):
    
                    #take the first four lines from each read and then remove them from the list
                    read_input_F = lines_input_F[0:4]
                    lines_input_F[0:4] = []
                    read_input_R = lines_input_R[0:4]
                    lines_input_R[0:4] = []
                    
                    if read_input_F[0].rsplit('/')[0] != read_input_R[0].rsplit('/')[0]:
                        print 'Error, read IDs not equal'
                        return 1
                    
                    if read_input_F[0].rsplit('/')[0] in IDs:
                        outputlines_F = outputlines_F + '\n'.join(read_input_F) + '\n'
                        outputlines_R = outputlines_R + '\n'.join(read_input_R) + '\n'
            
            if rev_read == 0:
            
                for i in range(0, ((len(lines_input_F)/4)-lenbit)):
    
                    #take the first four lines from each read and then remove them from the list
                    read_input_F = lines_input_F[0:4]
                    lines_input_F[0:4] = []
                    
                    if read_input_F[0].rsplit('/')[0] in IDs:
                        outputlines_F = outputlines_F + '\n'.join(read_input_F) + '\n'
        
            #print a chunk of filtered reads to disk        
            if len(outputlines_F) > 0:
                print >> f_output_F, outputlines_F.rstrip() 
                
                if rev_read == 1:
                    print >> f_output_R, outputlines_R.rstrip()
    
            #since the input files are being read in chunks of 100000 bytes, the end of each chunk contains an incomplete illumina read.  these remainder line(s) are prepended to the next chunk at the beginning of the while loop
            remainder_input_F = '\n'.join(lines_input_F)
            
            if rev_read == 1:
                remainder_input_R = '\n'.join(lines_input_R)

        #close the input and output files
        f_infile_F.close()
        f_output_F.close()
        
        if rev_read == 1:
            f_infile_R.close()
            f_output_R.close()
        
    #this is the case for when an index read is not present                         
    elif 'index' not in mode:
        #open the appropriate files for input and output, define variables
        #open files, see if one or two input files are specified
        
        try:
            infile_F, infile_R = infile.split(',')
            wtseq_F = wtseq
            complement_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'} #dictionary for complementing a nucleotide sequence
            wtseq_R = ''.join([complement_dict[base] for base in wtseq])
            
            rev_read = 1 #indicate that a reverse read is present
            
        except:
            infile_F = infile
            wtseq_F = wtseq
            rev_read = 0

        try:
            handle_F = infile_F.rstrip('.fq')
            f_infile_F = open(path + infile_F, 'U')
            f_output_F = open(path + '../tmp/' + handle_F + '_' + mode + '_filtered.fq', 'w')

        except:
            print 'Error: input files could not be opened'
            return 1
        
        if rev_read == 1:
            handle_R = infile_R.rstrip('.fq')
            
            try:
                f_infile_R = open(path + infile_R, 'U')
                
            except:
                print 'Error: input files could not be opened'
                return 1
                
            f_output_R = open(path + '../tmp/' + handle_R + '_' + mode + '_filtered.fq', 'w')
            
        #only the first 20 bases will be examined
        lenwt = 20
        
        #reset lenbit
        lenbit = 1
        remainder_F = ''
        remainder_R = ''
        
        while lenbit:
            #define loop-specific variables
            outputlines_F = ''
            outputlines_R = ''
                        
            #read in chunks of each input file and prepend the remainder, which is the incomplete set of lines left over from the previous loop
            raw_F = remainder_F + f_infile_F.read(lenread)
            
            if rev_read == 1:
                raw_R = remainder_R + f_infile_R.read(lenread)
            
            #check to see if the end of file has arrived (i.e. if the length of the read from the input file is less than the specified length
            if len(raw_F) < lenread:
                lenbit = 0
            
            #make sure that the reads from the file are of equal length 
            if rev_read == 1:
            
                if len(raw_F) != len(raw_R):
                    print 'Error, file read lengths unequal'
                    return 1
            
            #split the read into lines 
            lines_F = raw_F.split('\n')
            
            if rev_read == 1:
                lines_R = raw_R.split('\n')
            
            #iterate over the split reads in four line chunks (since each illumina read consists of four lines)
            if rev_read == 1:
            
                for i in range(0, ((len(lines_F)/4)-lenbit)):
                    mismatch_count_F = 0
                    mismatch_count_R = 0
                    
                    #take the first four lines from each read and then remove them from the list
                    read_F = lines_F[0:4]
                    lines_F[0:4] = []
                    read_R = lines_R[0:4]
                    lines_R[0:4] = []
                    
                    #pull the actual DNA sequences out of the each read
                    seq_F = read_F[1]
                    seq_R = read_R[1]
                    
                    #the 'N_exclude' case does not compare the sequence to the wild type if the sequence is an N
                    if mode == 'N_exclude':
                        
                        #iterate over the length of the read checking to see if each position matches the wild type (xrange is used here for speed)
                        for i in xrange(0, lenwt):
                            
                            if seq_F[i] != 'N' and seq_F[i] != wtseq_F[i]:
                                mismatch_count_F += 1
                            
                            if seq_R[i] != 'N' and seq_R[i] != wtseq_R[i]:
                                mismatch_count_R += 1
                    
                    #the 'N_include' option treats Ns just like any other base          
                    elif mode == 'N_include':
                        
                        #iterate over the length of the read checking to see if each position matches the wild type (xrange is used here for speed)
                        for i in xrange(0, lenwt):
                            
                            if seq_F[i] != wtseq_F[i]:
                                mismatch_count_F += 1
                            
                            if seq_R[i] != wtseq_R[i]:
                                mismatch_count_R += 1
                                
                    else:
                        print 'Error: choose a valid mode'
                        return 1
                    
                    #check to see if the maximum number of mismatches to the wildtype has been exceeded for either the forward or reverse read.  if not, write these reads to the filtered file     
                    if mismatch_count_F <= mismatch_threshold or mismatch_count_R <= mismatch_threshold:
                        outputlines_F = outputlines_F + '\n'.join(read_F) + '\n'
                        outputlines_R = outputlines_R + '\n'.join(read_R) + '\n'
            
            if rev_read == 0:
            
                for i in range(0, ((len(lines_F)/4)-lenbit)):
                    mismatch_count_F = 0
                    
                    #take the first four lines from each read and then remove them from the list
                    read_F = lines_F[0:4]
                    lines_F[0:4] = []
                    
                    #pull the actual DNA sequences out of the each read
                    seq_F = read_F[1]
                    
                    #the 'N_exclude' case does not compare the sequence to the wild type if the sequence is an N
                    if mode == 'N_exclude':
                        
                        #iterate over the length of the read checking to see if each position matches the wild type (xrange is used here for speed)
                        for i in xrange(0, lenwt):
                            
                            if seq_F[i] != 'N' and seq_F[i] != wtseq_F[i]:
                                mismatch_count_F += 1
                            
                    #the 'N_include' option treats Ns just like any other base          
                    elif mode == 'N_include':
                        
                        #iterate over the length of the read checking to see if each position matches the wild type (xrange is used here for speed)
                        for i in xrange(0, lenwt):
                            
                            if seq_F[i] != wtseq_F[i]:
                                mismatch_count_F += 1
                                
                    else:
                        print 'Error: choose a valid mode'
                        return 1
                    
                    #check to see if the maximum number of mismatches to the wildtype has been exceeded for either the forward or reverse read.  if not, write these reads to the filtered file     
                    if mismatch_count_F <= mismatch_threshold:
                        outputlines_F = outputlines_F + '\n'.join(read_F) + '\n'

        
            #print a chunk of output
            if len(outputlines_F) > 0:
                print >> f_output_F, outputlines_F.rstrip()
                
                if rev_read == 1:
                    print >> f_output_R, outputlines_R.rstrip()
            
            #since the input files are being read in chunks of 100000 bytes, the end of each chunk contains an incomplete illumina read.  these remainder line(s) are prepended to the next chunk at the beginning of the while loop
            remainder_F = '\n'.join(lines_F)
            
            if rev_read == 1:
                remainder_R = '\n'.join(lines_R)
        
        #close files
        f_infile_F.close()
        f_output_F.close()
        
        if rev_read == 1:
            f_infile_R.close()
            f_output_R.close()

    return 0
    
if __name__ == '__main__':

    print time.asctime(time.localtime())
    
    #set up the option parser to take command line parameters
    parser = optparse.OptionParser()
    parser.add_option('--path', action = 'store', type = 'string', dest = 'path', help = 'path from script to files')
    parser.add_option('--infile', action = 'store', type = 'string', dest = 'infile', help = 'FASTQ file, could be two comma separated files (forward and reverse reads)')
    
    parser.add_option('--mismatch_threshold', action = 'store', type = 'int', dest = 'mismatch_threshold', help = 'maximum number of mismatches that can occur before a read is discarded')
    parser.add_option('--wtseq', action = 'store', type = 'string', dest = 'wtseq', help = 'WT DNA sequence, ex: CAGTACGAAACCCTGCTGATCGAAACCGCTTCTTCTCTGGTTAAAAACGCT, could be two comma separated sequences(forward and reverse reads) if using N_include or N_exclude mode')
    parser.add_option('--index_file', action = 'store', type = 'string', dest = 'index_file', help = 'Solexa formatted read file with indices or NONE, if filtering is to be accomplished using just the wt sequence')
    parser.add_option('--index_sequence', action = 'store', type = 'string', dest = 'index_sequence', help = 'index sequence')
    parser.add_option('--mode', action = 'store', type = 'string', dest = 'mode', help = 'mode:  currently valid options are N_include (includes Ns as mismatches), N_exclude (excludes Ns as mismatches) and index (to use an index read to look up relevant reads)')
    parser.add_option('--local', action = 'store', type = 'string', dest = 'local', help = 'Is this a local run or should an SGE job be scheduled?')
    (option, args) = parser.parse_args()
    
    main(option.path, option.infile, option.mismatch_threshold, option.wtseq, option.index_file, option.index_sequence, option.mode, option.local)
     
    print time.asctime(time.localtime())
     
