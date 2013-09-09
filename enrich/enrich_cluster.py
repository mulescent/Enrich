#!/usr/bin/env python
'''
enrich_cluster: this module provides several functions, based on the DRMAA distributed resource management API, necessary for enrich to generate and execute grid engine jobs
'''

import drmaa, optparse, os, time, sys

__author__ = "Douglas M. Fowler"
__copyright__ = "Copyright 2011"
__credits__ = ["Douglas M Fowler", "Carlos L. Araya"]
__license__ = "FreeBSD"
__version__ = "0.2"
__maintainer__ = "Douglas M. Fowler"
__email__ = "dfowler@uw.edu"

def single_job(jt, job_path, arguments, output_path):
    '''
    single_job: generates a job template for a single DRMAA job
    '''
    
    jt.remoteCommand = job_path
    jt.outputPath = ':' + output_path
    jt.args = map(str, arguments) #the args attribute CANNOT accept lists with integers in them
    jt.jobEnvironment = os.environ #sets the job environment equal to the current execution environment
    jt.jobName = job_path.split('/')[-1]
    jt.blockEmail = True
    jt.joinFiles = True
    return(jt)      

def job_array(jt, job_path, arguments, output_path):
    '''
    job_array: generates a job template for a job array
    '''
    jt.jobEnvironment = os.environ
    jt.remoteCommand = job_path
    jt.outputPath = ':' + output_path   
    jt.args = [str(argument) for argument in arguments] #the args attribute CANNOT accept lists with integers in them
    jt.jobName = job_path.split('/')[-1]
    jt.blockEmail = True
    jt.joinFiles = True
    return(jt)
    
def fastq_file_break(project_directory, file_list):
    '''
    fastq_file_break: breaks input fastq files into chunks of 2000000 lines to parallelize the function of read_fuser
    '''

    if os.path.exists(project_directory + 'data/tmp/parts/') == 0:
        os.mkdir(project_directory + 'data/tmp/parts/')
    
    num_files = {}
    dirlist = os.listdir(project_directory + 'data/tmp/parts/')
    processed_file_list = set([item.split('.')[0] for item in dirlist])
    
    for item in file_list:

        if item.rstrip('.fq') not in processed_file_list and 'NONE' not in item and 'qc' not in item:
        
            print item + ' has not been atomized, this may take a few minutes'
            
            f = open(project_directory + 'data/tmp/' + item, 'U')
    
            x = 0
            l = f.readline().rstrip()
            c = []
            part_counter = 1
            o = open(project_directory + 'data/tmp/parts/' + item.rstrip('.fq') + '.' + str(part_counter), 'w')
            
            while l:
                if len(c) == 4:
                    print >> o, '\n'.join(c)
                    c = []
    
                elif x == 2000000:
                    part_counter = part_counter + 1
                    o.close()
                    o = open(project_directory + 'data/tmp/parts/' + item.rstrip('.fq') + '.' + str(part_counter), 'w')
                    x = 0
                    print 'part ' + str(part_counter) + ' completed'
                
                else:
                    c.append(l)
                    l = f.readline().rstrip()
                    x = x + 1
    
            if len(c) == 4:
                print >> o, '\n'.join(c)
            
            o.close()
            f.close()
            num_files[item] = part_counter
            
        elif item.rstrip('.fq') in processed_file_list:
            print item + ' appears to have been atomized already.  If the atomization process was interrupted, delete all parts files from project_directory/tmp/parts/'
            part_numbers = []
            
            for atomized in dirlist:
                
                if 'qc' not in atomized:            
                    parts = atomized.split('.')
                
                    if item.rstrip('.fq') in parts[0]:
                        
                        part_numbers.append(int(parts[1]))
                    
            num_files[item] = max(part_numbers)
        
    return(num_files)

def fastq_file_concatenate(project_directory, file_list):
    '''
    fastq_file_concatenate: takes the atomized read_fuser output files and generates a fused file
    ''' 

    dirlist = os.listdir(project_directory + 'data/tmp/parts/')

    for gap in [0,1]:
    
        for i in range(0, len(file_list)):
                                            
            if gap == 0 and 'NONE' not in file_list[i]:
                f_outfile = open((project_directory + 'data/tmp/' + file_list[i]), 'w') 
                header = 1
                
                for filename in dirlist:
                    if file_list[i] == filename.split('.')[0] and 'gap' not in filename:
                        
                        f_infile = open((project_directory + 'data/tmp/parts/' + filename), 'U')
                        l = f_infile.readline().rstrip()
                       
                        if header == 0:
                            l = f_infile.readline().rstrip()
                        
                        elif header == 1:
                            print >> f_outfile, l
                            header = 0
                            l = f_infile.readline().rstrip()
                        
                        while l:
                            print >> f_outfile, l
                            l = f_infile.readline().rstrip()
                        
                        f_infile.close()
                    
            if gap == 1 and 'NONE' not in file_list[i]:
                firstrun = 1
                                
                for filename in dirlist:

                    if file_list[i] in filename and 'gap' in filename:
                        
                        if firstrun == 1:
                            f_outfile = open((project_directory + 'data/tmp/' + file_list[i] + '_gap'), 'w')    
                            header = 1
                            firstrun = 0
                            
                        f_infile = open((project_directory + 'data/tmp/parts/' + filename), 'U')
                        l = f_infile.readline().rstrip()
                        
                        if header == 0:
                            l = f_infile.readline().rstrip()
                        
                        elif header == 1:
                            print >> f_outfile, l
                            header = 0
                            l = f_infile.readline().rstrip()
                        
                        while l:
                            print >> f_outfile, l
                            l = f_infile.readline().rstrip()
                        
                        f_infile.close()
    
def main():
    sys.exit()
    
if __name__=='__main__':
    parser = optparse.OptionParser()
    parser.add_option('--job_path', action = 'store', type = 'string', dest = 'job_path', help = 'job, including full path')
    parser.add_option('--args', action = 'store', type = 'string', dest = 'args', help = 'arguments')
    parser.add_option('--output_path', action = 'store', type = 'string', dest = 'output_path', help = 'output file, including full path')
    option, args = parser.parse_args()
    
    jobid = single_job(option.job_path, option.args, option.output_path)
    print jobid
    
    main()
