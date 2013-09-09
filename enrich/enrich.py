#!/usr/bin/env python
'''
enrich: the main implementation of Enrich

This script serves as the main user interface for the collection of enrich pipeline elements. It can be run in three difrerent modes: interactive, user command line, and raw command line. Interactive mode is launched from the command line and provides a series of text menus to guide the user through the process of initializing a project directory, creating a configuration file and running the pipeline. The user command line mode is for more advanced users and allows all elements of the pipeline to be run from the command line using --option style parameters (e.g. with a shell script for a large project). The raw command line mode is used by enrich itself when generating and running cluster jobs, and is not intended for direct user access.
'''

import sys, os, time, shutil, optparse, subprocess #imports of required standard modules

try:
    import drmaa, enrich_cluster #optional libraries, required for cluster
    
except:
    print 'Warning: DRMAA library not installed, cluster operations not available'

try:
    from scipy import stats

except:
    print 'Warning: SciPy.Stats library not installed, p-values will be reported as NaN'

try:
    import matplotlib

except:
    print 'Warning: Matplotlib not installed, plots will not be generated'


try:
    import qvalue

except:
    print 'Warning: Qvalue not installed, q-values will be reported as NaN'

try:
    import numpy

except:
    print 'Warning: Numpy not installed, plots will not be generated'

import enrich_config, enrich_xml_parser, fastq_filter, read_fuser, read_aligner, map_counts, map_ratios, map_unlink, enrich_plot #imports of enrich-specific libraries

__author__ = "Douglas M. Fowler"
__copyright__ = "Copyright 2011"
__credits__ = ["Douglas M Fowler", "Carlos L. Araya"]
__license__ = "FreeBSD"
__version__ = "0.2"
__maintainer__ = "Douglas M. Fowler"
__email__ = "dfowler@uw.edu"
     
def main():
    parser = optparse.OptionParser()
    parser.add_option('--mode', type = 'string', action = 'store', dest = 'mode', help = 'Currently valid modes are: initialize, configure, fastq_filter, read_fuse, read_align, map_counts, map_ratios, map_unlink, plot, run_all')
    parser.add_option('--config_file', type = 'string', action = 'store', dest = 'config_file', help = 'Path to the configuration file, which will specify all options for Enrich operations.  Note that configuration files can be created using the Enrich interactive mode or can be manually edited.  Choosing a configuration file is optional in interactive mode')
    option, args = parser.parse_args()
   
    valid_modes = ['initialize', 'configure', 'example', 'fastq_filter', 'read_fuse', 'read_align', 'map_counts', 'map_ratios', 'map_unlink', 'plot', 'run_all']
    if option.mode in valid_modes: # Look to see if a valid command line mode has been invoked and, if so, run it
        if option.mode == 'configure':
            enrich_config.main('main')
        
        elif option.mode == 'initialize':
            if len(args) == 1:
                project_directory_tool(args[0])
                
            else:
                sys.exit('Error: specify a path to a project directory to be validated or created')
                    
        elif option.mode == 'example':
            print 'Initializing example project directory...'
            if len(args) == 1:
                project_directory_tool(args[0], 'y')
                
            else:
                sys.exit('Error: specify a path to a project directory to be validated or created')   
       
        else: #all the following modes require a parsed configuration file, which is why they're grouped together
            try:
                cfg_data = enrich_xml_parser.main(option.config_file)
                cfg_data_flattened = enrich_config.enrich_config_flatten(cfg_data) #flatten configuration file to make passing arguments to modules easier
        
            except:
                sys.exit('Error: specified configuration file does not exist or is improperly formatted')
            
            if option.mode == 'fastq_filter':
                fq_filter(cfg_data_flattened, args)
                
            if option.mode == 'read_fuse':
                read_fuse(cfg_data_flattened)
            
            if option.mode == 'read_align':
                read_align(cfg_data_flattened, args)
            
            if option.mode == 'map_counts':
                counts(cfg_data_flattened, args)
                
            if option.mode == 'map_ratios':
                ratios(cfg_data_flattened, args)
                
            if option.mode == 'map_unlink':
                unlink(cfg_data_flattened, args)
                
            if option.mode == 'plot':
                plot(cfg_data_flattened)
                
            if option.mode == 'run_all':
                run_all(cfg_data_flattened)
    
    elif len(args) != 0: # Look to see if a bunch of arguments have been passed to the command line with no appropriate flags.  This option exists so that power users (or the grid engine) can call enrich easily 
        mode = args[0] # Get the intended mode
        args.pop(0) # Remove the mode from the arglist
        
        #try:
        if mode == 'fq_filter':
            fastq_filter.main(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7])
                    
        if mode == 'read_fuse':
            read_fuser.main(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10], args[11], args[12], args[13])
        
        if mode == 'read_align':
            read_aligner.main(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10], args[11])
        
        if mode == 'map_counts':
            map_counts.main(args[0], args[1], args[2])
            
        if mode == 'map_ratios':
            map_ratios.main(args[0], args[1], args[2], args[3])
            
        if mode == 'map_unlink':
            map_unlink.main(args[0], args[1], args[2], args[3], args[4], args[5])
                        
        #except:
#           sys.exit('Error: Command line mode failed')

    else: #if no valid command line options have been speficied, enter interactive mode
        cfg_order = ['general_info', 'fastq_filter', 'read_fuser', 'read_aligner', 'map_counts', 'map_ratios', 'map_unlink', 'enrich_plot'] #This list is for convenience - it enables the ordered printing of configuration information UPDATE
        
        print 'Enrich v0.2\n'
        while True: #loop to guide the choice or creation of a project directory
            project_directory = enrich_config.ainput('Choose a project directory, or type "example" to install an example project directory ', 1, '/path/to/directory/').getString()
            example = 'n'
            
            if project_directory == 'example':
                example = 'y'
                project_directory = enrich_config.ainput('Choose the example project directory ', 1, '/path/to/example/directory/').getString()
            
            if project_directory[-1] != '/': #append a / to the path if it does not exist
                project_directory = project_directory + '/'
            
            retcode = project_directory_tool(project_directory, example)
            
            if retcode == 0:
                break
                    
        while True: #main menu loop
            
            choice = enrich_config.ainput('Enrich menu\n\n1) Generate or modify a configuration file\n2) Load a configuration file\n3) View loaded configuration file\n4) Run a portion of the pipeline\n5) Run the entire pipeline\n6) Exit\n\nMake a choice: ', 1).getInteger()
            
            if choice == 1:
                cfg_file = enrich_config.main('notmain', project_directory) #run enrich_config to generate a configuration file
                print 'Created FX1\n'.replace('FX1', cfg_file)
                time.sleep(1)
                
            elif choice == 2:
                cfg_file = enrich_config.ainput('Please enter the name of the configuration file (located in project/input): ', 1).getString()
                
                try:
                    cfg_data = enrich_xml_parser.main(project_directory + 'input/' + cfg_file)
                    cfg_data_flattened = enrich_config.enrich_config_flatten(cfg_data) #flatten configuration file to make passing arguments to modules easier
                    print 'Configuration file FX1 loaded\n'.replace('FX1', cfg_file)
                
                except:
                    print 'Please choose a valid configuration file\n'
                
                time.sleep(1)
                
            elif choice == 3:
                try:
                    cfg_print(cfg_order, cfg_data) #print the loaded configuration file in the order specified by cfg_order    
                
                except:
                    print 'Please load a valid configuration file\n'
                
                time.sleep(1)
                    
            elif choice == 4:
                    
                while True: #run specific enrich element menu loop
                    try:
                        cfg_data_flattened.keys() #check to make sure that a configuration file has been loaded
                    
                    except:
                        print 'Load a valid configuration file first\n'
                        time.sleep(1)
                        break
                    
                    choice = enrich_config.ainput('Enrich pipeline menu\n\n1) Filter FASTQ files\n2) Fuse paired reads\n3) Quality filter, translate and align processed reads\n4) Quantitate unique protein sequences\n5) Calculate protein sequence ratios\n6) Unlink mutations\n7) Generate plots\n8) Exit to main menu\n\nMake a choice: ', 1).getInteger() 
                    
                    if choice == 1:
                        print 'Processing...'   
                        retcode = fq_filter(cfg_data_flattened)
                        
                        if retcode == 0:
                            print 'Input sequence filtering complete.  Filtered files can be found in FX1data/tmp/\n'.replace('FX1', cfg_data_flattened['path'])
                        
                        else:
                            print 'Error in fastq_filter, files not processed\n'
                        
                        time.sleep(1)
                                                
                    elif choice == 2:
                        print 'Processing...'
                        retcode = read_fuse(cfg_data_flattened)
                        
                        if retcode == 0:
                            print 'Paired read fusing complete.  Filtered files can be found in FX1data/tmp/\n'.replace('FX1', cfg_data_flattened['path'])
                        
                        else:
                            print 'Error in read_fuser, files not processed\n'
                        
                        time.sleep(1)
                        
                    elif choice == 3:
                        print 'Processing...'
                        retcode = read_align(cfg_data_flattened)
                        
                        if retcode == 0:
                            print 'Input translation and quality filtering complete.  Filtered files have a _DNA_qc (untranslated) or _PRO_qc (translated) suffix and can be found in FX1data/tmp/\n'.replace('FX1', cfg_data_flattened['path'])
                            
                        else: 
                            print 'Error in read_aligner, files not processed\n'
                            
                        time.sleep(1)
                        
                    elif choice == 4:
                        print 'Processing...'
                        retcode = counts(cfg_data_flattened)
                        
                        if retcode == 0:
                            print 'Counting of unique sequences complete.  Output files have a counts prefix and can be found in FX1data/output/\n'.replace('FX1', cfg_data_flattened['path'])
                            
                        else: 
                            print 'Error in map_counts, files not processed\n'
                            
                        time.sleep(1)
                    
                    elif choice == 5:
                        print 'Processing...'
                        retcode = ratios(cfg_data_flattened)
                        
                        if retcode == 0:
                            print 'Ratio calculations complete.  Output files have a ratio prefix and can be found in FX1data/output/\n'.replace('FX1', cfg_data_flattened['path'])
                            
                        else: 
                            print 'Error in map_ratios, files not processed\n'

                        
                    elif choice == 6:
                        print 'Processing...'
                        retcode = unlink(cfg_data_flattened)
                        
                        if retcode == 0:
                            print 'Unlink calculations complete.  Output files have an unlink prefix and can be found in FX1data/output/\n'.replace('FX1', cfg_data_flattened['path'])
                            
                        else: 
                            print 'Error in map_unlink, files not processed\n'
                    
                        time.sleep(1)
                        
                    elif choice == 7:
                        print 'Processing...'
                        retcode = plot(cfg_data_flattened)
                        
                        if retcode == 0:
                            print 'Plotting complete.  Plots can be found in FX1plots/\n'.replace('FX1', cfg_data_flattened['path'])
                            
                        else: 
                            print 'Error in enrich_plot, some files not processed\n'
                        
                        time.sleep(1)
                        
                    elif choice == 8:
                        break
                        
                    else:
                        'Make a valid numerical (e.g. 1) choice'
                            
            elif choice == 5:
                retcode = run_all(cfg_data_flattened)
                
            elif choice == 6:
                sys.exit()
                
            else:
                'Make a valid numerical (e.g. 1) choice'
                

def fq_filter(cfg_data_flattened, args=[]):
    '''
    fq_filter: provides an interface to the fastq_filter module
    '''
    
    if cfg_data_flattened['index_mode'] == 'NA': #if fastq_filter is not to be run, exit
        return 1
        
    retcodes = [] #a list to hold the return codes

    if cfg_data_flattened['local'] == 'L': #local mode case
        if cfg_data_flattened['input_read1_filename'] != 'NONE' or cfg_data_flattened['input_read2_filename'] != 'NONE':
            if cfg_data_flattened['input_read1_filename'] != 'NONE': #check to see if input read 1 exists
                infile = cfg_data_flattened['input_read1_filename'] #infile holds the set of reads to be filtered
                
                if cfg_data_flattened['input_read2_filename'] != 'NONE': #check to see if input_read 2 exists
                    infile = infile + ',' + cfg_data_flattened['input_read2_filename'] 
            
            elif cfg_data_flattened['input_read2_filename'] != 'NONE': #check to see if input_read 2 exists
                    infile = cfg_data_flattened['input_read2_filename'] 

            retcodes.append(fastq_filter.main(cfg_data_flattened['path'] + 'data/raw/', infile, int(cfg_data_flattened['index_mismatch_threshold']), cfg_data_flattened['wtDNA'], cfg_data_flattened['input_index_file'], cfg_data_flattened['index_sequence'], cfg_data_flattened['index_mode'], cfg_data_flattened['local'])) #call fastq_filter 
            
        
        if cfg_data_flattened['sel_read1_filename'] != 'NONE' or cfg_data_flattened['sel_read2_filename'] != 'NONE':
            if cfg_data_flattened['sel_read1_filename'] != 'NONE': #check to see if selected read 1 exists
                infile = cfg_data_flattened['sel_read1_filename']   #infile holds the set of reads to be filtered
                
                if cfg_data_flattened['sel_read2_filename'] != 'NONE': #check to see if selected read 2 exists
                    infile = infile + ',' + cfg_data_flattened['sel_read2_filename'] 
            
            elif cfg_data_flattened['sel_read2_filename'] != 'NONE': #check to see if input_read 2 exists
                infile = cfg_data_flattened['sel_read2_filename']  
           
            retcodes.append(fastq_filter.main(cfg_data_flattened['path'] + 'data/raw/', infile, int(cfg_data_flattened['index_mismatch_threshold']), cfg_data_flattened['wtDNA'], cfg_data_flattened['sel_index_file'], cfg_data_flattened['index_sequence'], cfg_data_flattened['index_mode'], cfg_data_flattened['local'])) #call fastq_filter 
    
    
    elif cfg_data_flattened['local']  == 'SGE': #SGE grid case
        try:
            jobids = [] #list of job ids
            s = drmaa.Session()
            s.initialize()
            
            if cfg_data_flattened['input_read1_filename'] != 'NONE' or cfg_data_flattened['input_read2_filename'] != 'NONE':
                if cfg_data_flattened['input_read1_filename'] != 'NONE': #check to see if input read 1 exists
                    infile = cfg_data_flattened['input_read1_filename'] #infile holds the set of reads to be filtered
                    
                    if cfg_data_flattened['input_read2_filename'] != 'NONE': #check to see if input_read 2 exists
                        infile = infile + ',' + cfg_data_flattened['input_read2_filename'] 
                
                elif cfg_data_flattened['input_read2_filename'] != 'NONE': #check to see if input_read 2 exists
                    infile = cfg_data_flattened['input_read2_filename'] 

                job = enrich_cluster.single_job(s.createJobTemplate(), 'enrich', ['fq_filter', cfg_data_flattened['path'] + 'data/raw/', infile, int(cfg_data_flattened['index_mismatch_threshold']), cfg_data_flattened['wtDNA'], cfg_data_flattened['input_index_file'], cfg_data_flattened['index_sequence'], cfg_data_flattened['index_mode'], cfg_data_flattened['local']], cfg_data_flattened['path'] + 'log/' + 'fastq_filter.log') #submit a fastq_filter job to the cluster 
                jobids.append(s.runJob(job))
                
            if cfg_data_flattened['sel_read1_filename'] != 'NONE' or cfg_data_flattened['sel_read2_filename'] != 'NONE':
                if cfg_data_flattened['sel_read1_filename'] != 'NONE': #check to see if selected read 1 exists
                    infile = cfg_data_flattened['sel_read1_filename']   #infile holds the set of reads to be filtered
                    
                    if cfg_data_flattened['sel_read2_filename'] != 'NONE': #check to see if selected read 2 exists
                        infile = infile + ',' + cfg_data_flattened['sel_read2_filename'] 
                
                elif cfg_data_flattened['sel_read2_filename'] != 'NONE': #check to see if input_read 2 exists
                    infile = cfg_data_flattened['sel_read2_filename'] 

                job = enrich_cluster.single_job(s.createJobTemplate(), 'enrich', ['fq_filter', cfg_data_flattened['path'] + 'data/raw/', infile, int(cfg_data_flattened['index_mismatch_threshold']), cfg_data_flattened['wtDNA'], cfg_data_flattened['sel_index_file'], cfg_data_flattened['index_sequence'], cfg_data_flattened['index_mode'], cfg_data_flattened['local']], cfg_data_flattened['path'] + 'log/' + 'fastq_filter.log') #submit a fastq_filter job to the cluster
                jobids.append(s.runJob(job))
            
            print 'Job(s) submitted to cluster'
            print ', '.join(jobids)
            
            if 'nowait' not in args: #if the nowait option has not been specified, wait until jobs have been completed
                print 'Waiting for job completion...'
                s.synchronize(jobids, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)
            
            s.deleteJobTemplate(job)    
            s.exit()
            return 0

        except:
            return 1
            
    if 1 in retcodes:
        return 1
        
    else:
        return 0
                

def read_fuse(cfg_data_flattened):
    '''
    read_fuse: provides an interface to the read_fuser module
    '''
    
    retcodes = [] #a list to hold the return codes
    
    if cfg_data_flattened['index_mode'] == 'NA': #if fastq_filter is not to be run, create a symbolic link back to the data in data/raw
        for item in [cfg_data_flattened['input_read1_filename'], cfg_data_flattened['input_read2_filename'], cfg_data_flattened['sel_read1_filename'], cfg_data_flattened['sel_read2_filename']]:
            command = 'ln -s PATHdata/raw/NAME PATHdata/tmp/NAME'.replace('PATH', cfg_data_flattened['path']).replace('NAME', item)
           
            try:
                goo = ''
                foo = ''
                retcode = subprocess.call(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            
            except:
                print 'File %s not linked' % item
                
    if cfg_data_flattened['local'] == 'L': #local mode case
    
        if cfg_data_flattened['input_read1_filename'] != 'NONE' or cfg_data_flattened['input_read2_filename'] != 'NONE': 
            retcodes.append(read_fuser.main(cfg_data_flattened['path'] + 'data/tmp/', cfg_data_flattened['input_read1_filename_fuser'], cfg_data_flattened['input_read2_filename_fuser'], int(cfg_data_flattened['paired_mismatch_threshold']), int(cfg_data_flattened['read1_overlap_start']), int(cfg_data_flattened['read2_overlap_start']), int(cfg_data_flattened['read1_overlap_end']), int(cfg_data_flattened['read2_overlap_end']), cfg_data_flattened['include_nonoverlap_region'], cfg_data_flattened['wtDNA'], cfg_data_flattened['fuser_mode'], cfg_data_flattened['run_aligner'], cfg_data_flattened['chaste'], cfg_data_flattened['local']))
            
        if cfg_data_flattened['sel_read1_filename'] != 'NONE' or cfg_data_flattened['sel_read2_filename'] != 'NONE': 
            retcodes.append(read_fuser.main(cfg_data_flattened['path'] + 'data/tmp/', cfg_data_flattened['sel_read1_filename_fuser'], cfg_data_flattened['sel_read2_filename_fuser'], int(cfg_data_flattened['paired_mismatch_threshold']), int(cfg_data_flattened['read1_overlap_start']), int(cfg_data_flattened['read2_overlap_start']), int(cfg_data_flattened['read1_overlap_end']), int(cfg_data_flattened['read2_overlap_end']), cfg_data_flattened['include_nonoverlap_region'], cfg_data_flattened['wtDNA'], cfg_data_flattened['fuser_mode'], cfg_data_flattened['run_aligner'], cfg_data_flattened['chaste'], cfg_data_flattened['local']))
    
    elif cfg_data_flattened['local'] == 'SGE':
        try:
            print 'Checking for atomized FASTQ files (found in project/data/tmp/)'
            num_files = enrich_cluster.fastq_file_break(cfg_data_flattened['path'], [cfg_data_flattened['input_read1_filename_fuser'], cfg_data_flattened['input_read2_filename_fuser'], cfg_data_flattened['sel_read1_filename_fuser'], cfg_data_flattened['sel_read2_filename_fuser']]) #fastq_file_break returns a dictionary of input filenames and the corresponding number of parts they have been broken into
            
            s = drmaa.Session()
            s.initialize() 
            jobids = []
            
            if cfg_data_flattened['input_read1_filename'] != 'NONE' or cfg_data_flattened['input_read2_filename'] != 'NONE':
                job = enrich_cluster.job_array(s.createJobTemplate(), 'enrich', ['read_fuse', cfg_data_flattened['path'] + 'data/tmp/parts/', cfg_data_flattened['input_read1_filename_fuser'].rstrip('.fq'), cfg_data_flattened['input_read2_filename_fuser'].rstrip('.fq'), int(cfg_data_flattened['paired_mismatch_threshold']), int(cfg_data_flattened['read1_overlap_start']), int(cfg_data_flattened['read2_overlap_start']), int(cfg_data_flattened['read1_overlap_end']), int(cfg_data_flattened['read2_overlap_end']), cfg_data_flattened['include_nonoverlap_region'], cfg_data_flattened['wtDNA'], cfg_data_flattened['fuser_mode'], cfg_data_flattened['run_aligner'], cfg_data_flattened['chaste'], cfg_data_flattened['local']], cfg_data_flattened['path'] + 'log/' + 'read_fuser.log') #submit an array of read_fuser jobs to the cluster
    
                jobids.extend(s.runBulkJobs(job, 1, num_files[cfg_data_flattened['input_read1_filename_fuser']], 1))
                s.deleteJobTemplate(job)
                
            if cfg_data_flattened['sel_read1_filename'] != 'NONE' or cfg_data_flattened['sel_read2_filename'] != 'NONE': 
                job = enrich_cluster.job_array(s.createJobTemplate(), 'enrich', ['read_fuse', cfg_data_flattened['path'] + 'data/tmp/parts/', cfg_data_flattened['sel_read1_filename_fuser'].rstrip('.fq'), cfg_data_flattened['sel_read2_filename_fuser'].rstrip('.fq'), int(cfg_data_flattened['paired_mismatch_threshold']), int(cfg_data_flattened['read1_overlap_start']), int(cfg_data_flattened['read2_overlap_start']), int(cfg_data_flattened['read1_overlap_end']), int(cfg_data_flattened['read2_overlap_end']), cfg_data_flattened['include_nonoverlap_region'], cfg_data_flattened['wtDNA'], cfg_data_flattened['fuser_mode'], cfg_data_flattened['run_aligner'], cfg_data_flattened['chaste'], cfg_data_flattened['local']], cfg_data_flattened['path'] + 'log/' + 'read_fuser.log') #submit an array of read_fuser jobs to the cluster
                
                jobids.extend(s.runBulkJobs(job, 1, num_files[cfg_data_flattened['sel_read1_filename_fuser']], 1))
                s.deleteJobTemplate(job)
                
            print 'Job(s) submitted to cluster'
            print ', '.join(jobids)
            print 'Waiting for job completion...'
            s.synchronize(jobids, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)
            s.exit()
            print 'Jobs completed, concatenating output files'
    
            enrich_cluster.fastq_file_concatenate(cfg_data_flattened['path'], [cfg_data_flattened['input_filename_aligner'], cfg_data_flattened['sel_filename_aligner']])
            return 0
        
        except:
            return 1
            
    if 1 in retcodes:
        return 1
        
    else:
        return 0
        

def read_align(cfg_data_flattened, args=[]):
    '''
    read_align: provides an interface to the read_aligner module
    '''
   
    retcodes = [] #a list to hold the return codes
    
    if cfg_data_flattened['local'] == 'L':
        if cfg_data_flattened['input_read1_filename'] != 'NONE' or cfg_data_flattened['input_read2_filename'] != 'NONE':
            retcodes.append(read_aligner.main(cfg_data_flattened['path'], cfg_data_flattened['input_filename_aligner'], cfg_data_flattened['wtDNA'], cfg_data_flattened['wtPRO'], cfg_data_flattened['gap_max'], cfg_data_flattened['unresolvable_max'], cfg_data_flattened['maximum_mutation_run'], cfg_data_flattened['avg_quality'], cfg_data_flattened['chaste'], cfg_data_flattened['Ncount_max'], cfg_data_flattened['fuser_mode'], cfg_data_flattened['local'])) #call read_aligner 
            
        
        if cfg_data_flattened['sel_read1_filename'] != 'NONE' or cfg_data_flattened['sel_read2_filename'] != 'NONE':
            retcodes.append(read_aligner.main(cfg_data_flattened['path'], cfg_data_flattened['sel_filename_aligner'], cfg_data_flattened['wtDNA'], cfg_data_flattened['wtPRO'], cfg_data_flattened['gap_max'], cfg_data_flattened['unresolvable_max'], cfg_data_flattened['maximum_mutation_run'], cfg_data_flattened['avg_quality'], cfg_data_flattened['chaste'], cfg_data_flattened['Ncount_max'], cfg_data_flattened['fuser_mode'], cfg_data_flattened['local'])) #call read_aligner 
    
    elif cfg_data_flattened['local']  == 'SGE':
        try:
            jobids = [] #list of job ids
            s = drmaa.Session()
            s.initialize()
            
            if cfg_data_flattened['input_read1_filename'] != 'NONE' or cfg_data_flattened['input_read2_filename'] != 'NONE':
                job = enrich_cluster.single_job(s.createJobTemplate(), 'enrich', ['read_align', cfg_data_flattened['path'], cfg_data_flattened['input_filename_aligner'], cfg_data_flattened['wtDNA'], cfg_data_flattened['wtPRO'], cfg_data_flattened['gap_max'], cfg_data_flattened['unresolvable_max'], cfg_data_flattened['maximum_mutation_run'], cfg_data_flattened['avg_quality'], cfg_data_flattened['chaste'], cfg_data_flattened['Ncount_max'], cfg_data_flattened['fuser_mode'], cfg_data_flattened['local']], cfg_data_flattened['path'] + 'log/' + 'read_aligner.log') #submit a read_aligner job to the cluster 
                jobids.append(s.runJob(job))
                
            if cfg_data_flattened['sel_read1_filename'] != 'NONE' or cfg_data_flattened['sel_read2_filename'] != 'NONE':
                job = enrich_cluster.single_job(s.createJobTemplate(), 'enrich', ['read_align', cfg_data_flattened['path'], cfg_data_flattened['sel_filename_aligner'], cfg_data_flattened['wtDNA'], cfg_data_flattened['wtPRO'], cfg_data_flattened['gap_max'], cfg_data_flattened['unresolvable_max'], cfg_data_flattened['maximum_mutation_run'], cfg_data_flattened['avg_quality'], cfg_data_flattened['chaste'], cfg_data_flattened['Ncount_max'], cfg_data_flattened['fuser_mode'], cfg_data_flattened['local']], cfg_data_flattened['path'] + 'log/' + 'read_aligner.log') #submit a read_aligner job to the cluster
                jobids.append(s.runJob(job))
            
            print 'Job(s) submitted to cluster'
            print ', '.join(jobids)
           
            if 'nowait' not in args:
                print 'Waiting for job completion...'
                s.synchronize(jobids, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

            s.deleteJobTemplate(job)
            s.exit()
            return 0

        except:
            return 1
            
    if 1 in retcodes:
        return 1
        
    else:
        return 0

def counts(cfg_data_flattened, args=[]):
    '''    
    counts: provides an interface to the map_counts module
    '''
    
    retcodes = [] #a list to hold the return codes
    
    input_files = cfg_data_flattened['input_filename_map_counts'].split(',')
    sel_files = cfg_data_flattened['sel_filename_map_counts'].split(',')
    
    if cfg_data_flattened['local'] == 'L':
        if cfg_data_flattened['input_read1_filename'] != 'NONE' or cfg_data_flattened['input_read2_filename'] != 'NONE':
            retcodes.append(map_counts.main(cfg_data_flattened['path'], input_files[0], cfg_data_flattened['local'])) #call map_counts
            retcodes.append(map_counts.main(cfg_data_flattened['path'], input_files[1], cfg_data_flattened['local'])) #call map_counts
        
        if cfg_data_flattened['sel_read1_filename'] != 'NONE' or cfg_data_flattened['sel_read2_filename'] != 'NONE':
            retcodes.append(map_counts.main(cfg_data_flattened['path'], sel_files[0], cfg_data_flattened['local'])) #call map_counts
            retcodes.append(map_counts.main(cfg_data_flattened['path'], sel_files[1], cfg_data_flattened['local'])) #call map_counts

    elif cfg_data_flattened['local']  == 'SGE':
        try:
            jobids = [] #list of job ids
            s = drmaa.Session()
            s.initialize()
            
            if cfg_data_flattened['input_read1_filename'] != 'NONE' or cfg_data_flattened['input_read2_filename'] != 'NONE':
                job = enrich_cluster.single_job(s.createJobTemplate(), 'enrich', ['map_counts', cfg_data_flattened['path'], input_files[0], cfg_data_flattened['local']], cfg_data_flattened['path'] + 'log/' + 'map_counts.log') #submit a map_counts job to the cluster
                jobids.append(s.runJob(job))
                
                job = enrich_cluster.single_job(s.createJobTemplate(), 'enrich', ['map_counts', cfg_data_flattened['path'], input_files[1], cfg_data_flattened['local']], cfg_data_flattened['path'] + 'log/' + 'map_counts.log') #submit a map_counts job to the cluster
                jobids.append(s.runJob(job))
                
            if cfg_data_flattened['sel_read1_filename'] != 'NONE' or cfg_data_flattened['sel_read2_filename'] != 'NONE':
                job = enrich_cluster.single_job(s.createJobTemplate(), 'enrich', ['map_counts', cfg_data_flattened['path'], sel_files[0], cfg_data_flattened['local']], cfg_data_flattened['path'] + 'log/' + 'map_counts.log') #submit a map_counts job to the cluster
                jobids.append(s.runJob(job))
                
                job = enrich_cluster.single_job(s.createJobTemplate(), 'enrich', ['map_counts', cfg_data_flattened['path'], sel_files[1], cfg_data_flattened['local']], cfg_data_flattened['path'] + 'log/' + 'map_counts.log') #submit a map_counts job to the cluster
                jobids.append(s.runJob(job))

            print 'Job(s) submitted to cluster'
            print ', '.join(jobids)

            if 'nowait' not in args:
                print 'Waiting for job completion...'
                s.synchronize(jobids, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)
                
            s.deleteJobTemplate(job)
            s.exit()
            return 0

        except:
            return 1
            
    if 1 in retcodes:
        return 1
        
    else:
        return 0

def ratios(cfg_data_flattened, args=[]):
    '''
    ratios: provides an interface to the map_ratios module
    '''

    retcodes = [] #a list to hold the return codes
    
    input_files = cfg_data_flattened['input_filename_map_ratios'].split(',')
    sel_files = cfg_data_flattened['sel_filename_map_ratios'].split(',')
    
    if cfg_data_flattened['local'] == 'L':
        retcodes.append(map_ratios.main(cfg_data_flattened['path'], sel_files[0], input_files[0], cfg_data_flattened['local'])) #call map_ratios
        retcodes.append(map_ratios.main(cfg_data_flattened['path'], sel_files[1], input_files[1], cfg_data_flattened['local'])) #call map_ratios
            
    elif cfg_data_flattened['local']  == 'SGE':
        try:
            jobids = [] #list of job ids
            s = drmaa.Session()
            s.initialize()
            
            job = enrich_cluster.single_job(s.createJobTemplate(), 'enrich', ['map_ratios', cfg_data_flattened['path'], sel_files[0], input_files[0], cfg_data_flattened['local']], cfg_data_flattened['path'] + 'log/' + 'map_ratios.log') #submit a map_ratios job to the cluster 
            jobids.append(s.runJob(job))
                
            job = enrich_cluster.single_job(s.createJobTemplate(), 'enrich', ['map_ratios', cfg_data_flattened['path'], sel_files[1], input_files[1], cfg_data_flattened['local']], cfg_data_flattened['path'] + 'log/' + 'map_ratios.log') #submit a map_ratios job to the cluster 
            jobids.append(s.runJob(job))
                
            print 'Job(s) submitted to cluster'
            print ', '.join(jobids)
            
            if 'nowait' not in args:
                print 'Waiting for job completion...'
                s.synchronize(jobids, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)
            
            s.deleteJobTemplate(job)    
            s.exit()
            return 0

        except:
            return 1
            
    if 1 in retcodes:
        return 1
        
    else:
        return 0
        
def unlink(cfg_data_flattened, args=[]):
    '''
    unlink: provides an interface to the map_unlink module
    '''

    retcodes = [] #a list to hold the return codes
    molecule = '' #a string to hold the molecule type
    length = '' #length of the wtPRO sequence
    files = []
    files.extend(cfg_data_flattened['input_filename_map_unlink'].split(','))
    files.extend(cfg_data_flattened['sel_filename_map_unlink'].split(','))
    unlink_modes = cfg_data_flattened['unlink_modes'].split(',')
    
    if cfg_data_flattened['local'] == 'L':  
        for infile in files:
            if 'DNA' in infile:
                molecule = 'DNA'
                length = len(cfg_data_flattened['wtDNA'])
                
            if 'PRO' in infile:
                molecule = 'PRO'
                
                try:
                    length = len(cfg_data_flattened['wtDNA'])/3
                
                except:
                    print 'Error, length of DNA not divisible by 3 (i.e. not a valid coding sequence)'
                    retcodes.append(1)
            
            for unlink_mode in unlink_modes:
                    retcodes.append(map_unlink.main(cfg_data_flattened['path'], infile, molecule, unlink_mode, length, cfg_data_flattened['local'])) #call map_unlink
                    
    if cfg_data_flattened['local']  == 'SGE':   
        jobids = [] #list of job ids
        s = drmaa.Session()
        s.initialize()
        
        for infile in files:
            if 'DNA' in infile:
                molecule = 'DNA'
                length = len(cfg_data_flattened['wtDNA'])
                
            if 'PRO' in infile:
                molecule = 'PRO'
                
                try:
                    length = len(cfg_data_flattened['wtDNA'])/3
                
                except:
                    print 'Error, length of DNA not divisible by 3 (i.e. not a valid coding sequence)'
                    retcodes.append(1)
            
            for unlink_mode in unlink_modes:
                try:
                    job = enrich_cluster.single_job(s.createJobTemplate(), 'enrich', ['map_unlink', cfg_data_flattened['path'], infile, molecule, unlink_mode, length, cfg_data_flattened['local']], cfg_data_flattened['path'] + 'log/' + 'map_unlink.log') #submit a map_unlink job to the cluster
                    jobids.append(s.runJob(job))
                                
                except:
                    return 1    
                    
        print 'Job(s) submitted to cluster'
        print ', '.join(jobids)
        
        if 'nowait' not in args:
            print 'Waiting for job completion...'
            s.synchronize(jobids, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

        s.deleteJobTemplate(job)
        s.exit()
        return 0

    if 1 in retcodes:
        return 1
        
    else:
        return 0

def plot(cfg_data_flattened):
    '''
    plot: provides an interface to the enrich_plot module
    '''

    retcodes = [] #a list to hold the return codes
    files = cfg_data_flattened['input_filename_freqplot'].split(',')
    
    for infile in files:
        if 'NONE' not in infile:
            retcodes.append(enrich_plot.freqplot(cfg_data_flattened['path'], infile)) #call enrich_plot.freqplot
    
    ratioplot_infiles = cfg_data_flattened['input_filename_ratioplot'].split(';')
    modes = ['AA', 'Pos']
     
    for item in ratioplot_infiles:
        ratioplot_infile = item.split(',')
        if 'NONE' not in ratioplot_infile[0] and 'NONE' not in ratioplot_infile[1]:
            retcodes.append(enrich_plot.ratioplot(cfg_data_flattened['path'], ratioplot_infile[0], ratioplot_infile[1])) #call enrich_plot.ratioplot
        
        for mode in modes:
            retcodes.append(enrich_plot.all_residue_plot(cfg_data_flattened['path'], ratioplot_infile[0], ratioplot_infile[1], cfg_data_flattened['wtPRO'], mode)) #call enrich_plot.all_residue_plot
    
    if 1 in retcodes:
        return 1
        
    else:
        return 0
        
def run_all(cfg_data_flattened):
    '''
    run_all: runs all pipeline elements in sequence
    '''

    retcodes = {} # a dictionary to hold the return codes of each pipeline element
    elements = ['fastq_filter', 'read_fuser', 'read_aligner', 'map_counts', 'map_ratios', 'map_unlink', 'enrich_plot']
    print 'Pipeline run initiated'
    print 'Running fastq_filter...'
    retcodes['fastq_filter'] = fq_filter(cfg_data_flattened)
    print 'Running read_fuser...'
    retcodes['read_fuser'] = read_fuse(cfg_data_flattened)
    print 'Running read_aligner...'
    retcodes['read_aligner'] = read_align(cfg_data_flattened)
    print 'Running map_counts...'
    retcodes['map_counts'] = counts(cfg_data_flattened)
    print 'Running map_ratios...'
    retcodes['map_ratios'] = ratios(cfg_data_flattened)
    print 'Running map_unlink...'
    retcodes['map_unlink'] = unlink(cfg_data_flattened)
    print 'Running enrich_plot...'
    retcodes['enrich_plot'] = plot(cfg_data_flattened)
    
    for key in elements: #check to see if each element succeded of failed
        if retcodes[key] == 0:
            status = 'Successful'
            
        else:
            status = 'Failed'
            
        print key + '\t' + status
    
    print 'Pipeline run complete.  Output files can be found in FX1data/output/ and FX1plots/\n'.replace('FX1', cfg_data_flattened['path'])
    

def cfg_print(cfg_order, cfg_data): 
    '''
    cfg_print: a function for printing enrich configuration files on screen
    '''
        
    try:    
        for item in cfg_order:
            print item
            elements = cfg_data[item].keys()
            elements.sort()
            
            for cfg in elements:
                print cfg + '\t' + cfg_data[item][cfg]
            
            print '\n'
                
    except:
        print 'Please load a configuration file'

def project_directory_tool(project_directory, example = 'n'): 
    '''
    project_directory_tool: a function that creates a enrich project directory, or if given a pre-existing directory, validates that it is a proper enrich project directory 
    '''
        
    if project_directory[-1] != '/': #append a / to the path if it does not exist
        project_directory = project_directory + '/'
    
    if project_directory[0] != '/': #check to see if an absolute path has been provided
        print 'Error: a complete absolute path (e.g. /path/to/project/directory/) must be specified'
        return(1)
        
    if os.path.exists(project_directory):
        dirlist = os.listdir(project_directory)
        
        if 'plots' in dirlist and 'data' in dirlist and 'input' in dirlist and 'log' in dirlist: #check to see that all appropriate directories exist, note that this function does not check for subdirectories or the enrich_default_config file 
            print 'Project directory exists and appears to be valid\n\n'
            return(0)
            
        else:
            sys.exit('Error: project directory is missing critical elements')
            return(1)
            
    else:
        if example == 'y':
            choice = 'y'
        
        else:
            choice = enrich_config.ainput('Would you like to create a new project at FX1 (y/n)? '.replace('FX1', project_directory)).getYN()
        
        if choice == 'y':
            try:
                os.mkdir(project_directory)
                os.mkdir(project_directory + 'data/')
                os.mkdir(project_directory + 'data/output/')
                os.mkdir(project_directory + 'data/tmp/')
                os.mkdir(project_directory + 'data/raw/')
                os.mkdir(project_directory + 'input/')
                os.mkdir(project_directory + 'log/')
                os.mkdir(project_directory + 'plots/')
            
            except:
                print 'Error: could not create project directory'
                return(1)
               
            try:
                from pkg_resources import resource_filename #pkg_resources is a setuptools API that allows you to interact with resources inside a package (such as a data file)
                from pkg_resources import resource_string #pkg_resources is a setuptools API that allows you to interact with resources inside a package (such as a data file)

                shutil.copy(resource_filename(__name__, 'enrich_default_config'), project_directory + 'input/enrich_default_config') #copy the default config file in place

            except:
                print 'Error: could not copy default configuration file'
            
            if example == 'y':
                try: #copy the example files in place
                    shutil.copy(resource_filename(__name__, 'example/unsel_example_F'), project_directory + 'data/raw/') 
                    shutil.copy(resource_filename(__name__, 'example/unsel_example_R'), project_directory + 'data/raw/') 
                    shutil.copy(resource_filename(__name__, 'example/unsel_example_index'), project_directory + 'data/raw/') 
                    shutil.copy(resource_filename(__name__, 'example/sel_example_F'), project_directory + 'data/raw/') 
                    shutil.copy(resource_filename(__name__, 'example/sel_example_R'), project_directory + 'data/raw/') 
                    shutil.copy(resource_filename(__name__, 'example/sel_example_index'), project_directory + 'data/raw/') 
                    
                    #install the example config files, updating the project directory element to the specified location of the example directory
                    outfile = open(project_directory + 'input/example_local_config', 'w')
                    print >> outfile, resource_string(__name__, 'example/example_local_config').replace('FX_PATH', project_directory) 
                    outfile.close()
                    
                    outfile = open(project_directory + 'input/example_SGE_config', 'w')
                    print >> outfile, resource_string(__name__, 'example/example_SGE_config').replace('FX_PATH', project_directory) 
                    
                    print 'Example project directory created\n'
                    
                except:
                    print 'Error: could not copy example files to example project directory'
                    
            if example != 'y':
                print 'Project directory created... please put your input FASTQ files in FX1data/raw/\n\n'.replace('FX1', project_directory)
            
            return(0)
                    
if __name__ == '__main__':
    main()
