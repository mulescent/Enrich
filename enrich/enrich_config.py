#!/usr/bin/env python
'''
enrich_config: creates and edits enrich configuration files

This module provides an interactive interface for generating and editing enrich xml configuration files. Configuration information is stored as a dictionary and then written by the enrich_config_write function. Some values are automatically generated based on previous values. The default configuration file is used both to validate files that are being edited (e.g. to check that no configuration elements are missing) and to provide sensible defaults.
'''

import os, re, sys #import standard modules
import enrich_xml_parser, enrich #import project-specific modules 

__author__ = "Douglas M. Fowler"
__copyright__ = "Copyright 2011"
__credits__ = ["Douglas M Fowler", "Carlos L. Araya"]
__license__ = "FreeBSD"
__version__ = "0.2"
__maintainer__ = "Douglas M. Fowler"
__email__ = "dfowler@uw.edu"
    
def main(mode, project_directory = ''):
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
    
    #set up the configuration file and get general information about the run
    config = {} #initialize a dictionary to hold all configuration values
    section_config = {} #initialize a dictionary to hold configuration values within this section
    section_name = 'general_info'

    #find out if we are modifying an existing file or generating a new one
    new = ainput('Would you like to modify an existing configuration file (y/n)? ').getYN()
    cwd = os.getcwd().rstrip('py') #get current working directory, to be used as a guess
    if mode == 'main':
        while True:
            section_config['path'] = ainput('Exact path to project directory (e.g. /path/to/my/project/): ', 1, cwd).getString()
            
            try:
                retcode = enrich.project_directory_tool(section_config['path'])
                
                if retcode == 0:
                    if section_config['path'][-1] != '/': #append a / to the path if it does not exist
                        section_config['path'] = section_config['path'] + '/'
                    break
                    
            except:
                print 'Error: you must enter a path to a project directory'
            
    else:
        section_config['path'] = project_directory
        
    if new == 'y':
        section_config['filename'] = ainput('Enter the name of the configuration file to be modified (found in the input directory): ', 1).getString()
        
        try:
            existing_cfg = enrich_xml_parser.main(section_config['path'] + 'input/' + section_config['filename'])
        
        except: 
            print 'Error: specified configuration file does not exist or is not properly formatted'
            return(1)
            
        try:
            default_cfg = enrich_xml_parser.main(section_config['path'] + 'input/enrich_default_config')        
        
        except:
            print 'Error: default configuration file (enrich_default_config) does not exist'
            return(1)
            
        if default_cfg.keys() != existing_cfg.keys():
            sys.exit('Error: the specified configuration file is incomplete, please create a new configuration file')
            
    else:
        print 'Initializing new configuration file'     
        section_config['filename'] = ainput('What would you like to name the configuration file? ', 1).getString()
        #load the default configuration values
        existing_cfg = enrich_xml_parser.main(section_config['path'] + 'input/enrich_default_config')       
        
    section_config['run_name'] = ainput('Run name: ', 1, existing_cfg[section_name]['run_name']).getString()
    section_config['local'] = ainput('Will this run be conducted locally (L) or in an SGE environment (SGE)? ', 1).getString()
    
    while True: #this loop ensures that we get a valid, translatable DNA sequence
        section_config['wtDNA'] = ainput('Final wild type sequence (e.g. after read pair fusion):', 1, existing_cfg[section_name]['wtDNA']).getSEQ()
        
        try:
            translated = ''
            for i in xrange(0, len(section_config['wtDNA']), 3):
                    translated += codon_table[section_config['wtDNA'][i:i+3]]
            
            section_config['wtPRO'] = translated
            break 
            
        except:
            print 'Error: a valid, translatable DNA sequence was not entered'
            
    section_config['input_read1_filename'] = ainput('Input Read 1 FASTQ file (this file must reside in data/raw). If none exists, enter NONE. ', 1, existing_cfg[section_name]['input_read1_filename']).getString() #this is the forward (e.g. sense) input library read
    section_config['input_read2_filename'] = ainput('Input Read 2 FASTQ file (this file must reside in data/raw). If none exists, enter NONE. ', 1, existing_cfg[section_name]['input_read2_filename']).getString() #this is the reverse input library read
    section_config['sel_read1_filename'] = ainput('Selected Read 1 FASTQ file (this file must reside in data/raw). If none exists, enter NONE. ', 1, existing_cfg[section_name]['sel_read1_filename']).getString() #this is the forward (e.g. sense) selected library read
    section_config['sel_read2_filename'] = ainput('Selected Read 2 FASTQ file (this file must reside in data/raw). If none exists, enter NONE. ', 1, existing_cfg[section_name]['sel_read2_filename']).getString() #this is the reverse selected library read
    
    output_file = section_config['path'] + 'input/' + section_config['filename'] #save the output file and path since we will need it after we have re-initialized config
    f_output = open(output_file, 'w') #define the output file and print some preliminary xml-related lines
    print >> f_output, '<?xml version = "2.0" ?>'
    print >> f_output, '<enrich_config root_element="TRUE">'
    
    enrich_config_write(f_output, section_config, section_name)
    config.update(section_config) #add the values from this section to the master dictionary
            
    #get fastq_filter information
    section_config = {} #initialize adictionary to hold configuration values for this section
    section_name = 'fastq_filter'

    section_config['run_filter'] = ainput('Do the input FASTQ files need to be filtered (y/n)? ', 1, existing_cfg[section_name]['run_filter']).getYN()
        
    if section_config['run_filter'] == 'y':
        section_config['input_index_file'] = ainput('If there is an index read for the input library, input the filename it now (this file must reside in data/raw).  If no index read was acquired, enter NONE.  If NONE is entered, the first 20 bases of the wild type sequence will be used to filter read 1 and the last 20 bases of the wild type sequence will be used to filter read 2. ', 1, existing_cfg[section_name]['input_index_file']).getString()
        
        if section_config['input_index_file'] != 'NONE':
            section_config['sel_index_file'] = ainput('Enter the index read file for the selected library: ', 1, existing_cfg[section_name]['sel_index_file']).getString()
            section_config['index_sequence'] = ainput('Index sequence: ', 1, existing_cfg[section_name]['index_sequence']).getSEQ()
            section_config['index_mode'] = 'index'
        
        if section_config['input_index_file'] == 'NONE':
            section_config['sel_index_file'] = 'NONE'
            section_config['index_mode'] = ainput('Should Ns be counted as mismatches (y/n)? ', 1, existing_cfg[section_name]['index_mode']).getYN()
            section_config['index_sequence'] = 'NA'
            
            if section_config['index_mode'] == 'y':
                section_config['index_mode'] = 'N_include'
                
            elif section_config['index_mode'] == 'n':
                section_config['index_mode'] = 'N_exclude'
            
            else:
                sys.exit('Error: choose a valid index mode')
            
        section_config['index_mismatch_threshold'] = ainput('Enter the mismatch threshold for the filter.  If you have an index read, this will be the number of mismatches allowed between the index read and its wild type sequence.  If you do not have an index read, this will be the number of mismatcehs allowed between the read and its wild type sequence ', 1, existing_cfg[section_name]['index_mismatch_threshold']).getInteger()
    
    else:
        section_config['input_index_file'] = 'NA'
        section_config['sel_index_file'] = 'NA'
        section_config['index_sequence'] = 'NA'
        section_config['index_mode'] = 'NA'
        section_config['index_mismatch_threshold'] = 'NA'
        
    enrich_config_write(f_output, section_config, 'fastq_filter')
    config.update(section_config) #add the values from this section to the master dictionary
    
    #get read_fuser information
    section_config = {}
    section_name = 'read_fuser'
    
    if config['index_mode'] == 'index':
        section_config['input_read1_filename_fuser'] = config['input_read1_filename'].rstrip('.fq') + '_' + config['index_mode'] + '_filtered.fq'
        section_config['input_read2_filename_fuser'] = config['input_read2_filename'].rstrip('.fq') + '_' + config['index_mode'] + '_filtered.fq'
        section_config['sel_read1_filename_fuser'] = config['sel_read1_filename'].rstrip('.fq') + '_' + config['index_mode'] + '_filtered.fq'
        section_config['sel_read2_filename_fuser'] = config['sel_read2_filename'].rstrip('.fq') + '_' + config['index_mode'] + '_filtered.fq'
        
    elif config['index_mode'] == 'N_include':
        section_config['input_read1_filename_fuser'] = config['input_read1_filename'].rstrip('.fq') + '_' + config['index_mode'] + '_filtered.fq'
        section_config['input_read2_filename_fuser'] = config['input_read2_filename'].rstrip('.fq') + '_' + config['index_mode'] + '_filtered.fq'
        section_config['sel_read1_filename_fuser'] = config['sel_read1_filename'].rstrip('.fq') + '_' + config['index_mode'] + '_filtered.fq'
        section_config['sel_read2_filename_fuser'] = config['sel_read2_filename'].rstrip('.fq') + '_' + config['index_mode'] + '_filtered.fq'
        
    elif config['index_mode'] == 'N_exclude':
        section_config['input_read1_filename_fuser'] = config['input_read1_filename'].rstrip('.fq') + '_' + config['index_mode'] + '_filtered.fq'
        section_config['input_read2_filename_fuser'] = config['input_read2_filename'].rstrip('.fq') + '_' + config['index_mode'] + '_filtered.fq'
        section_config['sel_read1_filename_fuser'] = config['sel_read1_filename'].rstrip('.fq') + '_' + config['index_mode'] + '_filtered.fq'
        section_config['sel_read2_filename_fuser'] = config['sel_read2_filename'].rstrip('.fq') + '_' + config['index_mode'] + '_filtered.fq'

    else:
        section_config['input_read1_filename_fuser'] = config['input_read1_filename']
        section_config['input_read2_filename_fuser'] = config['input_read2_filename']
        section_config['sel_read1_filename_fuser'] = config['sel_read1_filename']
        section_config['sel_read2_filename_fuser'] = config['sel_read2_filename']
        
    if config['input_read1_filename'] == 'NONE':
        section_config['fuser_mode'] = 'R2'
    
    elif config['input_read2_filename'] == 'NONE':
        section_config['fuser_mode'] = 'R1'
        
    else:
        section_config['fuser_mode'] = 'B'
        
    if config['input_read1_filename'] == 'NONE' and config['input_read2_filename'] == 'NONE':
        sys.exit('Error: both filenames cannot be NONE')
        
    section_config['read1_overlap_start'] = ainput('Read 1 overlap start: ', 1, existing_cfg[section_name]['read1_overlap_start']).getInteger()
    section_config['read1_overlap_end'] = ainput('Read 1 overlap end: ', 1, existing_cfg[section_name]['read1_overlap_end']).getInteger()
    section_config['read2_overlap_start'] = ainput('Read 2 overlap start: ', 1, existing_cfg[section_name]['read2_overlap_start']).getInteger()
    section_config['read2_overlap_end'] = ainput('Read 2 overlap end: ', 1, existing_cfg[section_name]['read2_overlap_end']).getInteger()
    section_config['include_nonoverlap_region'] = ainput('Should the non-overlapping 5 prime ends of read 1 and read 2 be appended to the overlapped region (y/n)? ', 1, existing_cfg[section_name]['include_nonoverlap_region']).getYN()
    section_config['run_aligner'] = ainput('Should a Needleman-Wunsch aligner be used to identify gapped read 1/read 2 pairs (warning: this increases analysis time significantly) (y/n)? ', 1, existing_cfg[section_name]['run_aligner']).getYN()
    
    if section_config['run_aligner'] == 'y':
        section_config['paired_mismatch_threshold'] = ainput('What is the number of mismatches in a read pair above which a Needleman-Wunsch alignment should be performed? ', 1, existing_cfg[section_name]['paired_mismatch_threshold']).getInteger()
        
    else:
        section_config['paired_mismatch_threshold'] = 0
    
    enrich_config_write(f_output, section_config, 'read_fuser')
    config.update(section_config) #add the values from this section to the master dictionary
    
    #get read_aligner information
    section_config = {}
    section_name = 'read_aligner'
    
    if config['fuser_mode'] == 'B' or config['fuser_mode'] == 'R1':
        section_config['input_filename_aligner'] = config['input_read1_filename_fuser'].rstrip('.fq') + '_' + config['fuser_mode']
        section_config['sel_filename_aligner'] = config['sel_read1_filename_fuser'].rstrip('.fq') + '_' + config['fuser_mode']
    
    elif config['fuser_mode'] == 'R2':
        section_config['input_filename_aligner'] = config['input_read2_filename_fuser'].rstrip('.fq') + '_' + config['fuser_mode']
        section_config['sel_filename_aligner'] = config['sel_read2_filename_fuser'].rstrip('.fq') + '_' + config['fuser_mode']
        
    if config['run_aligner'] == 'y':
        section_config['gap_max'] = ainput('Maximum number of paired read alignment gaps allowed in the finished sequence: ', 1, existing_cfg[section_name]['gap_max']).getInteger()
        
    else:
        section_config['gap_max'] = 0
        
    section_config['unresolvable_max'] = ainput('Maximum number of unresolved bases (e.g. read pair disagreements with equal quality scores) allowed in a finished sequence: ', 1, existing_cfg[section_name]['unresolvable_max']).getInteger()
    section_config['maximum_mutation_run'] = ainput('Maximum number of consecutive mutations allowed in a finished sequence: ', 1, existing_cfg[section_name]['maximum_mutation_run']).getInteger()
    section_config['avg_quality'] = ainput('Minimum read average quality score for each finished sequence: ', 1, existing_cfg[section_name]['avg_quality']).getInteger()
    section_config['chaste'] = ainput('Reject reads that fail the Illumina chastity filter (y/n)? ', 1, existing_cfg[section_name]['chaste']).getYN()
    section_config['Ncount_max'] = ainput('Maximum number of "N" bases allowed in a finished sequence: ', 1, existing_cfg[section_name]['Ncount_max']).getInteger()
    
    enrich_config_write(f_output, section_config, 'read_aligner')
    config.update(section_config) #add the values from this section to the master dictionary
    
    #get map_counts information
    section_config = {}
    section_name = 'map_counts'
    
    section_config['input_filename_map_counts'] = ','.join([config['input_filename_aligner'] + '_PRO_qc', config['input_filename_aligner'] + '_DNA_qc'])
    section_config['sel_filename_map_counts'] = ','.join([config['sel_filename_aligner'] + '_PRO_qc', config['sel_filename_aligner'] + '_DNA_qc'])

    enrich_config_write(f_output, section_config, 'map_counts')
    config.update(section_config) #add the values from this section to the master dictionary
    
    #get map_ratios information
    section_config = {}
    section_name = 'map_ratios'
    
    section_config['input_filename_map_ratios'] = ','.join(['counts_' + config['input_filename_aligner'] + '_PRO_qc', 'counts_' + config['input_filename_aligner'] + '_DNA_qc'])
    section_config['sel_filename_map_ratios'] = ','.join(['counts_' + config['sel_filename_aligner'] + '_PRO_qc', 'counts_' + config['sel_filename_aligner'] + '_DNA_qc'])

    enrich_config_write(f_output, section_config, 'map_ratios')
    config.update(section_config) #add the values from this section to the master dictionary

    #get map_unlink information
    section_config = {}
    section_name = 'map_unlink'
    section_config['unlink_modes'] = ainput('Modes in which to run map_unlink, separate multiple values by commas (valid options = wild_counts,counts): ', 1, existing_cfg[section_name]['unlink_modes']).getString()
    section_config['input_filename_map_unlink'] = ','.join(['counts_' + config['input_filename_aligner'] + '_PRO_qc', 'counts_' + config['input_filename_aligner'] + '_DNA_qc', 'counts_' + config['input_filename_aligner'] + '_PRO_qc.m1', 'counts_' + config['input_filename_aligner'] + '_DNA_qc.m1', 'counts_' + config['input_filename_aligner'] + '_PRO_qc.m2', 'counts_' + config['input_filename_aligner'] + '_DNA_qc.m2'])
    section_config['sel_filename_map_unlink'] = ','.join(['counts_' + config['sel_filename_aligner'] + '_PRO_qc', 'counts_' + config['sel_filename_aligner'] + '_DNA_qc', 'counts_' + config['sel_filename_aligner'] + '_PRO_qc.m1', 'counts_' + config['sel_filename_aligner'] + '_DNA_qc.m1', 'counts_' + config['sel_filename_aligner'] + '_PRO_qc.m2', 'counts_' + config['sel_filename_aligner'] + '_DNA_qc.m2'])

    enrich_config_write(f_output, section_config, 'map_unlink')
    config.update(section_config) #add the values from this section to the master dictionary
    
    #get enrich_plot information
    section_config = {}
    section_name = 'enrich_plot'
    
    modes = config['unlink_modes'].split(',')
    input_files = config['input_filename_map_counts'].split(',')
    sel_files = config['sel_filename_map_counts'].split(',')
    freqplot_list = []
    
    for item in modes:
        for i in xrange(0,1):
            if 'DNA' not in input_files[i] and 'wild' not in input_files[i]:
                freqplot_list.append('unlink_' + item + '_' + input_files[i])
            
            if 'DNA' not in sel_files[i] and 'wild' not in sel_files[i]:
                freqplot_list.append('unlink_' + item + '_' + sel_files[i])

    section_config['input_filename_freqplot'] = ','.join(freqplot_list)
    section_config['input_filename_ratioplot'] = ';'.join(['unlink_counts_' + config['sel_filename_aligner'] + '_PRO_qc,' + 'unlink_counts_' + config['input_filename_aligner'] + '_PRO_qc','unlink_counts_' + config['sel_filename_aligner'] + '_PRO_qc.m1,' + 'unlink_counts_' + config['input_filename_aligner'] + '_PRO_qc.m1', 'unlink_counts_' + config['sel_filename_aligner'] + '_PRO_qc.m2,' + 'unlink_counts_' + config['input_filename_aligner'] + '_PRO_qc.m2'])

    enrich_config_write(f_output, section_config, 'enrich_plot')
    config.update(section_config) #add the values from this section to the master dictionary
    
    #finish printing the configuration file and, if the function is being called as a module return the path to the newly-created configuration file
    print >> f_output, '</enrich_config>'
    f_output.close()    
    
    if mode == 'main':
        print 'Created configuration file ' + output_file
        
    if mode != 'main':
        return(output_file)
        
class ainput:
    '''
    ainput: creates an object for acquiring console input with one of several specific properties (e.g. is an integer, is y/n, is a string, etc)
    '''

    def __init__(self, msg='', req=0, default = ''):
      ''' This will create a instance of ainput object'''
      self.data = ''   #Initialize a empty data variable
      if not msg == '':
         self.ask(msg, req, default)
    
    def ask(self, msg, req, default):
      ''' This will display the prompt and retrieve the user input.'''
      
      if default != '':
          if req == 0:
             self.data = raw_input(msg + '(' + default + ', type DEFAULT to use) ')   #Save the user input to a local object variable
          else:
             self.data = raw_input(msg + '(' + default + ', type DEFAULT to use) ')
        
              #Check to see if the default is desired (e.g. if DEFAULT was entered)
          if self.data == 'DEFAULT':
              self.data = default
             
      else:
          
        if req == 0:
         self.data = raw_input(msg)   #Save the user input to a local object variable
        
        else:
         self.data = raw_input(msg)
      
      #Verify that the information was entered and its not empty.
      if req == 1 and self.data == '':
         print 'Error: response required'
         self.ask(msg, req, default)
    
    def getString(self):
      ''' Returns the user input as String'''
      return self.data
    
    def getInteger(self):
      ''' Returns the user input as a Integer'''
      while True:
        try:
            return int(self.data)
        
        except:
            self.data = raw_input('You must enter an integer: ')
    
    def getNumber(self):
      ''' Returns the user input as a Float number'''
      while True:
        try:
            return float(self.data)
        
        except:
            self.data = raw_input('You must enter a number: ')      
    
    def getYN(self):
      '''Returns a y or n'''
      while True:
        if self.data in ['Y', 'y', 'N', 'n']:
            return(self.data.lower())
            
        else:
            self.data = raw_input('You must enter "y" or "n": ')
            
    def getSEQ(self): 
        '''returns a DNA sequence'''
        while True:
            try:
                sequence = self.data.upper()
       
                if re.search('[^AGCT]', sequence):
                    self.data = raw_input('You must enter a valid DNA sequence (i.e. ATCG only) ')
                
                else:
                   return(sequence) 
            
            except:
                 self.data = raw_input('You must enter a valid DNA sequence (i.e. ATCG only) ')
                            
def enrich_config_write(file_handle, config_dict, section_name = ''):
    '''
    enrich_config_write: this function writes a chunck of config values to an open config file
    '''

    config_keys = sorted(config_dict.keys()) #sort the keys in alphabetical order
    print >> file_handle, '<' + section_name + ' section="TRUE">'
    
    for key in config_keys:
        print >> file_handle, '<' + key + '>' + str(config_dict[key]) + '</' + key + '>'
    
    print >> file_handle, '</' + section_name + '>'
    
def enrich_config_flatten(enrich_config):
    '''
    enrich_config_flatten: this function takes a enrich XML configuration file and flattens it into a dictionary with a depth of one (e.g. it discards section information).
    '''

    flattened_config = {}
    for section in enrich_config:
        flattened_config.update(enrich_config[section])
    
    return(flattened_config)
    
if __name__ == '__main__':
    main('main')
    
