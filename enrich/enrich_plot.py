#!/usr/bin/env python
'''
enrich_plot: contains a number of functions that use matplotlib to generate visualizations of enrich data.  This module is intended to be extensible, containing all enrich-related visuals
'''

__author__ = "Douglas M. Fowler"
__copyright__ = "Copyright 2011"
__credits__ = ["Douglas M Fowler", "Carlos L. Araya"]
__license__ = "FreeBSD"
__version__ = "0.2"
__maintainer__ = "Douglas M. Fowler"
__email__ = "dfowler@uw.edu"

import sys, os, time, math

try: 
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pylab as pylab
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt

except:
    print 'Error: Matplotlib not present. Cannot create visuals' 

try:
    import numpy as np

except:
    print 'Error: Numpy not present. Cannot create visuals'
    
def freqplot(path, infile):
    '''
    freqplot: takes a mapunlink file and produces a heat map for visualizing diversity
    '''

    '''
    Basic I/O checks - making sure files exist, etc.
    '''
    try:
        # Check to make sure the input file exists
        f_infile = open(path + 'data/output/' + infile, 'U')
        f_infile.close()
        
    except:
        print 'Error: could not open input file'
        return 1

    try:
        plt_test = plt
        np_test = np.array

    except: # already warned the user about this above, so no message
        return 1
    
    f_infile = open(path + 'data/output/' + infile, 'U')
    header_line = f_infile.readline().split()

    # Make sure we can locate the position column
    try:
        pos_index = header_line.index('position')
    except:
        print 'Could not find the position index for plotting. Was the input file labeled properly?'
        return 1

    '''
    Reading in the data

    xs: the positions
    ys: the AA substitutions/PRO substitutions
    zs: frequencies
    values: a len(y) x len(x) matrix containing all the z values.
            Necessary for plotting the data as a heat map.
    '''
        
    xs = []
    ys = [index for index in xrange(len(header_line)-1)]
    zs = []

    # Grab all the data from the file and dump it into xs and zs
    for line in f_infile:
        line = line.split()
        xs.append(int(line[pos_index]))
        zs.extend([float(line[index]) for index in xrange(len(line)) if index is not pos_index])

    values = np.array(zs)
    values = values.reshape(len(xs), len(ys))
    values = values.transpose()

    '''
    Plotting the data

    Setting up graphical parameters for pyplot
    '''

    # Labels and label fonts
        
    xlab = 'Position'

    if 'PRO' in infile:
        ylab = 'Amino acid'
    elif 'DNA' in infile:
        ylab = 'Nucleotide'
    else:
        print 'Could not determine type of data present for plotting a heat map. Is the file name appropriately labeled?'
        return 1
    
    zlab = 'Frequency'

    title_font = {'fontsize': 18}

    # Tick marks - the padding is necessary so heat map boxes are not cut off

    plotys = ['']
    plotys.extend([header_line[index] for index in xrange(len(header_line)) if index != pos_index])
    plotys.append('')


    # Plotting the actual data
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(values, cmap=cm.Blues)

    # Resizing the figure to be slightly larger
    
    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0]*1.5, DefaultSize[1]*1.5) )
    
    # Set the labels for the axes
    
    cax.axes.xaxis.set_label_text(xlab, **title_font)
    cax.axes.yaxis.set_label_text(ylab, **title_font)


    # Plotting the tick marks - the padding is necessary so boxes are not cut off
    len_y = [-1]
    len_y.extend(range(len(plotys)))
    len_y.append(len(plotys))

    pylab.yticks(len_y, plotys)

    # Drawing the color bar an appropriate size
    
    colorbar_args = {'shrink': 0.62}
    pylab.colorbar(cax, **colorbar_args)

    '''
    Saving

    Save the graph out as a pdf to disk
    '''
    fileout = infile + '_diversity.pdf'
    fileout_total = path + 'plots/' + fileout
    plt.title(fileout)
    
    pylab.savefig(fileout_total, dpi=100)
    
    return 0


def ratioplot(path, infile1, infile2):
    '''
    ratioplot: takes two mapunlink files, calculates and calculates the position-averaged ratio of the frequency of mutation between them
    '''

    '''
    Basic I/O checks - making sure files exist, etc.
    '''
    
    try: #check to make sure the input file exists
        f_infile = open(path + 'data/output/' + infile1, 'U')
        f_infile.close()
        f_infile = open(path + 'data/output/' + infile2, 'U')
        f_infile.close()
        
    except:
        print 'Error: could not open input file'
        return 1

    
    try:
        plt_test = plt
        np_test = np.array

    except: # already warned the user about this above, so no message
        return 1
    
    '''
    Reading the data: Calculate the ratio of frequencies for both files.
                      Assign the ratios appropriate colors based on whether the
                      ratio is above 1 or <= 1.

    file1_values/file2_values: dictionaries mapping {'position': frequency}
    '''
    f_infile = open(path + 'data/output/' + infile1, 'U')
    file1_values = read_position_file(f_infile)
        
    f_infile = open(path + 'data/output/' + infile2, 'U')
    file2_values = read_position_file(f_infile)


    ratios = []
    colors = []
    positions = []
    for position in file1_values:
        x = float(sum(file1_values[position]))
        y = float(sum(file2_values[position]))
        #rat = (x - position) / (y - position) WHY?
        if y != 0:
            rat = x / y
        
        else:
            rat = 1 #this is highly unsatisfactory, there should be some indication that this data is missing
            
        ratios.append(rat)
        positions.append(position)
        if (rat <= 1):
            colors.append('b')
        else:
            colors.append('r')


    '''
    Plotting the data: plot the ratios
    '''
            
    fig = plt.figure()
    fig.add_subplot(111)
    
    ax = plt.plot( [positions[index] for index in xrange(len(ratios)) if ratios[index] > 1],
             [ratios[index] for index in xrange(len(ratios)) if ratios[index] > 1], 'ro')
    ax = plt.plot( [positions[index] for index in xrange(len(ratios)) if ratios[index] <= 1],
             [ratios[index] for index in xrange(len(ratios)) if ratios[index] <= 1], 'bo')


    # Set the labels and other graphical parameters
    ylab = 'Position-averaged mutation frequency in selected/input'
    xlab = 'Position'

    plt.gca().set_xlim([min(positions) - 1, max(positions) + 1])
    plt.xlabel(xlab)
    plt.ylabel(ylab)

    # Draw the vertical lines that connect the baseline to each point
    plt.vlines(positions, [1] * len(positions), ratios, color=colors, linestyles='solid', lw=2)
    plt.axhline(y=1, color='gray')

    # Resize the image to be a little larger
    DefaultSize = fig.get_size_inches()    
    fig.set_size_inches( (DefaultSize[0]*1.5, DefaultSize[1]*1.5) )

    # Draw the title manually because otherwise the default offset doesn't work
    fileout =  infile1 + "_" + infile2 + '_position_enrichment.pdf'
    fig.text(0.52, 0.97, fileout, ha="center", va="center", size="medium")

    '''
    Saving to disk
    '''
    fileout_total = path + 'plots/' + fileout

    pylab.savefig(fileout_total, dpi=100)

    return 0


def read_position_file(f_infile):
    '''
    reads a file into a dictionary mapping positions -> frequencies
    '''

    header = f_infile.readline().split()
    position_index = header.index('position')

    file_values = {}

    for line in f_infile:
        line = line.split()
        file_values[int(line[position_index])] = [float(line[index]) for index in xrange(len(line)) if index is not position_index]

    f_infile.close()
    return file_values


def read_substitution_file(f_infile):
    '''
    Reads in a file mapping 
    '''
    header = f_infile.readline().split()
    position_index = header.index('position')
    
    file_values = {'positions' : []}
    
    for index in xrange(len(header)):
        if index is not position_index:
            file_values[header[index]] = {'values': []}

    for line in f_infile:
        line = line.split()
        [file_values['positions'].append(int(line[position_index]))]
        [file_values[header[index]]['values'].append(float(line[index])) for index in xrange(len(line)) if index is not position_index]

    f_infile.close()
    return file_values


def separate_data(substitution, positions, ratio, wt_data, axis_num, cols, rows):
    '''
    Another small helper function to follow DRY principles with the AA plot below.
    '''

    values = ratio[substitution]
    xr = xrange(len(values))

    # Separate out the values based on their category
    wt_x = [positions[index] for index in xr if (values[index] == 0 and wt_data['variable'][index] == substitution)]
    wt_y = [values[index] for index in xr if (values[index] == 0 and wt_data['variable'][index] == substitution)]
    nan_x = [positions[index] for index in xr if (values[index] == 0 and wt_data['variable'][index] != substitution)]
    nan_y = [values[index] for index in xr if (values[index] == 0 and wt_data['variable'][index] != substitution)]
    rest_x = [positions[index] for index in xr if values[index] != 0]
    rest_y = [values[index] for index in xr if values[index] != 0]

    ax = plt.subplot(rows, cols, axis_num)

    header_tick_values = range(0, len(values), (len(values)/4))

    all_residue_plot_data(wt_x, wt_y, nan_x, nan_y, rest_x, rest_y) 
    adjust_axes(ax, str(substitution), header_tick_values, header_tick_values, axis_num, cols, 18)
    
def AAPlot(wt_data, ratio, axis_num, positions):

    '''
    Takes in ratios and maps how values should be plotted. Residues are
    plotted based on 3 categories: WT residues, residues with ratios of 1,
    and all other residues.

    Since this is a faux lattice plot, each substitution must be iterated through.
    '''

    rows = int(math.ceil(float(len(ratio.keys())) / 5.0))
    cols = 5

    keys = ratio.keys()
    keys = sorted(keys)
    
    for substitution in keys:
        if (substitution != '*'): # Save the stop index for last
            axis_num += 1
            # a : {0, 1, 2, 0, ... }
            separate_data(substitution, positions, ratio, wt_data, axis_num, cols, rows)


    '''
    Plotting the final substitution, the stop codon/index
    '''
    substitution = '*'
    axis_num += 1
    separate_data(substitution, positions, ratio, wt_data, axis_num, cols, rows)
    
def PosPlot(wt_data, ratio, axis_num, positions):
    
    '''
    Takes in ratios and maps how values should be plotted. Positions are
    plotted based on 3 categories: WT residues, residues with ratios of 1,
    and all other residues.

    Rather than re-read the data based on positions instead of substitutions, it makes sense to just iterate through the substitution dictionary and gather the necessary data.
    '''
    position_d = {'header' : []}

    '''
    Dealing with the stop position/codon, as well as mapping positions -> frequencies
    '''
    aa_sorted = ratio.keys()
    aa_sorted.sort()
    aa_sorted.pop(0)
    aa_sorted.append('*')
    
    for substitution in aa_sorted:
        if (substitution != '*'):
            position_d['header'].append(substitution)
            for index in xrange(len(ratio[substitution])):
                if position_d.get(index,None) == None:
                    position_d[index] = []
                position_d[index].append(ratio[substitution][index])

    position_d['header'].append('*')
    for index in xrange(len(ratio['*'])):
        if position_d.get(index, None) == None:
            position_d[index] = []
        position_d[index].append(ratio['*'][index])


    rows = int(math.ceil(float(len(position_d.keys())) / 5.0))
    cols = 5
    
    keys = [key for key in position_d if key != 'header']
    keys.sort()
    
    for position in keys:
        
        axis_num += 1
        # 1: {0.1, ...}
        
        values = position_d[position]
        xr = xrange(len(values))

        '''
        Separate categorizing from the AA plot above because we're indexing on position.
        '''
        wt_x = [positions[index] for index in xr if (values[index] == 0 and wt_data['variable'][position] == position_d['header'][index])]
        wt_y = [values[index] for index in xr if (values[index] == 0 and wt_data['variable'][position] == position_d['header'][index])]
        nan_x = [positions[index] for index in xr if (values[index] == 0 and wt_data['variable'][position] != position_d['header'][index])]
        nan_y = [values[index] for index in xr if (values[index] == 0 and wt_data['variable'][position] != position_d['header'][index])]
        rest_x = [positions[index] for index in xr if values[index] != 0]
        rest_y = [values[index] for index in xr if values[index] != 0]

        ax = plt.subplot(rows, cols, axis_num)

        all_residue_plot_data(wt_x, wt_y, nan_x, nan_y, rest_x, rest_y)
        adjust_axes(ax, str(position), range(len(position_d['header'])), position_d['header'], axis_num, cols, 13)

def adjust_axes(ax, title, header_ticks, header_list, axis_num, cols, axesfontsize):
    '''
    A small helper function that adjusts the axes on the plots
    '''
    
    frame = plt.gca()
    
    axes_font = {'fontsize':axesfontsize}
    title_font = {'fontsize':18}
    
    plt.text(0.1,3, title, **title_font)
    
    plt.axhline()
    
    if (axis_num <= cols) and (axis_num % 2 == 0):
        frame.xaxis.set_ticks_position('top')
        frame.xaxis.set_label_position('top')
        frame.xaxis.tick_top()
        ax.set_xticks(header_ticks)
        ax.set_xticklabels(header_list, **axes_font)
    else:
        frame.axes.get_xaxis().set_visible(False)
        
    if (((axis_num - 1) / cols) % 2 == 0) and axis_num % cols == 1:
        ax.set_yticks((-4,-2,0,2,4))
        ax.set_yticklabels([-4,-2,0,2,4], **axes_font)
    else:
        frame.axes.get_yaxis().set_visible(False)

def all_residue_plot_data(wt_x, wt_y, nan_x, nan_y, rest_x, rest_y):
    '''
    Another small helper function that assists with plotting the actual data.
    '''
    plt.plot(wt_x, wt_y, "rs")
    plt.plot(nan_x, nan_y, color="gray", marker="s")
    plt.plot(rest_x, rest_y, "bo")
    plt.ylim(-4, 4)
    xs = wt_x
    xs.extend(nan_x)
    xs.extend(rest_x)
    plt.xlim(min(xs) - 1, max(xs) + 1)
    
def all_residue_plot(path, infile1, infile2, wtseq, mode):
    '''
    all_residue_plot: takes two mapunlink files, calculates and calculates the ratio of the frequency of mutation between them for each mutation-position combination
    '''

    '''
    Basic I/O checks - check that files exist, etc.
    '''
    
    try:
        # Check to make sure the input file exists
        f_infile = open(path + 'data/output/' + infile1, 'U')
        f_infile.close()
        f_infile = open(path + 'data/output/' + infile2, 'U')
        f_infile.close()
        
    except:
        print 'Error: could not open input file'
        return 1


    try:
        plt_test = plt
        np_test = np.array

    except: # already warned the user about this above, so no message
        return 1
    '''
    Read in the data
    file1_values/file2_values: a dictionary maaping {'position': 'frequency'}
    wt_data: a dictionary containing the wildtype positions and AA/proteins
    '''
    
    f_infile = open(path + 'data/output/' + infile1, 'U')
    file1_values = read_substitution_file(f_infile)

    f_infile = open(path + 'data/output/' + infile2, 'U')
    file2_values = read_substitution_file(f_infile)

    # Obtain wildtype residues and set at the enrichment ratio for the wildtype sequence:
    wt_data = {'position': [], 'variable': []}
    for position in xrange(len(wtseq)):
        wt_data['position'].append(position)
        wt_data['variable'].append(wtseq[position])


    '''
    Calculating the log2 enrichment ratio for each mutation
    '''

    ratio = {}
    positions = file1_values['positions']
    
    for substitution in file1_values:
        if substitution != 'positions':
            f1p = file1_values[substitution] # Caching for faster lookup
            f2p = file2_values[substitution]

            if (ratio.get(substitution, None) == None):
                ratio[substitution] = []

            for index in xrange(len(f1p['values'])):
                if f2p['values'][index] != 0.0 and f1p['values'][index] != 0.0:
                    ratio[substitution].append(math.log(f1p['values'][index] / f2p['values'][index],2))
                else:
                    ratio[substitution].append(0)


    '''
    Setting up the graph to accomodate multiple sub plots and look similar
    to a lattice plot
    '''
                    
    axis_num = 0
    fig = pylab.figure()
    fig.subplots_adjust(wspace=0.0001, hspace=0.0001)
    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0]*2.5, DefaultSize[1]*2.5) )    

    '''
    Plot the data - see individual functions for more details
    '''
    ylab = ""
    if mode is 'AA':
        AAPlot(wt_data, ratio, axis_num, positions)
        ylab = "Position"
    elif mode is 'Pos':
        PosPlot(wt_data, ratio, axis_num, positions)
        ylab = "Amino Acid Substitution"


    # Drawing axis labels and other titles

    fileout = infile1 + "_" + infile2 + '_all_residue_by_' + mode + '.pdf'
    fileout_total = path + 'plots/' + fileout

    fig.text(0.1, 0.52, "log2(Mutation frequency in selected/input)", ha="right", va="center", size="xx-large", rotation="vertical")
    fig.text(0.52, 0.08, ylab, ha="center", va="bottom", size="xx-large")
    fig.text(0.52, 0.99, fileout, ha="center", va="top", size="xx-large")

    '''
    Saving to disk
    '''
    
    pylab.savefig(fileout_total, dpi=100)
    
    return 0
