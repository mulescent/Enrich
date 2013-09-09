#!/usr/bin/env python
'''
qvalue: A python port of the qvalue package for R
'''

__author__ = "Wayne A. Gerard"
__copyright__ = "Copyright 2011"
__credits__ = ["Douglas M Fowler", "Carlos L. Araya"]
__license__ = "FreeBSD"
__version__ = "0.2"
__maintainer__ = "Douglas M. Fowler"
__email__ = "dfowler@uw.edu"


try:
    from scipy import interpolate
    from scipy.interpolate import interp1d
    from random import choice
except:
    print 'ERROR: Scipy could not be imported. Is scipy installed?'

def qvalue_rank(x):
    '''
    Ranking function which returns number of observations less than or equal
    '''
    idx = range(len(x))
    idx.sort(lambda y,z: cmp(x[y], x[z]))
    factor_dict = {}
    index = 0
    sorted_list = sorted(x)

    for item in sorted_list:
        if item not in factor_dict:
            factor_dict[item] = index
            index += 1
    bin = []    
    [bin.append(factor_dict[item]) for item in x]
    tbl = [bin.count(index) for index in range(0, max(bin)+1)]
    cs = map(lambda x: sum(tbl[0:x+1]), range(0,len(tbl))) # x+1 because python slices are exclusive on the end    
    new_tbl = []
    [new_tbl.extend([cs[index]] * tbl[index]) for index in xrange(len(cs))] 
    
    final_table = []
    max_val = max(idx)
    for index in xrange(max_val+1):
        final_table.append(float('NaN'))

    for index in xrange(len(new_tbl)):
        final_table[idx[index]] = new_tbl[index]

    return final_table
    

def qvalue(p, lambda_val=None, pi0_method='smoother', fdr_level=None, smooth_df=3, smooth_log_pi0=False, robust=False):
    '''
    Necessary arguments:
    p: a vector of p-values

    Optional arguments:
    fdr_level: a level at which to control the FDR
    lambda_val: the value of the tuning parameter to estimate pi0
    pi0_method: either "smoother" or "bootstrap", the method for automatically choosing tuning parameter in the estimation of pi0, the proportion of true null hypotheses
    robust: an indicator of whether it is desired to make the estimate more robust for small p-values and a direct finite sample estimate of pFDR
    smooth_df: degrees of freedom to use in smoother
    smooth_log_pi0: should smoothing be done on a log scale?
    '''

    if lambda_val == None:
        lambda_val = xrange(0, 95, 5)
        lambda_val = [float(x) / 100.0 for x in lambda_val]

    # error checking

    if (min(p) < 0 or max(p) > 1):
        print 'Error: p-values not in valid range.'
        return 0
    if (len(lambda_val) > 1 and len(lambda_val) < 4):
        print 'Error: if length of lambda greater than 1, you need at least 4 values.'
        return 0
    if (len(lambda_val) > 1 and (min(lambda_val) < 0 or max(lambda_val) >= 1)):
        print 'Error: Lambda must be within [0, 1).'
    m = len(p)


    # various ways to estiamte pi0
    if (len(lambda_val) == 1):
        if (lambda_val < 0 or lambda_val >= 1):
            print 'Error: Lambda must be within [0, 1).'
        
        pi0 = float(len(filter(lambda x: x >= lambda_val, p))) / float(len(p))
        pi0 = pi0 / (1 - lambda_val)
        pi0 = min(pi0, 1)
    
        
    else:
        pi0 = [0] * len(lambda_val)
        
        for index in xrange(len(lambda_val)):
            pi0[index] = float(len(filter(lambda x: x >= lambda_val[index], p))) / float(len(p))
            pi0[index] = pi0[index] / (1 - lambda_val[index])
        
        if (pi0_method == 'smoother'):
            if (smooth_log_pi0):
                pi0 = Math.log(pi0)

            spi0 = interpolate.InterpolatedUnivariateSpline(lambda_val,pi0)
            pi0 = spi0(max(lambda_val))
            
            if (smooth_log_pi0):
                pi0 = Math.exp(pi0)
                
            pi0 = min(pi0, 1)
            
        elif (pi0_method == 'bootstrap'):

            minpi0 = min(pi0) 
            mse = [0] * len(lambda_val)
            pi0_boot = [0] * len(lambda_val)
            for i in xrange(1, 101):
                p_boot = [choice(p) for k in xrange(m)]
                for j in xrange(len(lambda_val)):
                    pi0_boot[j] = float(len(filter(lambda x: x >= lambda_val[j], p))) / float(len(p))
                
                for index in xrange(len(mse)):
                    mse[index] = float(mse[index]) + ((float(pi0_boot[index]) - float(minpi0)) ** 2)

            mse_indices = [index for index in xrange(len(mse)) if mse[index] == min(mse)]

            pi0_real = []
            for index in mse_indices:
                pi0_real.append(pi0[index])

            pi0 = min(pi0_real)
            pi0 = min(pi0, 1)

        else:
            print 'pi0_method must be on of "smoother" or "bootstrap"'
        
        if (pi0 <= 0):
            print 'The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.'

        if ((fdr_level != None) and (fdr_level <= 0 or fdr_level > 1)):
            print 'ERROR: fdr_level must be within (0, 1].'
    
        # The estimated q-values calculated here
        u = range(len(p))
        u.sort(lambda x,y: cmp(p[x], p[y]))

        v = qvalue_rank(p)
        qvalue = []

        for i in xrange(len(p)):
            qvalue.append(pi0 * m * p[i] / v[i])

        qvalue[u[m-1]] = min(qvalue[u[m-1]], 1)

        for index in xrange(m-1):
            i = -index-1
            nextindex = min(i+1, -1)
            qvalue[u[i]] = min(qvalue[u[i]], qvalue[u[nextindex]], 1)

        # The results are returned
        if (fdr_level != None):
            sig_values = filter(lambda x: x <= fdr_level, qvalue)
            retval = {'pi0': pi0, 'qvalues': qvalue, 'pvalues': p, 'fdr_level': fdr_level,
                      'significant': sig_values, 'lambda' : lambda_val}
        else:
            retval = {'pi0': pi0, 'qvalues': qvalue, 'pvalues': p, 'lambda' : lambda_val}
        return retval
