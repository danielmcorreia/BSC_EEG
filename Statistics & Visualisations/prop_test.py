#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 17:02:24 2022

@author: Daniel
"""

from statsmodels.stats.proportion import proportions_ztest
import numpy as np

def prop_ztest(count, nobs, value, alternative='two-sided'):
    
    '''count: The number of successes (count sucesses drug1)
    nobs: The number of trials (total numb datapoints, usually 195 for 1sec, 0.8 overlap)
    value: The hypothesized population proportion (prop. drug2)
    p_vals: list to store pvals
    alternative: The alternative hypothesis (two tailed for = or =!)
    
    H0: p = p0 (population proportion is equal to hypothesized proportion p0)
    H1 (two-tailed): p ≠ p0 (population proportion is not equal to some hypothesized value p0)'''
    
    stat, pval = proportions_ztest(count, nobs, value, alternative)

    return pval


def eval_pval(pval):
    
    '''evaluates level os significance of pval for reporting puposes'''
    
    if pval < .001:
        return ('***', 'p < .001')
    if pval < .01:
        return ('**', 'p < .01')
    if pval < .05:
        return ('*', 'p < .05')
    if pval > .05:
        return ('n.s.', 'p > .05')


# Custom function to draw the diff bars

'''def label_diff(i,j,text,X,Y):
    x = 0.5
    y = 75
    dx = abs(X[i]-X[j])

    props = {'connectionstyle':'bar','arrowstyle':'-',\
                 'shrinkA':20,'shrinkB':20,'linewidth':2}
    ax.annotate(text, xy=(X[i],y+7), zorder=10)
    ax.annotate('', xy=(X[i],y), xytext=(X[j],y), arrowprops=props)

# Call the function
label_diff(0,1,'p=0.0370',np.array([0,1]),75)'''
#label_diff(0,1,'p<0.0001',ind,menMeans)
#label_diff(0,1,'p=0.0025',ind,menMeans)