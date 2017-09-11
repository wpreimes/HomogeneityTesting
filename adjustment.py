# -*- coding: utf-8 -*-
"""
Created on Aug 18 12:27 2017

@author: wpreimes
"""
from datetime import datetime
import numpy as np
import scipy.stats as stats
import pandas as pd

def compare_B(B, adjust_param, tolerance = 0.05):
    p1_B1_after = B[0,0] #add
    p1_B2_after = B[1,0] #mult

    p2_B1_after = B[0,1] #add
    p2_B2_after = B[1,1] #mult

    if adjust_param == 'both':
        if not np.isclose(p1_B1_after, p2_B1_after, atol=tolerance) or \
                not np.isclose(p1_B2_after, p2_B2_after, atol=tolerance):
            return False
        else:
            return True
    if adjust_param == 'add':
        if not np.isclose(p1_B1_after, p2_B1_after, atol=tolerance):
            return False
        else:
            return True
    if adjust_param == 'mult':
        if not np.isclose(p1_B2_after, p2_B2_after, atol=tolerance):
            return False
        else:
            return True

def adjustment(data_daily, B, breaktime,
               adjust_param = 'both', adjust_part='first', return_part=None):
    '''
    :param data: pd.DataFrame
        from Homogeneity Testing, must not contain nan
    :param B: np.matrix
    :param breaktime: string or datetime
    :param adjust_param: 'add' for adjusting only additive bias, 'mult' for adjusting only multiplicative bias
                         'both' for adjusting both model parameters
    :param adjust_part: string
        'first' to adjust values BEFORE break time
        'last' to adjust values AFTER break time
    :param return_part: string
        'all' to return adjusted dataframe over whole time frame
        'first' to return only data before breaktime
        'last' to return only data after breaktime
    :return: pd.DataFrame, dict
    '''
    if not return_part:
        return_part = adjust_part

    dataframe = data_daily.copy()

    B = np.matrix([[B['b1'][0], B['b2'][0]],
                   [B['b1'][1], B['b2'][1]]])

    part1 = dataframe[['testdata']][:breaktime]
    part2 = dataframe[['testdata']][breaktime + pd.DateOffset(1):]

    if adjust_part == 'last':
        # Perform the linear adjustment for testdata after the breaktime

        cc = B[0, 0] / B[0, 1]  # k1/k2
        dd = B[1, 0] - cc * B[1, 1] # d1 - cc*d2


        if adjust_param == 'd':
            adjusted_part2 = (cc * part2 + dd) / cc
        else:
            adjusted_part2 = cc * part2 + dd

        dataframe['adjusted'] = pd.concat([part1, adjusted_part2])

    elif adjust_part == 'first':
        # Perform the linear adjustment for testdata before the breaktime

        # B =   | k1  k2|
        #       | d1  d2 |

        cc = B[0, 1] / B[0, 0]  # k2/k1
        dd = B[1, 1] - cc * B[1, 0] # d2 - cc*d1

        if adjust_param == 'd':
            meandiff = part2.mean() - part1.mean()
            adjusted_part1 = part1 + meandiff
        else:
            adjusted_part1 = cc * part1 + dd

        dataframe['adjusted'] = pd.concat([adjusted_part1, part2])
    else:
        raise Exception("select 'first' or 'last' for part to adjust")

    adj_setting = {'slope': cc, 'intercept': dd}

    if return_part == 'all':
        return dataframe['adjusted'], adj_setting
    elif return_part == 'last':
        return dataframe['adjusted'][breaktime+pd.DateOffset(1):], adj_setting
    elif return_part == 'first':
        return dataframe['adjusted'][:breaktime], adj_setting
    else:
        raise Exception("Select 'first' , 'last' or 'all' for part to return")




def adjustment_params(data, breaktime):
    '''
    Model Data before/after breaktime via linear model, return Model paramters
    :param data: pd.DataFrame
    :param breaktime: datetime.datetime or str
    :return: np.matrix, dict
    '''
    dataframe = data.copy()
    dataframe = dataframe.dropna()
    if isinstance(breaktime, str):
        breaktime = datetime.strptime(breaktime, '%Y-%m-%d')


    i1 = dataframe[:breaktime]
    i2 = dataframe[breaktime+pd.DateOffset(1):]

    # TODO: What happens if timeframe before / after breaktime is not equally long

    # if i1.index.size != i2.index.size:
    #     if i1.index.size < i2.index.size:
    #         i2 = i2[:i1.index.size]
    #     elif i1.index.size > i2.index.size:
    #         i1 = i1[-1 * i2.index.size:]

    B = {'b1':None, 'b2':None}
    rval = []
    pval = []
    for i,data in enumerate([i1, i2], start=1):

        testdata = data['testdata'].values
        refdata = data['refdata'].values

        w = -1 * testdata # using approximate solutions 0 and 0 for k and d
        # w = refdata - testdata # using approximate solutions 1 and 0 for k and d --> different B
        l = -w
        A = np.stack(( data['refdata'].values, np.ones(data['refdata'].values.size)))
        b = np.dot(np.dot(np.linalg.inv(np.dot(A, np.transpose(A))), A), l)

        B['b%i' %i] = b

        r, p = stats.pearsonr(refdata, testdata)
        rval.append(r)
        pval.append(p)

    regtest = {'R': rval, 'p': pval}

    #TODO: activate tests

    if any(r < 0 for r in rval):
        raise Exception('2: negative correlation found (%s)' % str(rval))
    if any(p > 0.1 for p in pval): # todo: war 0.05
        raise Exception('3: positive correlation not significant (%s)' % str(pval))

    return B, regtest



