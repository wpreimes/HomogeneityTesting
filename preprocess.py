# -*- coding: utf-8 -*-
"""
Created on Tue Aug 09 10:12:54 2016

@author: wpreimes
"""

import numpy as np
import pandas as pd


def check_continuity(tmp_ts, max_gap=5):
    ts = tmp_ts.copy()

    istart = []
    iend = []
    lens = []
    oklens = []

    in_seg = False
    nan_cnt = 0

    for ind in ts.index:

        if not np.isnan(ts[ind]):
            nan_cnt = 0
            not_nan = ind
            if in_seg == False:
                istart.append(np.where(ts.index == ind)[0][0].astype('int'))
                in_seg = True

        else:
            nan_cnt += 1
            if (in_seg == True) & (nan_cnt == max_gap):
                iend.append(np.where(ts.index == not_nan)[0][0].astype('int'))
                in_seg = False

    if len(iend) < len(istart):
        iend.append(len(ts) - 1)

    for ind in np.arange(len(istart)):
        seg = ts.iloc[istart[ind]:iend[ind] + 1]
        lens.append(len(seg))
        oklens.append(len(seg[~np.isnan(seg)]))

    return istart, iend, lens, oklens


def add_nan(dataframe):
    '''
    Add nan for missing values of a time series saved as pandas dataframe
    
    Parameters
    -------
    dataframe: pandas dataframe with datetime as index
    '''
    idx = pd.date_range(dataframe.index[0], dataframe.index[-1])
    dataframe.index = pd.DatetimeIndex(dataframe.index)
    dataframe = dataframe.reindex(idx, fill_value=np.nan)
    return dataframe


def fill_gaps(tmp_data, max_gap=10):
    from mathnew.smoothn import smoothn

    data = tmp_data.copy()

    istart, iend, lens, oklens = check_continuity(data, max_gap=max_gap)

    for ind in np.arange(len(istart)):
        smoothed = smoothn(data.iloc[istart[ind]:iend[ind] + 1].values.copy(), s=0)
        # DF_Time['sm_nogaps'] = smoothed
        data.iloc[istart[ind]:iend[ind] + 1] = smoothed
        data[data < 0] = 0

    return data


def confirm(prompt=None, resp=False):
    '''
    From: http://code.activestate.com/recipes/541096-prompt-the-user-for-confirmation/
    prompts for yes or no response from the user. Returns True for yes and
    False for no.
    '''

    if prompt is None:
        prompt = 'Confirm'

    if resp:
        prompt = '%s [%s]|%s: ' % (prompt, 'y', 'n')
    else:
        prompt = '%s [%s]|%s: ' % (prompt, 'n', 'y')

    while True:
        ans = raw_input(prompt)
        if not ans:
            return resp
        if ans not in ['y', 'Y', 'n', 'N']:
            print 'please enter y or n.'
            continue
        if ans == 'y' or ans == 'Y':
            return True
        if ans == 'n' or ans == 'N':
            return False


def calc_subseries(data):
    '''
    Find series of consecutive data in a dataseries (not interrupted by nans)
    
    Return:
    subseries: List of consecutive data series, not containing nan
        
    startend: Information on the start and end position of the subseries
    in the original data series
    '''
    data = np.array(data)
    datamask = []
    datastart = np.NaN
    for index, value in enumerate(data):
        if np.isnan(value) == False and np.isnan(datastart) == True:
            datastart = index
        if np.isnan(value) == True and np.isnan(datastart) == False:
            datamask.append([datastart, index - 1])
            datastart = np.NaN
    if np.isnan(datastart) == False:
        datamask.append([datastart, data.size - 1])

    startend = []
    subseries = []
    for subset in datamask:
        startend.append(np.array([subset[0], subset[1]]))
        subseries.append(data[subset[0]:subset[1]])

    return subseries, startend



