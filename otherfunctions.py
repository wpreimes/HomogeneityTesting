# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 12:58:31 2016

@author: wpreimes
"""

from scipy import stats
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from netCDF4 import date2num
from scipy import io
import os
from datetime import datetime
import csv



def cci_extract(_string):
    cont = _string.upper().split('_')
    if len(cont) == 3:
        return {'prefix': cont[0], 'version': cont[1], 'type': cont[2], 'adjust': None}
    elif len(cont) == 4 and cont[3] == 'ADJUSTED':
        return {'prefix': cont[0], 'version': cont[1], 'type': cont[2], 'adjust': cont[3]}
    else:
        raise Exception('Unknown Product')

def cci_string_combine(info):
    if info['adjust']:
        return  "_".join([info.get(key) for key in ['prefix', 'version', 'type','adjust']])
    else:
        return "_".join([info.get(key) for key in ['prefix', 'version', 'type']])

def split(el, n):
    '''
    Split list of cells in n approx. equal parts for multiprocessing
    :param el: list of elements to split
    :param n: number of lists to split input up into
    :return: list
    '''
    k, m = divmod(len(el), n)
    return (el[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))


def create_workfolder(path):
    # type: (str) -> str
    i = 1
    while os.path.exists(os.path.join(path, 'v' + str(i))):
        i += 1
    else:
        os.makedirs(os.path.join(path, 'v' + str(i)))

    workfolder = os.path.join(path, 'v' + str(i))
    print('Create workfolder: %s' % str(workfolder))

    return workfolder


def regress(data, refdata_col_name):
    '''
    Perform regression of column refdata on column testdata
    '''
    dataframe = data.copy()
    out = dataframe[refdata_col_name]
    dataframe = dataframe.dropna()
    R, pval = stats.pearsonr(dataframe[refdata_col_name], dataframe.testdata) # Correlation between Refdata and Testdata


    if R < 0 or np.isnan(R):
        ress = [np.nan]
        return out, R, pval, ress

    testdata = dataframe.testdata.values
    refdata = dataframe[refdata_col_name].values
    refdata_ones = np.vstack([refdata, np.ones(len(refdata))]).T

    ress = np.linalg.lstsq(refdata_ones, testdata)[0][::-1]
    dataframe['ones'] = 1
    xm = np.matrix(dataframe.as_matrix(columns=['ones', refdata_col_name]))
    out = np.dot(xm, np.matrix(ress).transpose())

    return pd.Series(index=dataframe.index, data=np.squeeze(np.asarray(out))), R, pval, ress


def datetime2matlabdn(dt):
    ord = dt.toordinal()
    mdn = dt + timedelta(days=366)
    frac = (dt - datetime(dt.year, dt.month, dt.day, 0, 0, 0)).seconds / (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac


def dates_to_num(dates):
    calendar = 'standard'
    units = 'days since 1900-01-01 00:00:00'
    timestamps=[]
    for date in dates:
        timestamps.append(pd.Timestamp(date).to_datetime())

    return np.sort(date2num(timestamps, units, calendar))


def dt_to_dec(dt):
    # Datetime object to decimal year
    startyear = datetime(dt.year, 1, 1)
    endyear = startyear.replace(year=dt.year + 1)
    return dt.year + ((dt - startyear).total_seconds() / float((endyear - startyear).total_seconds()))


def save_as_mat(path, gpi, test_prod, ref_prod, anomaly, timeframe):

    # type: (int) -> None

    #Saves the SM data for active timeframe for testproduct and refproduct to the selected path
    # format: *gpi*.mat


    test_obj = HomogTest(test_prod,
                         ref_prod,
                         0.01,
                         anomaly)

    print('Exporting Testdata and Referencedata to .mat')
    exp_data = test_obj.read_gpi(gpi,
                                 test_obj.range[0],
                                 test_obj.range[1]) # type: pd.DataFrame
    exp_data = exp_data / 100

    matdate = []
    for dt in exp_data['refdata'].index:
        matdate.append(datetime2matlabdn(dt))

    x = {'tspan': matdate,
         'sm': exp_data['refdata'].values}

    y = {'tspan': matdate,
         'sm': exp_data['testdata'].values}

    timeframe = [datetime2matlabdn(timeframe[0]),
                 datetime2matlabdn(timeframe[1])]


    data_dict = {'X': x, 'Y': y, 'timeframe': timeframe}
    io.savemat(os.path.join(path, str(gpi)),
               data_dict, oned_as='column')

