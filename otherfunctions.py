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


def join_files(filefolder, filelist):
    merged_file_name = 'saved_points.csv'
    merged = []

    for file in filelist:
        filepath = os.path.join(filefolder, file)
        data = csv_read_write(filepath, 'read')
        for row in data:
            merged.append(row)
        os.remove(filepath)

    merged_int = map(int, [item for sublist in merged for item in sublist])
    path = csv_read_write(os.path.join(filefolder, 'saved_points.csv'), 'write', merged_int)

    return path, np.asarray(merged_int)


def csv_read_write(csv_path, mode, data=None):
    if mode == 'write':
        if not os.path.isfile(csv_path):
            with open(csv_path, 'wb') as file:
                wr = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
                wr.writerow(data)
        else:
            if os.path.isfile(csv_path):
                with open(csv_path, 'ab') as file:
                    wr = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
                    wr.writerow(data)
        return csv_path
    if mode == 'read':
        return_data = []
        with open(csv_path, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                return_data.append(row)
        return return_data


def regress(data):
    '''
    Perform regression of column refdata on column testdata
    '''
    dataframe = data.copy()
    out = dataframe.refdata
    dataframe = dataframe.dropna()
    R, pval = stats.pearsonr(dataframe.refdata, dataframe.testdata) # Correlation between Refdata and Testdata


    if R < 0 or np.isnan(R):
        ress = [np.nan]
        return out, R, pval, ress

    testdata = dataframe.testdata.values
    refdata = dataframe.refdata.values
    refdata_ones = np.vstack([refdata, np.ones(len(refdata))]).T

    ress = np.linalg.lstsq(refdata_ones, testdata)[0][::-1]
    dataframe['ones'] = 1
    xm = np.matrix(dataframe.as_matrix(columns=['ones', 'refdata']))
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

def fill_holes_lress(fill_series, other_series):
    '''
    Perform regression of column refdata on column testdata as used for combining
    multiple insitu measurements
    '''
    # Find holes in max_series where there is data in other_series
    df = pd.DataFrame(data={'fill_series': fill_series,
                            'other_series': other_series},
                      index=fill_series.index)
    # Drop values where both sets are nan, no interpolation possible here
    df = df.dropna(how='all')
    # Group blocks where contiguous values are Nan in the series to fill

    df['nanmask'] = df['fill_series'].isnull().astype(bool)
    df['nanmask_shift'] = df['nanmask'].shift(1)
    df['nan_change'] = df['nanmask'] != df['nanmask_shift']
    df['nangroup'] = df['nan_change'].cumsum()

    holes = df[['fill_series', 'other_series', 'nangroup', 'nanmask']].groupby(['nangroup', 'nanmask'])
    # calculate regression coefficients for holes from other_series
    for count, hole in holes:
        if hole.nanmask.all():
            slope, intercept, rvalue, pvalue, stderr = stats.linregress(range(len(hole.fill_series)),
                                                                        hole.other_series)

            df.ix[hole.index, 'fill_series'] = intercept + slope * hole.other_series

    return df['fill_series']

    # fill holes in max_series from regression coefficients found


def merge_ts(dataframe_in):
    '''
    Merge temporally coinciding timeseries via interpolation from
    bivariate linear regression.
    Input: Dataframe of temporally coinciding timeseries.
    -Find TS with most values. Use a reference to enhance.
    -Enhance missing values in reference via linear regression from other TS.
        Only use correlating (R>0.8) TS for interpolation.
    -Return gap filled reference TS
    '''
    dataframe = dataframe_in
    max_series_name = np.argmax(dataframe.notnull().sum().to_frame()[0])
    quick_corr = dataframe.corr().sort_values(max_series_name, ascending=False).index.values[1:]
    # Order other TSs via correlation to max_series and return list corr_ordered
    # Improve max_series first with best fitting other,than second best fitting other,...
    DF_other_series = pd.DataFrame()
    for other_series_name in dataframe[quick_corr]:
        corr, pval = stats.pearsonr(dataframe[[max_series_name, other_series_name]].dropna()[max_series_name],
                                    dataframe[[max_series_name, other_series_name]].dropna()[other_series_name])
        if corr >= 0.8 and pval < 0.01:
            DF_other_series[other_series_name] = dataframe[other_series_name]


    max_series = dataframe[max_series_name]
    other_series_merged = DF_other_series.mean(axis=1)
    if max_series.isnull().all():
        return None
    elif other_series_merged.isnull().all():
        # If there is no second TS, return the unchanged reference series
        return max_series
    else:
        # Else perform interpolation from linear regression
        merged_ts = fill_holes_lress(max_series, other_series_merged)
        return merged_ts


