# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 12:58:31 2016

@author: wpreimes
"""

from scipy import stats
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

from datetime import datetime


def regress(dataframe):
    '''
    Perform regression of column refdata on column testdata
    '''
    # Berechnet Korrelation zwischen Testdaten und Referenzdaten
    R, pval = stats.pearsonr(dataframe.refdata, dataframe.testdata)
    out = dataframe.refdata
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

    return out, R, pval, ress


def datetime2matlabdn(dt):
    ord = dt.toordinal()
    mdn = dt + timedelta(days=366)
    frac = (dt - datetime(dt.year, dt.month, dt.day, 0, 0, 0)).seconds / (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac


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


def merge_ts(dataframe):
    '''
    Merge temporally coinciding timeseries via interpolation from
    bivariate linear regression.
    Input: Dataframe of temporally coinciding timeseries.
    -Find TS with most values. Use a reference to enhance.
    -Enhance missing values in reference via linear regression from other TS.
        Only use correlating (R>0.8) TS for interpolation.
    -Return gap filled reference TS
    '''
    max_series_name = np.argmax(dataframe.notnull().sum().to_frame()[0])
    quick_corr = dataframe.corr().sort(max_series_name, ascending=False).index.values[1:]
    # Order other TSs via correlation to max_series and return list corr_ordered
    # Improve max_series first with best fitting other,than second best fitting other,...
    DF_other_series = pd.DataFrame()
    for other_series_name in dataframe[quick_corr]:
        corr, pval = stats.pearsonr(dataframe[[max_series_name, other_series_name]].dropna()[max_series_name],
                                    dataframe[[max_series_name, other_series_name]].dropna()[other_series_name])
        if corr >= 0.8 and pval < 0.01:
            DF_other_series[other_series_name] = dataframe[other_series_name]

            # dataframe=dataframe.rename(columns={max_series_name:'refdata',other_series_name:'testdata'})

            # out,R,pval,ress=regress(dataframe)
            # DF_Sensor[max_series_name+'_new']=out

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

def dt_to_dec(dt):
    # Datetime object to decimal year
    startyear = datetime(dt.year, 1, 1)
    endyear = startyear.replace(year=dt.year+1)
    return dt.year + ((dt - startyear).total_seconds() / float((endyear - startyear).total_seconds()))


def cci_timeframes(product, skip_times=None):
    # TODO: Add timeframes for active and passive products for 22 and 33
    # TODO: Check timeframe and breaktimes for CCI33

    #Order matters!!!!
    timeframes = {'cci_22_combined': np.array(
                                      [['2011-10-01', '2014-12-31'], ['2007-01-01', '2012-07-01'],
                                       ['2002-07-01', '2011-10-01'], ['1998-01-01', '2007-01-01'],
                                       ['1991-08-01', '2002-07-01'], ['1987-07-01', '1998-01-01']]),
                  'cci_31_combined': np.array(
                                      [['2011-10-01', '2015-05-01'], ['2010-07-01', '2012-07-01'],
                                       ['2007-10-01', '2011-10-01'], ['2007-01-01', '2010-07-01'],
                                       ['2002-07-01', '2007-10-01'], ['1998-01-01', '2007-01-01'],
                                       ['1991-08-01', '2002-07-01'], ['1987-09-01', '1998-01-01']]),
                  'cci_31_passive': np.array(
                                     [['2012-07-01', '2015-12-31'], ['2011-10-01', '2015-05-01'],
                                      ['2010-07-01', '2012-07-01'], ['2007-10-01', '2011-10-01'],
                                      ['2002-07-01', '2010-07-01'], ['1998-01-01', '2007-10-01'],
                                      ['1987-09-01', '2002-07-01']]),
                  'cci_31_active': np.array(
                                     [['2007-01-01', '2015-12-01'], ['2003-02-01', '2012-11-01'],
                                      ['1997-05-01', '2007-01-01'], ['1991-08-01', '2003-02-01']]),
                  'cci_33_combined': np.array(
                                     [['2012-07-01', '2016-12-31'], ['2011-10-01', '2015-05-01'],
                                      ['2010-07-01', '2012-07-01'], ['2007-10-01', '2011-10-01'],
                                      ['2007-01-01', '2010-07-01'], ['2002-07-01', '2007-10-01'],
                                      ['1998-01-01', '2007-01-01'], ['1991-08-01', '2002-07-01'],
                                      ['1987-09-01', '1998-01-01']])
               }

    breaktimes = {'cci_22_combined': np.array(
                                      ['2012-07-01', '2011-10-01', '2007-01-01', '2002-07-01', '1998-01-01',
                                       '1991-08-01']),
                  'cci_31_combined': np.array(
                                      ['2012-07-01', '2011-10-01', '2010-07-01', '2007-10-01', '2007-01-01',
                                       '2002-07-01', '1998-01-01', '1991-08-01']),
                  'cci_31_passive': np.array(
                                      ['2015-05-01', '2012-07-01', '2011-10-01', '2010-07-01', '2007-10-01',
                                       '2002-07-01', '1998-01-01']),
                  'cci_31_active': np.array(['2012-11-01', '2007-01-01', '2003-02-01', '1997-05-01']),
                  'cci_33_combined': np.array(
                                      ['2015-05-01', '2012-07-01', '2011-10-01', '2010-07-01', '2007-10-01',
                                       '2007-01-01', '2002-07-01', '1998-01-01', '1991-08-01'])
                }

    ranges = {'cci_22_combined': ['1980-01-01', '2014-12-31'],
              'cci_31_combined': ['1980-01-01', '2015-12-31'],
              'cci_31_passive': ['1980-01-01', '2015-12-31'],
              'cci_31_active': ['1991-08-01', '2015-12-31'],
              'cci_32_combined': ['1980-01-01', '2015-12-31'],
              'cci_33_combined': ['1980-01-01', '2016-12-31']
              }

    if product not in timeframes.keys() or product not in breaktimes.keys() or product not in ranges.keys():
        raise Exception('No test times for selected cci product')
    return_timeframes = timeframes[product]
    return_breaktimes = breaktimes[product]
    return_ranges = ranges[product]

    if return_timeframes.shape[0] != return_breaktimes.size:
        raise Exception('Number of timeframes %i and breaktimes %i not equal' % (return_timeframes.size,
                                                                                 return_breaktimes.size))

    if skip_times:
        return_timeframes = np.array([return_timeframes[i] for i in range(return_timeframes.shape[0]) if i not in skip_times])
        return_breaktimes = np.array([return_breaktimes[i] for i in range(return_breaktimes.size) if i not in skip_times])


    return {'breaktimes': return_breaktimes, 'timeframes': return_timeframes, 'ranges': return_ranges}

