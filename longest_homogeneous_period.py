# -*- coding: utf-8 -*-
"""
Created on Jul 13 12:22 2017

@author: wpreimes
"""
from scipy import stats
import numpy as np
import pandas as pd
import os
import glob
from datetime import datetime, timedelta
from points_to_netcdf import pd_from_2Dnetcdf, points_to_netcdf
from save_data import load_Log
import re
from datetime import datetime
from HomogeneityTesting.otherfunctions import cci_timeframes, dt_to_dec
from typing import Union

def calc_longest_homogeneous_period(workdir, create_netcdf=True):
    # type: (str) -> Union(pd.DataFrame,str)
    '''
    Calculate the duration of the longest homogeneouse period of the data in the working folder
    Needs testing over all break times of a certain product
    :param workdir: Path to netcdf data from homogeneity testing
    :return: Dataframe of all homogeneous Periods and their lengths and the start year / end year
            of the longest homogeneitouse period
    '''
    fileslist = glob.glob(os.path.join(workdir, "HomogeneityTest*.nc"))
    filenames = [afile.replace(workdir + '\\', '') for afile in fileslist]

    log = load_Log(workdir)
    products = log.get_products()
    ref_prod = products['ref_prod']
    test_prod = products['test_prod']

    mergetimes = cci_timeframes(test_prod)
    breaktimes = mergetimes['breaktimes']

    dates = [datetime.strptime(ts, "%Y-%m-%d") for ts in breaktimes]
    sort_order = np.argsort(np.asarray(dates))

    breaktimes = np.asarray(breaktimes)[sort_order]
    starttimes = np.array([timeframe[0] for timeframe in mergetimes['timeframes']])[sort_order]
    endtimes = np.array([timeframe[1] for timeframe in mergetimes['timeframes']])[sort_order]

    DF_Period = pd.DataFrame()
    for breaktime in breaktimes:
        breaktime_re = re.compile(".+%s.nc" % breaktime)
        for filename in filenames:
            if re.match(breaktime_re, filename):
                test_results = pd_from_2Dnetcdf(os.path.join(workdir, filename), grid='land')
                #TODO: h_all can be removed
                if 'test_results' not in test_results.columns.values:
                    DF_Period = pd.concat([DF_Period, test_results[['h_all']]],axis=1)
                    DF_Period = DF_Period.rename(columns={'h_all': datetime.strptime(breaktime, '%Y-%m-%d')})
                else:
                    DF_Period = pd.concat([DF_Period, test_results[['test_results']]],axis=1)
                    DF_Period = DF_Period.rename(columns={'test_results':  datetime.strptime(breaktime, '%Y-%m-%d')})

    starttimes = [datetime.strptime(time, '%Y-%m-%d') for time in starttimes]
    breaktimes = [datetime.strptime(time, '%Y-%m-%d') for time in breaktimes]
    endtimes = [datetime.strptime(time, '%Y-%m-%d') for time in endtimes]

    periodsizes = []
    i = 1
    for starttime, breaktime, endtime in zip(starttimes, breaktimes, endtimes):
        # DF_Period['%i.1'%i] = np.where(DF_Period[breaktime]==4.,(endtime-starttime).days, ((breaktime-starttime).days),(endtime-breaktime).days))

        DF_Period['%i' % i] = (breaktime - starttime).days
        periodsizes.append((breaktime - starttime).days)
        # DF_Period['%i.2'%i] = (endtime-breaktime).days
        i += 1
    DF_Period['%i' % i] = (endtimes[-1] - breaktimes[-1]).days


    periods = map(list, DF_Period[DF_Period.columns.values[-i:]].values)
    #breaktimes_str = [breaktime.strftime('%Y-%m-%d') for breaktime in breaktimes]
    calcs = map(list, DF_Period[breaktimes].values)
    results = []

    for values, ops in zip(periods, calcs):
        ops.append([4.0, np.nan])
        r = []
        while ops:
            if ops[0] != 4.0:
                r.append(values[0])
                del values[0]
                del ops[0]

            elif ops[0] == 4.0:
                values[1] += values[0]
                del ops[0]
                del values[0]
        # r.append(values[-1])
        results.append(r)


    max_period = []
    startdays = []
    i=0
    j=0
    for periods in results:
        #if periods[:-1] == periodsizes:
        if max(periods) in periodsizes:
            i+=1
            max_period.append(np.nan)
            startdays.append(np.nan)
        else:
            j+=1
            pos = periods.index(max(periods))
            max_period.append(periods[pos])
            startdays.append(sum(periods[:pos]))

    first_date = starttimes[0]
    startdates = []
    startyears = []
    enddates = []
    endyears = []

    for start, maxperiod in zip(startdays, max_period):
        if np.isnan(start) or np.isnan(maxperiod):
            startdates.append(np.nan)
            enddates.append(np.nan)
            startyears.append(np.nan)
            endyears.append(np.nan)
        else:
            startdate = first_date + timedelta(days=int(start))
            enddate = startdate + timedelta(days=int(maxperiod))
            startdates.append(startdate)
            enddates.append(enddate)
            startyears.append(dt_to_dec(startdate))
            endyears.append(dt_to_dec(enddate))

    startdates = np.asarray(startdates)
    enddates = np.asarray(enddates)
    max_period = np.asarray(max_period) / 365.
    DF_Period['max_Period'] = max_period
    DF_Period['startdate'] = startdates
    DF_Period['endate'] = enddates
    DF_Period['startyear'] = startyears
    DF_Period['endyear'] = endyears

    if create_netcdf:
        ncfilename = 'HomogeneousPeriod'
        var_meta_dicts = {'max_Period': {'description': 'Shows Length of the longest homogeneous Period of the pixel, no test = break'},
                          'startyear': {'description': 'Shows the year (decimal number) when the longest homogeneous period started'},
                          'endyear': {'description': 'Shows the year (decimal number) when the longest homogeneous period ended'},
                          }
        points_to_netcdf(DF_Period[['max_Period','startyear','endyear']],
                         workdir, None, ncfilename,
                         None, var_meta_dicts)

        return ncfilename+'.nc'
    else:
        return DF_Period

#DF_Period = calc_longest_homogeneous_period(r"H:\HomogeneityTesting_data\output\CCI33", create_netcdf=False)
#print(DF_Period.loc[242363])