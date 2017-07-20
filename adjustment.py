# -*- coding: utf-8 -*-
"""
Created on Thu Jun 08 10:16:42 2017

@author: wpreimes
"""
import pandas as pd
import numpy as np
from scipy import stats
import os
import pygeogrids.netcdf as nc
import rsdata.root_path as root
from datetime import datetime, timedelta
from pynetcf.time_series import OrthoMultiTs

from HomogeneityTesting.otherfunctions import dates_to_num, regress
from HomogeneityTesting.interface import HomogTest
from HomogeneityTesting.otherfunctions import cci_timeframes
from HomogeneityTesting.points_to_netcdf import pd_from_2Dnetcdf, create_cellfile_name
from HomogeneityTesting.grid_functions import cells_for_continent

class Adjust(HomogTest):
    def __init__(self, workdir, test_prod, ref_prod, anomaly):
        HomogTest.__init__(self, test_prod, ref_prod, None, None, None, anomaly)
        self.workdir = workdir
        mergetimes = cci_timeframes(self.test_prod)
        self.breaktimes = mergetimes['breaktimes'][::-1]
        self.timeframes = mergetimes['timeframes'][::-1]
        test_results = self._init_import_test_results()
        self.adjustment_timeframes = self._init_adjustment_times(test_results,
                                                                 self.breaktimes,
                                                                 self.timeframes)



    def _init_import_test_results(self):
        dataframes = []
        for breaktime in self.breaktimes:
            filename = 'HomogeneityTest_%s_%s.nc' % (self.ref_prod, breaktime)
            testresults = pd_from_2Dnetcdf(os.path.join(self.workdir, filename),
                                           grid='land')[['test_results']]
            testresults= testresults.rename(columns={'test_results': breaktime})
            dataframes.append(testresults)

        return pd.concat(dataframes, axis=1)

    def _init_adjustment_times(self, test_results, breaktimes, timeframes):
        adjustment_times = pd.DataFrame(index=test_results.index.values,
                                             columns = [timeframes[0][0]] + np.ndarray.tolist(breaktimes) + [timeframes[-1][-1]])
        adjustment_times[timeframes[0][0]] = 99 #Timeframes are special

        for breaktime in breaktimes:
            adjustment_times[breaktime] = 1
            no_break_gpis = test_results[breaktime].loc[(test_results[breaktime] == 4) | \
                                                        (np.isnan(test_results[breaktime]))].index.values
            adjustment_times[breaktime].ix[no_break_gpis] = 0

        adjustment_times[timeframes[-1][-1]] = 99  # Timeframes are special
        return adjustment_times




    def regression_adjustment(self, dataframe, breaktime, return_whole=False):

        i1 = dataframe[:breaktime]
        i2 = dataframe[breaktime + pd.DateOffset(1):]

        B = []
        rval = []
        pval = []
        for data in [i1, i2]:
            refdata = np.stack((np.ones(data['refdata'].values.size), data['refdata'].values))
            testdata = data['testdata'].values
            b = np.dot(np.linalg.inv(np.asmatrix(np.dot(refdata, np.transpose(refdata)))),
                       np.matrix.transpose(np.asmatrix(np.dot(refdata, testdata))))
            B.append(b)
            r, p = stats.pearsonr(refdata[1], testdata)
            rval.append(r)
            pval.append(p)

        B = np.squeeze(np.asarray(B))

        regtest = {'R': rval, 'p': pval}
        # print('Regression coefficients are: ' + str(regtest))

        if any(r < 0 for r in rval):
            raise Exception('negative correlation found, adjustment NOT performed')
        if any(p > 0.05 for p in pval):
            raise Exception('positive correlation not significant, adjustment NOT performed')

        # Perform the linear adjustment for testdata after the breaktime
        #cc = B[0][1] / B[1][1]
        #dd = B[0][0] - cc * B[1][0]
        # Perform the linear adjustment for testdata before the breaktime
        #TODO: Is this correct??
        cc = B[1][1] / B[0][1]
        dd = B[1][0] - cc * B[0][0]

        #adjusted_part1 = dataframe[['testdata']][:breaktime]
        #adjusted_part2 = cc * dataframe[['testdata']][breaktime:] + dd
        # TODO: Is this correct??
        adjusted_part1 = cc * dataframe[['testdata']][:breaktime] + dd
        adjusted_part2 = dataframe[['testdata']][breaktime:]
        dataframe['adjust'] = pd.concat([adjusted_part1, adjusted_part2])
        #dataframe['adjusted'] = dataframe.adjust[breaktime:]
        if return_whole:
            adjusted = dataframe.adjust
            #dataframe['adjusted'] = dataframe.adjust
        else:
            adjusted = dataframe.adjust[:breaktime]
            #dataframe['adjusted'] = dataframe.adjust[:breaktime]
        adj_setting = {'status':'Adjustment performed',
                       'slope': cc, 'intercept': dd, 'model': B}
        '''
        print 'Adjustment was performed using the following parameters:'
        print 'Slope: %f' %adj_setting['slope']
        print 'Intercept: %f' %adj_setting['intercept']
        print 'Model: {0} and {1}'.format(adj_setting['model'][0],adj_setting['model'][1])
        '''
        del dataframe
        return adj_setting, adjusted


    def adjust_ts(self, gpi):
        # Perform adjustment of Timeseries AFTER breaktime (if inhomogeneity exists data AFTER the breaktime
        # is matched to data BEFORE breaktime)
        adjustment_pattern = self.adjustment_timeframes.loc[gpi]
        starttime = datetime.strptime(adjustment_pattern.index.values[0],'%Y-%m-%d')
        endtime = datetime.strptime(adjustment_pattern.index.values[-1],'%Y-%m-%d')
        breaktimes = []
        for breaktime_string in adjustment_pattern.loc[adjustment_pattern==1].index.values:
            breaktimes.append(datetime.strptime(breaktime_string, '%Y-%m-%d'))
        breaktimes = breaktimes[::-1]
        data = self.read_gpi(gpi)[starttime:endtime]

        if len(breaktimes) == 0:
            self.fill_ncfile(gpi, data['testdata'], self.adjusted_data_path, ts_column_name=self.test_prod)
            return


        df_adjusted=pd.DataFrame(columns=['testdata'])

        for i, breaktime in enumerate(breaktimes):
            if i == 0:
                # first iteration never uses new data
                if len(breaktimes) > 1:
                    data_to_adjust = data[breaktimes[i+1]:] # ----BT-----|BT-----BT-----end|
                else:
                    data_to_adjust = data # |start-----BT-----end|
            elif i != len(breaktimes)-1: # common case
                if any(pd.notnull(df_adjusted[breaktime:])):
                    # Merge adjusted data (before break) and unadjusted data (after break)
                    testdata_merge = pd.concat([data['testdata'][breaktimes[i+1]:breaktime],
                                                df_adjusted['testdata'][breaktime:]])
                    # Add refdata (bias corrected) for adjustement
                    data_to_adjust = pd.DataFrame(index = testdata_merge.index,
                                                  data ={'testdata': testdata_merge,
                                                         'refdata': data['refdata'][breaktimes[i+1]:]})
                    #-----BT-------|BT--------BT--adjusted data--|
                else:
                    data_to_adjust = data[breaktimes[i+1]:]
            elif i == len(breaktimes)-1: # last iteration
                if any(pd.notnull(df_adjusted[breaktime:])):
                    testdata_merge = pd.concat([data['testdata'][:breaktime],
                                                df_adjusted['testdata'][breaktime:]])
                    data_to_adjust = pd.DataFrame(index = testdata_merge.index,
                                                  data ={'testdata': testdata_merge,
                                                         'refdata': data['refdata']})
                    #-----BT-------|BT--------BT--adjusted data--|
                else:
                    data_to_adjust = data # no adjusted data available and last breaktime, use whole frame as is

            try:
                bias_corr_refdata, rxy, pval, ress = regress(data_to_adjust)
                data_to_adjust['bias_corr_refdata'] = bias_corr_refdata
                data_to_adjust = data_to_adjust[['testdata', 'bias_corr_refdata']].rename(columns={'bias_corr_refdata': 'refdata'})

                if i == 0:
                    return_whole = True
                else:
                    return_whole = False
                adj_settings, adjusted_data = self.regression_adjustment(data_to_adjust, breaktime, return_whole)

            except Exception as e:
                adj_setting = {'status': e, 'slope': None,
                               'intercept': None, 'model': None}
                #TODO: SAve this also in the netcdf file
                #TODO: Adjust daily data!!!
                if i == len(breaktimes)-1:
                    adjusted_data = data_to_adjust['testdata'][:breaktime]
                else:
                    adjusted_data = data_to_adjust['testdata'][breaktimes[i+1]:breaktime]
            #ts_adjusted = df_adjusted['testdata'].append(adjusted_data)
            if df_adjusted['testdata'].index.size == 0:
                df_adjusted['testdata']=adjusted_data
            else:
                save = (df_adjusted['testdata'].append(adjusted_data)).sort_index()
                df_adjusted = pd.DataFrame(columns=['testdata'])
                df_adjusted['testdata'] = save

        data['adjusted']=df_adjusted['testdata']
        data[['adjusted','testdata']].plot()
        return data[['adjusted']]

        # self.
        #self.fill_ncfile(gpi, data[['testdata']], path=self.adjusted_data_path, ts_column_name='testdata')



def fill_ncfile(gpi, ts, path, ts_column_name = None, grid = None):
    if not grid:
        grid = nc.load_grid(os.path.join(root.r, 'Datapool_processed', 'GLDAS', 'GLDAS_NOAH025_3H.020',
                                         'ancillary', 'GLDASv2_025_land_grid.nc'))

    cell, filename = create_cellfile_name(gpi, grid)
    grid_points = grid.grid_points_for_cell(cell)[0]
    filepath = os.path.join(path, filename + '.nc')
    lonlat = grid.gpi2lonlat(gpi)
    if ts_column_name:
        data = {ts_column_name: ts[ts_column_name].values}
    else:
        data = {'sm_adjusted': ts.values}

    if os.path.isfile(filepath):
        with OrthoMultiTs(filepath, mode='a') as ncfile:
            # TODO: dates_to_num is slow and dates are same for all gpis, remove from iterating part of function
            dates = dates_to_num(ts.index)
            ncfile.write_ts(loc_id=gpi, data=data, dates=dates,
                            loc_descr=None, lon=lonlat[0], lat=lonlat[1], alt=None,
                            fill_values=None, attributes=None, dates_direct=True)
    else:
        # Create file and add dates and gpi dimension
        # Dates are all Days between start of first time frame and end of last time frame

        with OrthoMultiTs(filepath, mode='w', n_loc=grid_points.size) as ncfile:
            ncfile.variables['location_id'][:] = grid_points
            dates = dates_to_num(ts.index).tolist()
            ncfile.extend_time(dates, direct = True)
        # Retry
        fill_ncfile(gpi, ts, path, ts_column_name, grid)


adjust_obj = Adjust(r"H:\HomogeneityTesting_data\output\CCI31EGU",
                    'cci_31_combined',
                    'merra2',
                    0.1)

# random test data
#startdate = datetime.strptime(adjust_obj.timeframes[0][0], '%Y-%m-%d')
#enddate = datetime.strptime(adjust_obj.timeframes[-1][-1], '%Y-%m-%d')
#dates = [startdate + timedelta(days=x) for x in range((enddate - startdate).days + 1)]
#ts = pd.DataFrame(index = dates, data = {'rand_data': np.random.rand(len(dates))})

gpis=[242363,242364]
adjusted_data_path = r"D:\users\wpreimes\datasets\CCI_adjusted"
for gpi in gpis:
    adjusted_data = adjust_obj.adjust_ts(gpi)
    fill_ncfile(gpi=gpi, ts=adjusted_data, path=adjusted_data_path, ts_column_name='adjusted')

'''
adjust_obj.fill_ncfile(gpi=242363,
                       ts=ts,
                       ts_column_name='rand_data',
                       path=r"D:\users\wpreimes\datasets\CCI_adjusted")
'''