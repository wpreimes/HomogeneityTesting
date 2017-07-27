# -*- coding: utf-8 -*-
"""
Created on Thu Jun 08 10:16:42 2017

@author: wpreimes
"""
import sys
sys.path.append(r"H:\workspace")

import pandas as pd
import numpy as np
from scipy import stats
import os
import pygeogrids.netcdf as nc
from datetime import datetime, timedelta
from pynetcf.time_series import GriddedNcIndexedRaggedTs
from pynetcf.point_data import GriddedPointData
from points_to_netcdf import points_to_netcdf
from sklearn import linear_model
from otherfunctions import dates_to_num, regress
from interface import HomogTest
from otherfunctions import cci_timeframes
from points_to_netcdf import pd_from_2Dnetcdf, create_cellfile_name
from grid_functions import cells_for_continent


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
            testresults = testresults.rename(columns={'test_results': breaktime})
            dataframes.append(testresults)

        return pd.concat(dataframes, axis=1)

    def _init_adjustment_times(self, test_results, breaktimes, timeframes):
        adjustment_times = pd.DataFrame(index=test_results.index.values,
                                        columns=[timeframes[0][0]] + np.ndarray.tolist(breaktimes) + [
                                            timeframes[-1][-1]])
        adjustment_times[timeframes[0][0]] = 99  # Timeframes are special

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
            raise Exception('1: negative correlation found, adjustment NOT performed')
        if any(p > 0.05 for p in pval):
            raise Exception('2: positive correlation not significant, adjustment NOT performed')

        # Perform the linear adjustment for testdata after the breaktime
        # cc = B[0][1] / B[1][1]
        # dd = B[0][0] - cc * B[1][0]
        # Perform the linear adjustment for testdata before the breaktime
        # TODO: Is this correct??
        cc = B[1][1] / B[0][1]
        dd = B[1][0] - cc * B[0][0]
        # adjusted_part1 = dataframe[['testdata']][:breaktime]
        # adjusted_part2 = cc * dataframe[['testdata']][breaktime:] + dd
        adjusted_part1 = cc * dataframe[['testdata']][:breaktime] + dd
        adjusted_part2 = dataframe[['testdata']][breaktime:]

        #TODO: Test adjustment

        dataframe['adjust'] = pd.concat([adjusted_part1, adjusted_part2])
        # dataframe['adjusted'] = dataframe.adjust[breaktime:]
        if return_whole:
            adjusted = dataframe.adjust
            self.regression_adjustment(dataframe, breaktime)
            # dataframe['adjusted'] = dataframe.adjust
        else:
            adjusted = dataframe.adjust[:breaktime]
            # dataframe['adjusted'] = dataframe.adjust[:breaktime]
        adj_setting = {'status': '0: adjustment performed',
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
        starttime = datetime.strptime(adjustment_pattern.index.values[0], '%Y-%m-%d')
        endtime = datetime.strptime(adjustment_pattern.index.values[-1], '%Y-%m-%d')
        breaktimes = []
        for breaktime_string in adjustment_pattern.loc[adjustment_pattern == 1].index.values:
            breaktimes.append(datetime.strptime(breaktime_string, '%Y-%m-%d'))
        breaktimes = breaktimes[::-1]
        data = self.read_gpi(gpi)[starttime:endtime]

        if len(breaktimes) == 0:
            return data[['testdata']].rename(columns={'testdata': 'not_adjusted'})

        df_adjusted = pd.DataFrame(columns=['testdata'])
        adj_settings = {}
        for i, breaktime in enumerate(breaktimes):
            if i == 0:
                # first iteration never uses new data
                if len(breaktimes) > 1:
                    data_to_adjust = data[breaktimes[i + 1]:]  # ----BT-----|BT-----BT-----end|
                else:
                    data_to_adjust = data  # |start-----BT-----end|
            elif i != len(breaktimes) - 1:  # common case
                if any(pd.notnull(df_adjusted[breaktime:])):
                    # Merge adjusted data (before break) and unadjusted data (after break)
                    testdata_merge = pd.concat([data['testdata'][breaktimes[i + 1]:breaktime],
                                                df_adjusted['testdata'][breaktime:]])
                    # Add refdata (bias corrected) for adjustement
                    data_to_adjust = pd.DataFrame(index=testdata_merge.index,
                                                  data={'testdata': testdata_merge,
                                                        'refdata': data['refdata'][breaktimes[i + 1]:]})
                    # -----BT-------|BT--------BT--adjusted data--|
                else:
                    data_to_adjust = data[breaktimes[i + 1]:]
            elif i == len(breaktimes) - 1:  # last iteration
                if any(pd.notnull(df_adjusted[breaktime:])):
                    testdata_merge = pd.concat([data['testdata'][:breaktime],
                                                df_adjusted['testdata'][breaktime:]])
                    data_to_adjust = pd.DataFrame(index=testdata_merge.index,
                                                  data={'testdata': testdata_merge,
                                                        'refdata': data['refdata']})
                    # -----BT-------|BT--------BT--adjusted data--|
                else:
                    data_to_adjust = data  # no adjusted data available and last breaktime, use whole frame as is

            try:
                bias_corr_refdata, rxy, pval, ress = regress(data_to_adjust)
                data_to_adjust['bias_corr_refdata'] = bias_corr_refdata
                data_to_adjust = data_to_adjust[['testdata', 'bias_corr_refdata']].rename(
                    columns={'bias_corr_refdata': 'refdata'})

                if i == 0:
                    return_whole = True
                else:
                    return_whole = False
                settings, adjusted_data = self.regression_adjustment(data_to_adjust, breaktime, return_whole)

            except Exception as e:
                settings = {'status': e[0], 'slope': np.nan,
                            'intercept': np.nan, 'model': np.nan}
                # TODO: SAve this also in the netcdf file
                # TODO: Adjust daily data!!!
                if i == len(breaktimes) - 1:
                    adjusted_data = data_to_adjust['testdata'][:breaktime]
                else:
                    adjusted_data = data_to_adjust['testdata'][breaktimes[i + 1]:breaktime]
            # ts_adjusted = df_adjusted['testdata'].append(adjusted_data)
            if df_adjusted['testdata'].index.size == 0:
                df_adjusted['testdata'] = adjusted_data
            else:
                save = (df_adjusted['testdata'].append(adjusted_data)).sort_index()
                df_adjusted = pd.DataFrame(columns=['testdata'])
                df_adjusted['testdata'] = save

            adj_settings[breaktime.strftime('%Y-%m-%d')] = settings

        data['adjusted'] = df_adjusted['testdata']

        return adj_settings, data[['adjusted']]



def run_adjustment():

    adjust_obj = Adjust(r"H:\HomogeneityTesting_data\output\CCI31EGU",
                        'adjusted_cci',
                        'merra2',
                        0.1)


    grid_path = r"D:\users\wpreimes\datasets\grids\qdeg_land_grid.nc"
    cell_grid = nc.load_grid(grid_path)
    cells = cells_for_continent('Australia')

    adjusted_data_path = r"D:\users\wpreimes\datasets\CCI_31_D\adjusted_temp"
    # TODO: Add own list for points should have been adjusted but could not
    unadjusted_gps = {'gpi': [], 'slope': [], 'intercept': [], 'adjustment_class': []}  # 0= UNadjusted, 1=adjusted

    dataset = GriddedNcIndexedRaggedTs(path=adjusted_data_path, grid=cell_grid, mode='w')
    adjusted_gps = GriddedPointData(os.path.join(adjusted_data_path,'adjusted_gps.nc'), cell_grid, mode='w')

    for cell_index, cell in enumerate(cells):
        gpis = cell_grid.grid_points_for_cell(cell)[0]
        adjustment_status={'gpi':[], 'status':[]}
        adjustment_stats={}
        for index, gpi in enumerate(gpis):
            if index % 50 == 0:
                print('%i of %i' % (index, gpis.size))
            try:
                adj_settings, adjusted_data = adjust_obj.adjust_ts(gpi)
            except:
                adjustment_status['gpi'].append(gpi)
                adjustment_status['status'].append(0)
                continue
            for breaktime_str, settings in adj_settings.iteritems():
                if breaktime_str not in adjustment_stats.keys():
                    adjustment_stats[breaktime_str] = {'gpi' : [],
                                                       'intercept': [], 'slope': []}

                adjustment_stats[breaktime_str]['gpi'].append(gpi)
                adjustment_stats[breaktime_str]['intercept'].append(adj_settings[breaktime_str]['intercept'])
                adjustment_stats[breaktime_str]['slope'].append(adj_settings[breaktime_str]['slope'])

            if adjusted_data.columns.values[0] == 'not_adjusted':
                adjustment_status['gpi'].append(gpi)
                adjustment_status['status'].append(0)
            if adjusted_data.columns.values[0] == 'adjusted':
                adjustment_status['gpi'].append(gpi)
                adjustment_status['status'].append(1)


            dataset.write(gpi, adjusted_data)
        dataset.close()


        points_to_netcdf(pd.DataFrame(index=adjustment_status['gpi'],
                                      data={'status': adjustment_status['status']}),
                         path=adjusted_data_path,
                         filename='adjustment_status')
        for breaktime_str in adjustment_stats.keys():
            points_to_netcdf(pd.DataFrame(index=adjustment_stats[breaktime_str]['gpi'],
                                          data={'intercept': adjustment_stats[breaktime_str]['intercept'],
                                                'slope': adjustment_stats[breaktime_str]['slope']}),
                             path=adjusted_data_path,
                             filename=breaktime_str+'_adj_stats')


        '''
        adjusted_gps = {}
        adjusted_ts_data = pd.DataFrame()
        if os.path.isfile(os.path.join(adjusted_data_path, str(cell) + '.nc')): continue  # already processed
        print('cell %i of %i' % (cell_index, len(cells)))
        gpis = cell_grid.grid_points_for_cell(cell)[0]
        for index, gpi in enumerate(gpis):
            if index % 50 == 0:
                print('%i of %i' % (index, gpis.size))
            try:
                adj_settings, adjusted_data = adjust_obj.adjust_ts(gpi)
            except:
                continue
            column_name = adjusted_data.columns.values[0]
            if column_name == 'not_adjusted':
                unadjusted_gps['gpi'].append(gpi)
            if column_name == 'adjusted':
                for breaktime_str, settings in adj_settings.iteritems():
                    if breaktime_str not in adjusted_gps.keys():
                        adjusted_gps[breaktime_str] = pd.DataFrame(index=cell_grid.grid_points_for_cell(cell)[0],
                                                                   data={'intercept':np.nan, 'slope':np.nan})

                    adjusted_gps[breaktime_str].loc[gpi,['intercept','slope']]=[settings['slope'], settings['intercept']]

            if adjusted_ts_data.index.size == 0:
                adjusted_ts_data = adjusted_data.rename(columns={column_name: gpi})
            else:
                adjusted_ts_data[gpi] = adjusted_data

            dataset.write(gpi,adjusted_data)

        for breaktime_str, data in adjusted_gps.iteritems():
            points_to_netcdf(dataframe=data, path=adjusted_data_path,
                             filename='adjustment_stats_' + breaktime_str)

        '''
        '''
        dates = np.array([date.to_datetime() for date in pd.to_datetime(adjusted_ts_data.index.values)])
        data = {'sm_adjusted': np.array([adjusted_ts_data[gpi].values for gpi in adjusted_ts_data.columns.values])}
        loc_ids = adjusted_ts_data.columns.values
        write_gpis_to_file(data, dates, loc_ids, adjusted_data_path, str(cell) + '.nc')
        '''


#run_adjustment()
