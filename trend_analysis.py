# -*- coding: utf-8 -*-
"""
Created on Sep 16 11:33 2017

@author: wpreimes
"""

from scipy.stats import theilslopes
from interface import BreakTestData
import matplotlib.pyplot as plt
import numpy as np
from pygeogrids.netcdf import load_grid
from smecv_grid.grid import SMECV_Grid_v042
import pandas as pd
from datetime import datetime
from grid_functions import cells_for_continent
import os
import xarray as xr
from mktest import mk_test
from scipy import stats
from multiprocessing import Pool
from otherfunctions import split
from multiprocessing import Process, Queue
import matplotlib.dates as mdates
from save_data import RegularGriddedCellData
from grid_functions import split_cells
from cci_timeframes import CCITimes
from datetime import datetime

from scipy.stats import linregress
from nc_image_reading import RegularGriddedCellData_Reader

class PlotTrends(object):
    def __init__(self, trends_file, homogtest_file, testdata_name, refdata_name, only_significant_trends, outlier_filter_trends=None):
        self.trends_file = RegularGriddedCellData_Reader(trends_file)
        if homogtest_file:
            self.homogtest_file = RegularGriddedCellData_Reader(homogtest_file)
        else:
            self.homogtest_file = None
        self.outlier_filter_trends = outlier_filter_trends
        self.testdata_name = testdata_name
        self.refdata_name = refdata_name
        self.only_significant = only_significant_trends

    @staticmethod
    def remove_outliers(column_in, b, t, as_array=False):
        '''
        Removes all data from a dataframe column above and below the selected
        percentiles

        Parameters
        ---------
        column: pandas dataframe column to filter

        b: bottom percentile in % (0-100)
        t: top percentile in % (0-100)
        Returns
        ---------
        column: The filtered data column

        points: a list of points which where removed during filtering
        '''
        column = column_in.copy(True).dropna()
        toppercentile = np.nanpercentile(column, t)
        bottompercentile = np.nanpercentile(column, b)
        points = column[(column > toppercentile) | (column < bottompercentile)].index.values
        size = points.size
        column[(column > toppercentile) | (column < bottompercentile)] = np.nan
        print 'Removed %i values using percentiles(%i and %i): %f and %f' % (
        size, b, t, bottompercentile, toppercentile)
        return column


    def load_trends(self, time):
        if self.only_significant:
            slope_var = 'slope_lin_significant'
        else:
            slope_var = 'slope_lin'
        DF_Points_from_file = self.trends_file.read(time, ['%s_%s' % (self.testdata_name, slope_var),
                                                     '%s_%s' % (self.refdata_name, slope_var)]).dropna()
        if self.outlier_filter_trends:
            b, t = self.outlier_filter_trends[0], self.outlier_filter_trends[1]
            for col in DF_Points_from_file:
                DF_Points_from_file[col] = self.remove_outliers(DF_Points_from_file[col], b, t)

        return DF_Points_from_file.rename(columns={'%s_%s' % (self.testdata_name, slope_var): 'slope_%s' % self.testdata_name,
                                                   '%s_%s' % (self.refdata_name, slope_var): 'slope_%s' % self.refdata_name})


    def load_breaks(self, time):
        DF_Points_from_file = self.homogtest_file.read(time, ['test_results']).dropna()
        return DF_Points_from_file


    @staticmethod
    def calc_scatter_stats(Err0, Err1):

        N = Err0.size if Err0.size == Err1.size else np.nan
        RMSD = np.sqrt(np.nanmean(Err1 - Err0) ** 2)
        Bias = np.nanmean(Err1) - np.nanmean(Err0)

        mask1 = ~np.isnan(Err1) & ~np.isnan(Err0)
        R = stats.pearsonr(Err0[mask1], Err1[mask1])
        slope = linregress(Err0[mask1], Err1[mask1])
        return N, RMSD, Bias, slope, R

    @staticmethod
    def boxposition(box_pos):
        p = []
        if box_pos[0] == 't':
            p.append(0.85)
        elif box_pos[0] == 'b':
            p.append(0.15)
        else:
            return 'Position is not tl,tr,bl or br'

        if box_pos[1] == 'r':
            p.append(0.85)
        elif box_pos[1] == 'l':
            p.append(0.25)
        else:
            return 'Position is not tl,tr,bl or br'

        return p

    def scatter_trends(self, starttime, out_path=None,
                       regressionline=True, oneoneline=True, at='breaks', box_pos='tl', scale='log'):
        '''
        :param starttime: datetime
            Time to read image from netcdf file
        :param out_path: str
            Path to save the plot to
        :param regressionline: bool
            Add a regression line to the scatter plot
        :param oneoneline: bool
            Add a 1:1 line to the scatter plot
        :param at: str
            'breaks' to make scatter plot for trends where breaks were found
        :param box_pos: str
            Position of stats-box in the scatter plot
        :return:
        '''

        trends = self.load_trends(starttime)
        print starttime
        if at == 'breaks':
            if not self.homogtest_file: raise Exception('No break test file passed')
            test_results = self.load_breaks(starttime).dropna()
            index_pos = test_results.loc[test_results['test_results']!=4.,:].index
        elif at == 'nobreaks':
            if not self.homogtest_file: raise Exception('No break test file passed')
            test_results = self.load_breaks(starttime).dropna()
            index_pos = test_results.loc[test_results['test_results']==4.,:].index
        elif at == 'all':
            test_results = None
            index_pos = trends.index
        else:
            raise Exception("Select 'breaks' or 'all' for parameter 'at' ")

        trends = trends.loc[index_pos]
        x = trends['slope_%s' % self.testdata_name]
        y = trends['slope_%s' % self.refdata_name]

        p = self.boxposition(box_pos)
        fig = plt.figure()
        ax = fig.add_subplot(111)

        N, RMSD, Bias, slope, R = self.calc_scatter_stats(trends['slope_' + self.testdata_name].values,
                                                          trends['slope_'+ self.refdata_name].values)

        if scale == 'log':
            hb = ax.hexbin(x, y,
                           gridsize=100, bins='log', cmap='Reds', mincnt=1)
        else:
            hb = ax.hexbin(x, y,
                           gridsize=100, cmap='Reds', mincnt=1)


        '''
        max_x = np.nanmax(np.abs(x))
        max_y = np.nanmax(np.abs(y))
        max_x = 0.006
        max_y = 0.006
        ax.set_xlim([max_x * -1, max_x])
        ax.set_ylim([max_y * -1, max_y]) #TODO: change to +-0.1?? #TODO: change to +-0.1??
        #max = np.max(max_x, max_y)
        '''
        # ax.axis([np.nanmin(x), np.nanmax(x), np.nanmin(y), np.nanmax(y)])


        title = "%s Trends vs. %s Trends" % (self.testdata_name, self.refdata_name)
        ax.set_title(title)
        ax.set_xlabel(self.testdata_name)
        ax.set_ylabel(self.refdata_name)
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label('log10(N)')

        plt.axvline(0, color='black', alpha=0.3)
        plt.axhline(0, color='black', alpha=0.3)

        ax.text(p[1], p[0], 'RMSD:%f\n' % RMSD + \
                'N:%i\n' % N + \
                'R:%f' % R[0], horizontalalignment='center', backgroundcolor='white', verticalalignment='center',
                transform=ax.transAxes, fontsize=15)

        if regressionline == True:
            ax.plot(x, slope[0] * x + slope[1], linestyle='-', color='blue', linewidth=2)

        plt.tight_layout()
        if out_path:
            fig.savefig(os.path.join(out_path, title + '(from '+str(starttime.date())+').png'))
            plt.close()

        else:
            plt.show()

        return fig

class Trends(BreakTestData):
    def __init__(self, test_prod, ref_prod, anomaly, use_adjusted=False):
        BreakTestData.__init__(self, test_prod, ref_prod, anomaly)

    @staticmethod
    def TheilSenTrend(df_in, alpha=0.95, plotit=False):
        '''
        implements a method for robust linear regression.
        It computes the slope as the median of all slopes between paired values.
        :param df_in:
        :param alpha:
        :param plotit:
        :return:
        '''
        df = df_in.copy()
        return_df = pd.DataFrame(index=['medslope', 'medintercept', 'lo_slope', 'up_slope'],
                                 columns=df.columns.values)

        if plotit:
            n_col = df.columns.values.size
            ymin = np.min(df.min().values)
            ymax = np.max(df.max().values)
            fig, axs = plt.subplots(n_col, sharex=True)

        df['t'] = range(df.index.values.size)

        for i, col_name in enumerate(df.columns.values[df.columns.values != 't']):
            y = df[[col_name, 't']].dropna()[col_name].values
            t = df[[col_name, 't']].dropna()['t'].values
            medslope, medintercept, lo_slope, up_slope = theilslopes(y, t, alpha)

            return_df = return_df.set_value('medslope', col_name, medslope)
            return_df = return_df.set_value('medintercept', col_name, medintercept)
            return_df = return_df.set_value('lo_slope', col_name, lo_slope)
            return_df = return_df.set_value('up_slope', col_name, up_slope)

            if plotit:
                index = df[[col_name]].dropna().index.values
                # x = np.array(range(y.size))
                axs[i].plot(index, y, 'b.')
                axs[i].plot(index, medintercept + medslope * t, 'r-')
                axs[i].set_ylim([ymin, ymax])
                axs[i].set_title(col_name)

        if plotit:
            plt.tight_layout()
            plt.show()

        return return_df


    @staticmethod
    def MannKendall(df_in, alpha=0.95, plotit=False):

        df = df_in.copy()
        return_df = pd.DataFrame(index=['trend', 'slope', 'intercept', 'p'], columns=df.columns.values)

        if plotit:
            n_col = df.columns.values.size
            ymin = np.min(df.min().values)
            ymax = np.max(df.max().values)
            fig, axs = plt.subplots(n_col, sharex=True)

        df['t'] = range(df.index.values.size)

        for i, col_name in enumerate(df.columns.values[df.columns.values != 't']):
            y = df[[col_name, 't']].dropna()[col_name].values
            t = df[[col_name, 't']].dropna()['t'].values

            mk, m, c, p = mk_test(t=t, x=y, eps=0, alpha=alpha, Ha='upordown')
            return_df = return_df.set_value('trend', col_name, mk)
            return_df = return_df.set_value('slope', col_name, m)
            return_df = return_df.set_value('intercept', col_name, c)
            return_df = return_df.set_value('p', col_name, p)

            if plotit:
                index = df[[col_name]].dropna().index.values
                # x = np.array(range(y.size))
                axs[i].plot(index, y, 'b.')
                axs[i].plot(index, c + m * t, 'r-')
                axs[i].set_ylim([ymin, ymax])
                axs[i].set_title(col_name)

        if plotit and any(return_df.loc['trend'].values)==True:
            print return_df
            plt.show()
        else:
            plt.close()

        return return_df

    def adjusted_gpis(self, breaktest_grid_path):
        grid = load_grid(breaktest_grid_path)
        return grid.get_grid_points()[0]



def process_for_cells(q, grid, cells, Trends_obj, save_obj, starttimes, endtimes, testdata_name, refdata_name, plotit, process_no):
    for icell, cell in enumerate(cells):
        print 'Processing QDEG Cell %i (iteration %i of %i)' % (cell, icell + 1, len(cells))
        for iteration, gpi in enumerate(grid.grid_points_for_cell(cell)[0]):
            print gpi
            try:
                DF_Time = Trends_obj.read_gpi(gpi).rename(columns={'testdata_original':'testdata'})
            except:
                continue
            for starttime, endtime in zip(starttimes, endtimes):
                try:
                    df_time = DF_Time.loc[starttime:endtime]
                    df_time = Trends_obj.temp_resample(df_time, '3M')
                    # TODO: deactivate correction
                    df_time['refdata_original'] = Trends_obj.ref_data_correction(df_time[['testdata', 'refdata_original']],
                                                                        'refdata_original', False)
                    df_time = df_time.rename(columns={'testdata': testdata_name, 'refdata_original': refdata_name})

                    mka_results = Trends_obj.MannKendall(df_time, 0.05, plotit)
                    #theilsen_results = Trends_obj.TheilSenTrend(DF_Time, 0.05, plotit)

                    # if not mannkendall_params.loc['p', testdata_name]<p_min: continue
                    # theilsen_params = Trends_obj.TheilSenTrend(DF_Time.dropna() if plotit else DF_Time, p_min, plotit)

                # print theilsen_params
                # print mannkendall_params

                    data_dict = {}
                    for name in [testdata_name, refdata_name]:
                        data_dict['%s_significant_trend_mask' % name] = 1 if mka_results.loc['trend', name] == True else 0
                        #data_dict['%s_slope_theilsen' % name] = theilsen_results.loc['medslope', name]
                        data_dict['%s_slope_lin' % name] = mka_results.loc['slope', name]
                        data_dict['%s_p_mannkendall' % name] = mka_results.loc['p', name]
                        data_dict['%s_slope_lin_significant' % name] = mka_results.loc['slope', name] if mka_results.loc[
                                                                       'trend', name] == True else np.nan

                    save_obj.add_data(data_dict, gpi, time=datetime.strptime(starttime,'%Y-%m-%d'))
                except:
                    continue

    q.put(True)

def calc_trends(out_dir, testdata_name, refdata_name, starttimes, endtimes, cells_identifier, plotit, parallel_processes):

    if parallel_processes != 1 and plotit:
        plotit=False

    grid = SMECV_Grid_v042()

    save_obj = RegularGriddedCellData(os.path.join(out_dir, 'temp_cellfiles', testdata_name),
                                      grid,
                                      times=[datetime.strptime(starttime,'%Y-%m-%d') for starttime in starttimes])

    Trends_obj = Trends(testdata_name, refdata_name, False, False)


    cells = list(
        split(split_cells(cells_identifier, SMECV_Grid_v042()), parallel_processes))  # split cells for processes

    processes = []
    q = Queue()
    finished_processes = []

    for process_no in range(parallel_processes):
        cells_for_process = cells[process_no]
        p = Process(target=process_for_cells, args=(q,
                                                    grid,
                                                    cells_for_process,
                                                    Trends_obj,
                                                    save_obj,
                                                    starttimes,
                                                    endtimes,
                                                    testdata_name,
                                                    refdata_name,
                                                    plotit,
                                                    process_no))
        processes.append(p)
        p.start()

    for process in processes:
        finished_processes.append(q.get(True))

    for process in processes:
        process.join()

    global_file_name = save_obj.make_global_file(out_dir,
                                                 'GLOBAL_%s_Trends.nc' % testdata_name,
                                                 False, False,
                                                 keep_cell_files=True)
    return

def start_trends():
    out_dir = r'D:\users\wpreimes\datasets\CCITrends_output\cci_adjusted'
    starttime = datetime(1988, 1, 1)
    endtime = datetime(2010, 12, 31)

    times =CCITimes('CCI_41_COMBINED').get_times()
    timeframes = times['timeframes']
    timeframes = np.insert(timeframes, -2, [['1988-01-01','2016-12-31']], axis=0)
    refdata_names = ['merra2']
    testdata_names = ['CCI_41_COMBINED_ADJUSTED']

    for refdata_name in refdata_names:
        for testdata_name in testdata_names:
            calc_trends(out_dir=out_dir,
                        testdata_name=testdata_name,
                        refdata_name=refdata_name,
                        starttimes=timeframes[:,0],
                        endtimes=timeframes[:,1],
                        cells_identifier='Australia',
                        plotit=False,
                        parallel_processes=4)



if __name__ == '__main__':
    trends_file = r"D:\users\wpreimes\datasets\CCITrends_output\cci_unadjusted\GLOBAL_CCI_41_COMBINED_Trends.nc"
    homogtest_file = r"D:\users\wpreimes\datasets\HomogeneityTesting_data\AUS_tf_fit\GLOBAL_test_results_before_adjustment.nc"
    trends_obj = PlotTrends(trends_file,
                            homogtest_file,
                            'CCI_41_COMBINED',
                            'merra2',
                            only_significant_trends=False,
                            outlier_filter_trends=(1,99))


    # start time for full trend: time = datetime(1988,1,1)
    times =CCITimes('CCI_41_COMBINED').get_times(as_datetime=True)
    starttimes = times['timeframes'][:,0]
    #times = np.insert(times, -2, [datetime(1988,1,1)], axis=0)
    for starttime in starttimes:
        trends_obj.scatter_trends(starttime, r'D:\users\wpreimes\datasets\CCITrends_output\cci_unadjusted',
                           regressionline=True, oneoneline=True, at='nobreaks', box_pos='tl', scale='linear')