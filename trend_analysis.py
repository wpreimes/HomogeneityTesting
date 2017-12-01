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
from preprocess import remove_outliers
from multiprocessing import Pool
from otherfunctions import split
from multiprocessing import Process, Queue
import matplotlib.dates as mdates


def calc_scatter_stats(Err0, Err1):
    from scipy.stats import linregress
    RMSD = np.sqrt(np.nanmean(Err1 - Err0) ** 2)
    Bias = np.nanmean(Err1) - np.nanmean(Err0)

    mask1 = ~np.isnan(Err1) & ~np.isnan(Err0)
    R = stats.pearsonr(Err0[mask1], Err1[mask1])
    slope = linregress(Err0[mask1], Err1[mask1])
    return RMSD, Bias, slope, R


def scatter_trends(ncfiles, testdata_name, refdata_name,slope_name, save_path=None,
                   regressionline=True, oneoneline=True, box_pos='tl'):
    '''
    Plot scatter plots of trends from saved nc files
    :param nc_files:
    :return:
    '''
    p = []
    if box_pos[0] == 't':
        p.append(0.9)
    elif box_pos[0] == 'b':
        p.append(0.1)
    else:
        return 'Position is not tl,tr,bl or br'

    if box_pos[1] == 'r':
        p.append(0.9)
    elif box_pos[1] == 'l':
        p.append(0.2)
    else:
        return 'Position is not tl,tr,bl or br'
    concat_df = []
    for ncfile in ncfiles:
        ncfile = xr.open_dataset(ncfile)
        concat_df.append(ncfile.to_dataframe().reset_index(inplace=False).set_index(['location_id']))
    DF_Points_from_file = pd.concat(concat_df)
    DF_Points_from_file = DF_Points_from_file.loc[DF_Points_from_file.index.dropna()]

    DF_Points_from_file = DF_Points_from_file[['%s_slope_%s' % (testdata_name, slope_name),
                                               '%s_slope_%s' % (refdata_name, slope_name)]].dropna()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = remove_outliers(DF_Points_from_file['%s_slope_%s' % (testdata_name, slope_name)],
                        5, 95)[0].values
    y = remove_outliers(DF_Points_from_file['%s_slope_%s' % (refdata_name, slope_name)],
                        5, 95)[0].values
    RMSD, Bias, slope, R = calc_scatter_stats(x, y)
    hb = ax.hexbin(x, y,
                   gridsize=100, bins='log', cmap='Reds', mincnt=1)
    # ax.axis([np.nanmin(x), np.nanmax(x), np.nanmin(y), np.nanmax(y)])


    title = "%s Trends %s vs. refdata" % (slope_name, testdata_name)
    ax.set_title(title)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('log10(N)')

    plt.axvline(0, color='black', alpha=0.3)
    plt.axhline(0, color='black', alpha=0.3)

    ax.text(p[1], p[0], 'RMSD:%f\n' % RMSD + \
            'Bias:%f\n' % Bias + \
            'R:%f' % R[0], horizontalalignment='center', backgroundcolor='white', verticalalignment='center',
            transform=ax.transAxes)

    if regressionline == True:
        ax.plot(x, slope[0] * x + slope[1], linestyle='-', color='blue', linewidth=2)

    if save_path:
        fig.savefig(os.path.join(save_path, title + '.png'))
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

        return return_df

    def adjusted_gpis(self, breaktest_grid_path):
        grid = load_grid(breaktest_grid_path)
        return grid.get_grid_points()[0]


def calc_trends(out_dir, testdata_name, refdata_name, starttime, endtime, cells='global', gpis=False):
    plotit = True
    grid = SMECV_Grid_v042()


    if cells == 'global':
        cells = grid.get_cells()
    elif any(isinstance(cell,str) for cell in cells):
        cells = cells_for_continent(cells)
    grid = grid.subgrid_from_cells(cells)

    Trends_obj = Trends(testdata_name, refdata_name, False, False)
    concat_df = []

    if cells:
        grid = grid.subgrid_from_cells(cells)
    if gpis:
        grid = grid.subgrid_from_gpis(gpis)

    for i, cell in enumerate(grid.get_cells(), start=1):
        DF_Cell = pd.DataFrame(index=grid.grid_points_for_cell(cell)[0], columns=['lat', 'lon'])
        DF_Cell['lat'] = grid.grid_points_for_cell(cell)[2]
        DF_Cell['lon'] = grid.grid_points_for_cell(cell)[1]

        print('iteration %i of %i' % (i, len(cells)))

        for gpi in grid.grid_points_for_cell(cell)[0]:
            print gpi
            try:
                DF_Time = Trends_obj.read_gpi(gpi)
                DF_Time = DF_Time[starttime:endtime]
                DF_Time = Trends_obj.temp_resample(DF_Time, '3M')
                #TODO: deactivate correction
                DF_Time['refdata'] = Trends_obj.ref_data_correction(DF_Time[['testdata', 'refdata']])
                DF_Time = DF_Time.rename(columns={'testdata': testdata_name, 'refdata': refdata_name})
            except:
                continue
            try:
                mka_results = Trends_obj.MannKendall(DF_Time, 0.05, plotit)
                #theilsen_results = Trends_obj.TheilSenTrend(DF_Time, 0.05, plotit)

                # if not mannkendall_params.loc['p', testdata_name]<p_min: continue
                # theilsen_params = Trends_obj.TheilSenTrend(DF_Time.dropna() if plotit else DF_Time, p_min, plotit)
            except:
                continue

            # print theilsen_params
            # print mannkendall_params

            for name in [testdata_name, refdata_name]:
                DF_Points = DF_Points.set_value(gpi, '%s_significant_trend_mask' % name,
                                                1 if mka_results.loc['trend', name ] == True else 0)
                DF_Points = DF_Points.set_value(gpi, '%s_slope_theilsen' % name,
                                                theilsen_results.loc['medslope', name])
                DF_Points = DF_Points.set_value(gpi, '%s_slope_lin' % name,
                                                mka_results.loc['slope', name])
                DF_Points = DF_Points.set_value(gpi, '%s_p_mannkendall' % name,
                                                mka_results.loc['p', name])
                DF_Points = DF_Points.set_value(gpi, '%s_slope_lin_significant' % name,
                                                mka_results.loc['slope', name] if mka_results.loc['trend', name ] == True else np.nan)


        concat_df.append(DF_Points)

    DF_Points = pd.concat(concat_df)
    DF_Points['location_id'] = DF_Points.index.astype(int)

    DF_Points = DF_Points.sort_values(['lat', 'lon']) \
        .set_index(['lat', 'lon'])
    global_image = DF_Points.to_xarray()
    global_image.to_netcdf(os.path.join(out_dir, 'trends_%s_%s.nc' % (refdata_name, testdata_name)))



def main(multicore=True):
    cells = ['Australia']
    out_dir = r'D:\users\wpreimes\datasets\CCITrends_output'
    parallel_processes=1

    starttime = datetime(1980, 1, 1)
    endtime = datetime(2016, 12, 31)

    slope_names = ['mannkendall']
    refdata_names = ['merra2']
    testdata_names = ['CCI_41_COMBINED', 'CCI_41_COMBINED_ADJUSTED']
    files_path = r"D:\users\wpreimes\datasets\CCITrends_output"

    for refdata_name in refdata_names:
        for testdata_name in testdata_names:
            calc_trends(out_dir, testdata_name, refdata_name, starttime, endtime, cells='United_States', gpis=None)







if __name__ == "__main__":

    main()


    '''
    for refdata_name in refdata_names:
        for testdata_name in testdata_names:
            for slope_name in slope_names:
                path = os.path.join(files_path,"trends_merra2_%s.nc" %testdata_name)

                fig1 = scatter_trends(path, testdata_name, refdata_name, slope_name, save_path=files_path)
    '''
