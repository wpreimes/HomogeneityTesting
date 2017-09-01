# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 17:06:28 2017

@author: wpreimes
"""
import sys
from rsdata import root_path
import os

import numpy as np

import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from cci_timeframes import CCITimes
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime
from longest_homogeneous_period import calc_longest_homogeneous_period
from save_data import load_Log
from smecv_grid.grid import SMECV_Grid_v042
import xarray as xr
from interface import HomogTest
from otherfunctions import resample_to_monthly

import platform

def valid_months_plot(workdir, testproduct, refproduct, min_monthly_values):

    #NOT IN PROCESS

    # Count the number of valid months i.e, month that contain >=
    # the minimum amount of values

    #Count values BEFORE breaktime

    #Count value AFTER breaktime
    time_obj = CCITimes(testproduct)

    grid = SMECV_Grid_v042() # type: CellGrid()


    #grid = grid.subgrid_from_cells([564]) #TODO: Delete

    gpis = grid.get_grid_points()[0]
    lons = grid.get_grid_points()[1]
    lats = grid.get_grid_points()[2]


    test_obj = HomogTest(testproduct, refproduct,['wilkoxon', 'fligner_killeen'], 0.01, False,None)
    DF_Points = pd.DataFrame(index = gpis, data = {'lon': lons, 'lat': lats} )
    for breaktime in time_obj.get_times(ignore_conditions=True)['breaktimes']:
        DF_Points['before %s' % breaktime] = np.nan
        DF_Points['after %s' % breaktime] = np.nan

    for iteration, gpi in enumerate(gpis):
        times = time_obj.get_times(ignore_conditions=True, as_datetime=False)
        starts = np.flipud(np.array([timeframe[0] for timeframe in times['timeframes']]))
        ends = np.flipud(np.array([timeframe[1] for timeframe in times['timeframes']]))
        breaktimes = np.flipud(times['breaktimes'])

        if iteration %100 == 0:
            print 'processing gpi %i (%i of %i)' % (gpi, iteration, gpis.size)
        try:
            df_time = test_obj.read_gpi(gpi, starts[0], ends[-1])

            for start, breaktime, end in zip(starts, breaktimes, ends):
                start = datetime.strptime(start, '%Y-%m-%d')
                breaktime = datetime.strptime(breaktime, '%Y-%m-%d')
                end = datetime.strptime(end, '%Y-%m-%d')

                df_subtime = df_time[['testdata']][start:end].dropna()
                df_resample = resample_to_monthly(df_subtime, min_monthly_values)
                df_group, len_bef, len_aft = \
                    test_obj.group_by_breaktime(df_resample, breaktime, 3, ignore_exception=True)

                DF_Points.loc[gpi, 'before %s' % str(breaktime.date())] = len_bef
                DF_Points.loc[gpi, 'after %s' % str(breaktime.date())] = len_aft
        except:
            continue
    DF_Points = DF_Points.sort_values(['lat', 'lon']) \
                         .set_index(['lat', 'lon'])

    global_image = DF_Points.to_xarray()
    global_image.to_netcdf(os.path.join(workdir, '%s_valid_monthly_obs.nc' % testproduct))



def inhomo_plot_with_stats(workdir, filename):
    # type: (str) -> dict
    '''
    :param workdir: path to directory containing nc files from HomogeneityTesting
    :return: None
    '''
    ncpath = os.path.join(workdir, filename)
    ccigrid = SMECV_Grid_v042() # type: CellGrid()
    land_points = ccigrid.get_grid_points()[0]

    log = load_Log(workdir)
    products = log.get_products()
    ref_prod = products['ref_prod']
    test_prod = products['test_prod']

    ncfile = xr.open_dataset(ncpath)
    DF_Points = pd.DataFrame(index=range(0, 720 * 1440))
    DF_Points_from_file = ncfile.to_dataframe().reset_index(inplace=False).set_index(['location_id'])
    DF_Points_from_file = DF_Points_from_file.loc[DF_Points_from_file.index.dropna()]
    for name in DF_Points_from_file.columns.values:
        DF_Points[name] = DF_Points_from_file[name]

    splitname = filename.split('_')

    breaktime_str = splitname[1]

    all_gpis = land_points.size
    tested_gpis = DF_Points['test_results'].loc[DF_Points['test_results'].isin([1., 2., 3., 4.])].size
    hwk_gpis = DF_Points['test_results'].loc[DF_Points['test_results'].isin([1., 3.])].size
    hfk_gpis = DF_Points['test_results'].loc[DF_Points['test_results'].isin([2., 3.])].size
    hboth_gpis = DF_Points['test_results'].loc[DF_Points['test_results'].isin([3.])].size

    all_tested = (float(tested_gpis) / float(all_gpis)) * 100.
    try:
        wk_tested = (float(hwk_gpis) / float(tested_gpis)) * 100.
        fk_tested = (float(hfk_gpis) / float(tested_gpis)) * 100.
        both_tested = (float(hboth_gpis) / float(tested_gpis)) * 100.
    except:
        wk_tested, fk_tested, both_tested = 0, 0, 0

    textbox = 'Ground Points Tested: %.2f%%\n' % all_tested + \
              'WK (of Tested): %.2f%%\n' % wk_tested + \
              'FK (of Tested): %.2f%%\n' % fk_tested + \
              'Both (of Tested): %.2f%%' % both_tested

    colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (.5, .5, .5)]
    cmap = LinearSegmentedColormap.from_list('HomogeneityGroups',
                                             colors,
                                             N=4)

    img = np.empty(DF_Points.index.size, dtype='float32')
    img.fill(np.nan)
    img[DF_Points.index] = DF_Points['test_results'].values

    # mask array where invalid values (nans) occur
    img_masked = np.ma.masked_invalid(img.reshape((180 * 4, 360 * 4)))

    f = plt.figure(num=None, figsize=(20, 10), dpi=90, facecolor='w', edgecolor='k')

    m = Basemap(projection='mill', llcrnrlat=-60, urcrnrlat=80,
                llcrnrlon=-180, urcrnrlon=180, resolution='c', )

    m.drawcoastlines()
    m.drawcountries()
    lons = (np.arange(360 * 4) * 0.25) - 179.875
    lats = (np.arange(180 * 4) * 0.25) - 89.875
    lons, lats = np.meshgrid(lons, lats)

    im = m.pcolormesh(lons, lats, img_masked, cmap=cmap, latlon=True)

    im.set_clim(vmin=1, vmax=4)

    cb = m.colorbar(im, "bottom", size="5%", pad="8%")
    cb.outline.set_visible(False)
    cb.set_ticks([])

    for t in cb.ax.get_xticklabels():
        t.set_fontsize(20)

    cb.set_label('Break detected by Test:', fontsize=20, labelpad=-50)
    for j, lab in enumerate(['Wilkoxon', 'Fligner-Killeen', 'Both', 'None']):
        cb.ax.text((2 * j + 1) / 8.0, .5, lab, fontsize=15, ha='center', va='center')

    title = 'HomogeneityTesting \n %s|%s|Breaktime:%s' % (test_prod, ref_prod, breaktime_str)
    plt.title(title, fontsize=20)

    plt.annotate(textbox, fontsize=15, xy=(0.025, 0.05), xycoords='axes fraction',
                 bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 3})

    plotfile_name = 'HomogeneityTest_%s_%s' % (ref_prod, breaktime_str)
    plt.savefig(workdir + '\\' + plotfile_name + '.png', dpi=f.dpi)

    return {'tested_gps': all_tested, 'wk_of_tested' : wk_tested, 'fk_of_tested': fk_tested,
            'both_of_tested': both_tested}

'''
def compare_RTM_RTG(workdir, model_prod):

    times = {'ISMN-merge': {'starttimes': [], 'endtimes': [], 'breaktimes': []},
             model_prod: {'starttimes': [], 'endtimes': [], 'breaktimes': []}}
    files = {'ISMN-merge': [], model_prod: []}

    for ref_data in ['ISMN-merge', model_prod]:
        fileslist = glob.glob(os.path.join(workdir, "DF_Points_%s*.csv" % ref_data))
        files[ref_data] = fileslist
        for filename in fileslist:
            # TODO: Trollsift Parser
            filename = filename.replace(workdir, '')
            filename = filename.replace('.', '_')
            splitname = filename.split('_')
            times[ref_data]['starttimes'].append(splitname[3])
            times[ref_data]['endtimes'].append(splitname[5])
            times[ref_data]['breaktimes'].append(splitname[4])
     
    # if (times['ISMN_merge']['starttimes']==times[model_prod]['starttimes']) and \
    #    (times['ISMN_merge']['endtimes']==times[model_prod]['endtimes']) and \
    #    (times['ISMN_merge']['breaktimes']==times[model_prod]['breaktimes']):
    

    cols = int(math.ceil(float(len(files['ISMN-merge']))))
    fig, axs = plt.subplots(1, cols, figsize=(cols * 3, 5), facecolor='w', edgecolor='k', sharey=False)

    fig.suptitle('Comparison Break Detection RTM,RTG (USA)', fontsize=20, fontweight='bold')
    fig.subplots_adjust(top=0.85, wspace=0.2, hspace=0.1)
    axs = axs.ravel()
    comp_meas = []
    for i, (ismn_file, model_file) in enumerate(zip(files['ISMN-merge'], files[model_prod])):
        index = 0
        if times['ISMN-merge']['breaktimes'][i] == times[model_prod]['breaktimes'][i]:
            DF_Points_ismn = pd.read_csv(os.path.join(ismn_file), index_col=0)
            DF_Points_model = pd.read_csv(os.path.join(model_file), index_col=0)

            if (DF_Points_ismn.index == DF_Points_model.index).all():
                DF_Points_merged = DF_Points_model[['lat', 'lon', 'test_results']].rename(
                    columns={'test_results': model_prod})
                DF_Points_merged['ISMN-merge'] = DF_Points_ismn['test_results']
                DF_Points_merged = DF_Points_merged.dropna(how='any')
                DF_Points_merged['diff'] = np.nan
                if not DF_Points_merged.empty:
                    for gpi in DF_Points_merged.index.values:
                        if DF_Points_merged[model_prod].loc[gpi] == np.nan or DF_Points_merged['ISMN-merge'].loc[
                            gpi] == np.nan:
                            continue
                        elif DF_Points_merged[model_prod].loc[gpi] == DF_Points_merged['ISMN-merge'].loc[gpi]:
                            DF_Points_merged['diff'].loc[gpi] = 0
                        elif DF_Points_merged[model_prod].loc[gpi] != 4. and DF_Points_merged['ISMN-merge'].loc[
                            gpi] == 4.:
                            DF_Points_merged['diff'].loc[gpi] = 1
                        elif DF_Points_merged[model_prod].loc[gpi] == 4. and DF_Points_merged['ISMN-merge'].loc[
                            gpi] != 4.:
                            DF_Points_merged['diff'].loc[gpi] = -1
                        index += 1

                    data = DF_Points_merged['diff'].values
                    n, bins, patches = axs[i].hist(data, bins=3, range=(-1.5, 1.5), align='mid')
                    cm = plt.cm.get_cmap('hsv')
                    bin_centers = 0.5 * (bins[:-1] + bins[1:])
                    col = bin_centers - min(bin_centers)
                    col /= max(col)
                    for c, p in zip(col, patches):
                        plt.setp(p, 'facecolor', cm(c))
        comp_meas.append(index)
        axs[i].set_xlim([-2, 2])
        axs[i].set_ylim([0, 200])
        axs[i].set_xticks(np.arange(-2, 2, 1.0))
        axs[i].set_title('Break: ' + str(times[model_prod]['breaktimes'][i]))
    
    #fig.text(0.5, 0.03, 'Comparison Class', ha='center', va='center',fontsize=12)
    #fig.text(0.09, 0.5, 'Number of Tests', ha='center', va='center', rotation='vertical',fontsize=12)
    #fig.text(0.58, -0.04, '0=No difference between RTM and RTG \n 1=RTM found Break, RTG not \n -1=RTG found Break, RTM not', 
    #         style='italic',bbox={'facecolor':'grey', 'alpha':0.5, 'pad':10})     
    
    fig.savefig(os.path.join(workdir, 'RTM_vs_RTG'))
    print('Number of points to compare:' + str(comp_meas))
'''


def show_processed_gpis(workdir,plottype, filename):
    '''
    Calculate spatial plots for the areas where Homogeneity Tests were (not)
    performed
    '''
    ncpath = os.path.join(workdir, filename)
    ccigrid = SMECV_Grid_v042() # type: CellGrid()
    land_points = ccigrid.get_grid_points()[0]

    log = load_Log(workdir)
    products = log.get_products()
    ref_prod = products['ref_prod']
    test_prod = products['test_prod']
    ncfile = xr.open_dataset(ncpath)

    DF_Points = pd.DataFrame(index=range(0, 720 * 1440))
    DF_Points_from_file = ncfile.to_dataframe().reset_index(inplace=False).set_index(['location_id'])
    DF_Points_from_file = DF_Points_from_file.loc[DF_Points_from_file.index.dropna()]
    for name in DF_Points_from_file.columns.values:
        DF_Points[name] = DF_Points_from_file[name]


    splitname = filename.split('_')
    breaktime_str = splitname[1]


    if plottype=='test':
        statname = 'test_status'
        title = 'HomogeneityTesting Coverage \n %s|%s|%s' % (test_prod, ref_prod, breaktime_str)
        plotfile_name = 'RT_coverage_%s_%s' % (ref_prod, breaktime_str)
    else:
        statname = 'adj_status'
        title = 'BreakAdjustment Coverage \n %s|%s|%s' % (test_prod, ref_prod, breaktime_str)
        plotfile_name = 'adjustment_coverage_%s_%s' % (ref_prod, breaktime_str)

    status_var = ncfile.variables[statname].attrs['Values']




    colors = ["#2ecc71", "#9b59b6", "#3498db", "#95a5a6", "#e74c3c",
              "#34495e", "#FFC200", "#FFC0CB", "#FF1493", "#FFFF00"]

    N = len(colors)
    cmap = LinearSegmentedColormap.from_list(statname,
                                             colors,
                                             N=N)

    img = np.empty(DF_Points.index.size, dtype='float32')
    img.fill(np.nan)
    img[DF_Points.index] = DF_Points[statname].values

    # mask array where invalid values (nans) occur
    img_masked = np.ma.masked_invalid(img.reshape((180 * 4, 360 * 4)))

    f = plt.figure(num=None, figsize=(20, 10), dpi=90, facecolor='w', edgecolor='k')

    m = Basemap(projection='mill', llcrnrlat=-60, urcrnrlat=80,
                llcrnrlon=-180, urcrnrlon=180, resolution='c', )

    m.drawcoastlines()
    m.drawcountries()


    lons = (np.arange(360 * 4) * 0.25) - 179.875
    lats = (np.arange(180 * 4) * 0.25) - 89.875
    lons, lats = np.meshgrid(lons, lats)

    im = m.pcolormesh(lons, lats, img_masked, cmap=cmap, latlon=True)

    im.set_clim(vmin=0, vmax=N)

    cb = m.colorbar(im, "bottom", size="5%", pad="8%")
    cb.outline.set_visible(True)
    cb.set_ticks([])

    for t in cb.ax.get_xticklabels():
        t.set_fontsize(20)

    cb.set_label('Error occurred during:', fontsize=20, labelpad=-50)

    for j in range(N):
        cb.ax.text((float(j)/float(N)) + float(1./float(N*2)), .5, str(j), fontsize=15, ha='center', va='center')

    plt.title(title, fontsize=20)

    stats = DF_Points[statname].dropna().values
    all_gpis = float(land_points.size)

    # Grab the groups and their values from the meta data of the nc file
    meta ={}

    groups_strings = status_var.split(',')
    for string in groups_strings:
        string = string.split('=')
        meta.update({int(string[0]): string[1]})

    textbox = 'Total ground points: %i\n' % all_gpis

    for val, s in meta.iteritems():
        stat = float(stats[np.where(stats == float(val))].size / all_gpis) * 100.
        textbox += str(val)+': ' + s + ': ' + '%.2f%%\n' % stat

    plt.annotate(textbox, fontsize=10, xy=(0.025, 0.05), xycoords='axes fraction',
                 bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 3})


    plt.savefig(workdir + '\\' + plotfile_name + '.png', dpi=f.dpi)

    return meta

def longest_homog_period_plots(workdir, startyear_plot=True, endyear_plot=True):
    
    #Create Plot of the longest homogeneouse Period of all nc files in the working dir
    #:param workdir: str
    #    Path to the directory containing the nc files returned by Homogeneity Test
    #:param same_colorbar: bool
    #    True to use fixed colobar (for comparsion of results for differenct cci versions
    #    False to use a different colorbar for different cci versions
    #:return: None

    filename = calc_longest_homogeneous_period(workdir)
    ncfile = xr.open_dataset(os.path.join(workdir, filename))

    DF_Points = pd.DataFrame(index=range(0, 720 * 1440))
    DF_Points_from_file = ncfile.to_dataframe().reset_index(inplace=False).set_index(['location_id'])
    DF_Points_from_file = DF_Points_from_file.loc[DF_Points_from_file.index.dropna()]
    for name in DF_Points_from_file.columns.values:
        DF_Points[name] = DF_Points_from_file[name]


    log = load_Log(workdir)
    products = log.get_products()
    ref_prod = products['ref_prod']
    test_prod = products['test_prod']

    lons = (np.arange(360 * 4) * 0.25) - 179.875
    lats = (np.arange(180 * 4) * 0.25) - 89.875
    lons, lats = np.meshgrid(lons, lats)

    cmap = plt.cm.jet

    levels = [8, 11.5, 15.5, 21]
    cmaplist = [(0.5, 0, 0, 1), (0.7, 0, 0, 1), (1, 0.2, 0.2, 1), (1, 0.2, 0.2, 1),
                (0, 0.7, 0.7, 1), (0, 0.6, 0.7, 1), (0, 0.6, 0.7, 1), (0, 0.3, 1, 1),
                (0, 0, 0.5, 1), (0, 0, 0.5, 1), (0, 0, 0.5, 1)]

    cmap = cmap.from_list('LongestPeriodCMap', cmaplist, cmap.N)

    img = np.empty(DF_Points.index.size, dtype='float32')
    img.fill(np.nan)
    img[DF_Points.index] = DF_Points['max_Period'].values

    # mask array where invalid values (nans) occur
    img_masked = np.ma.masked_invalid(img.reshape((180 * 4, 360 * 4)))

    f = plt.figure(num=None, figsize=(20, 10), dpi=90, facecolor='w', edgecolor='k')

    m = Basemap(projection='mill', llcrnrlat=-60, urcrnrlat=80, \
                llcrnrlon=-180, urcrnrlon=180, resolution='c', )

    m.drawcoastlines()
    m.drawcountries()

    im = m.pcolormesh(lons, lats, img_masked, cmap=cmap, latlon=True)

    im.set_clim(vmin=levels[0], vmax=levels[-1])

    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin=levels[0], vmax=levels[-1]))
    sm._A = []
    cb = m.colorbar(sm, "bottom", size="5%", pad="8%")

    cb.set_ticklabels(["%.1f" % lev for lev in levels])
    for t in cb.ax.get_xticklabels():
        t.set_fontsize(20)

    cb.set_label('Period Length [years]:', fontsize=20, labelpad=-75)

    title = 'Length of Longest Homogeneous Period\n %s|%s' % (test_prod, ref_prod)
    plt.title(title, fontsize=20)

    filename = 'LongestHomogPeriod.png'
    plt.savefig(workdir + '\\' + filename, dpi=f.dpi)

    return filename
'''
def extreme_break(workdir, ref_prod, test_prod):
    gpis = [363662, 406520, 564947, 391490]
    for gpi in gpis:
        if gpi == 564947:
            ttime = ['1998-01-01', '2002-07-01', '2007-01-1']
        else:
            ttime = ['1992-01-01', '1998-01-01', '2002-01-01']
        test_data = QDEGdata_M(products=[test_prod, ref_prod])
        df = test_data.read_gpi(gpi, ttime[0], ttime[2])
        df = df.rename(columns={ref_prod: 'refdata', test_prod: 'testdata'})
        corr, pval = stats.spearmanr(df['testdata'], df['refdata'])
        df['bias_corr_refdata'], rxy, pval, ress = regress(df)

        fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
        df = df.dropna()
        df = df.rename(columns={'refdata': 'Model Reference Data', 'testdata': 'CCI SM data'})
        df[['CCI SM data']].plot(ax=ax1, color='blue')
        df[['Model Reference Data']].plot(ax=ax2, color='green')
        plt.setp(ax2.get_yticklabels()[-1], visible=False)

        fig.suptitle('GPI %i: Example Break' % gpi, fontsize=20, fontweight='bold')
        fig.text(0.5, 0.05, 'Date', ha='center', va='center', fontsize=12)
        fig.text(0.05, 0.5, 'Soil moisture [m3 m-3]', ha='center', va='center', rotation='vertical', fontsize=12)
        fig.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

        filename = 'ExtremeBreak_%s_%i' % (ref_prod, gpi)
        plt.savefig(workdir + '\\' + filename + '.png')
'''
if __name__ == '__main__':
    pass
    # df = pd_from_2Dnetcdf(r"H:\HomogeneityTesting_data\output\CCI33\HomogeneityTest_merra2_1991-08-01.nc", return_only_tested=False)
    # extreme_break(r'H:\HomogeneityTesting_data\output\CCI22EGU','gldas_v2','cci_22')
    #show_tested_gpis(r"H:\HomogeneityTesting_data\output\v5","HomogeneityTestResult_2007-01-01_image.nc")
    #longest_homog_period_plots(r"H:\HomogeneityTesting_data\output\v5")
    # show_processed_gpis(r"H:\HomogeneityTesting_data\output\CCI33", 'test',"HomogeneityTest_merra2_1991-08-01.nc" )
    #inhomo_plot_with_stats(r'H:\HomogeneityTesting_data\output\v5',"HomogeneityTestResult_2007-01-01_image.nc")
    # compare_RTM_RTG(r'H:\HomogeneityTesting_data\output\CCI31EGU','merra2')
    # show_tested_gpis(r"H:\HomogeneityTesting_data\output\v15","HomogeneityTest_merra2_2011-10-01")
    valid_months_plot(r'H:\HomogeneityTesting_data\output\CCI_available_data_plots', 'CCI_33_COMBINED','merra2', 10)