# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 17:06:28 2017

@author: wpreimes
"""
import os
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from longest_homogeneous_period import calc_longest_homogeneous_period
from save_data import LogFile
from smecv_grid.grid import SMECV_Grid_v042
import xarray as xr
from cci_timeframes import CCITimes
from csv_writer import dict_csv_wrapper
from matplotlib.offsetbox import AnchoredText



def ts_adjustment_plot(data_in, breaktime, timeframe, products, lonlat, resample_m=True):
    pass
    fig, ax = plt.subplots()
    data_in.plot(ax=ax, title='CCI SM Adjustment')

    textbox = 'lat: %.2f, lon: %.2f%%\n' % (lonlat[1], lonlat[0]) +  \
              'mean unadjusted before bt: %.2f%%\n' % data_stats['mean_before_break_original'] + \
              'mean unadjusted after bt: %.2f%%\n' % data_stats['mean_after_break_original'] + \
              'mean adjusted before bt: %.2f%%\n' % data_stats['mean_before_break_testdata'] + \
              'mean adjusted after bt: %.2f%%\n' % data_stats['mean_after_break_testdata'] + \
              'mean refdata before bt: %.2f%%\n' % data_stats['mean_before_break_refdata'] + \
              'mean refdata after bt: %.2f%%\n' % data_stats['mean_after_break_refdata'] + \
              'var unadjusted: %.2f%%\n' % var_unadjusted + \
              'var adjusted: %.2f%%\n' % var_adjusted

    plt.annotate(textbox, fontsize=15, xy=(0.025, 0.05), xycoords='axes fraction',
                 bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 3})
    anchored_text = AnchoredText(info_text, loc=2)
    ax.plot(x, y)
    ax.add_artist(anchored_text)



def combined_testresults_plot(workdir, file_before_adjustment, file_after_adjustment):
    '''
    Combine test results before and after adjustment to a single, global file
    :param workdir:
    :param file1:
    :param file2:
    :return:
    '''
    ncpath_before = os.path.join(workdir, 'global_images', file_before_adjustment)
    ncpath_after = os.path.join(workdir, 'global_images', file_after_adjustment)
    ccigrid = SMECV_Grid_v042()  # type: CellGrid()
    land_points = ccigrid.get_grid_points()[0]

    log = LogFile(workdir)
    test_prod, ref_prod = log.get_products()

    ncfile_before = xr.open_dataset(ncpath_before)
    ncfile_after = xr.open_dataset(ncpath_after)

    for time in ncfile_before['time']:
        DF_Points = pd.DataFrame(index=range(0, 720 * 1440))

        DF_Points_from_file1 = ncfile_before.sel(time=time).to_dataframe().reset_index(inplace=False).set_index(['gpi'])
        DF_Points_from_file1 = DF_Points_from_file1.loc[DF_Points_from_file1.index.dropna()]
        DF_Points['test_results_before'] = DF_Points_from_file1['test_results']
        DF_Points_from_file2 = ncfile_after.sel(time=time).to_dataframe().reset_index(inplace=False).set_index(['gpi'])
        DF_Points_from_file2 = DF_Points_from_file2.loc[DF_Points_from_file2.index.dropna()]
        DF_Points['test_results_after'] = DF_Points_from_file2['test_results']
        DF_Points['test_results_merged'] = DF_Points['test_results_before'].dropna()
        DF_Points.loc[DF_Points['test_results_after'].dropna().index, 'test_results_merged'] = DF_Points['test_results_after'].dropna()

        breaktime_str = str(np.datetime64(np.datetime_as_string(time.values,
                                                                timezone='local')[:10]))

        all_gpis = land_points.size
        tested_gpis = DF_Points['test_results_merged'].loc[DF_Points['test_results_merged'].isin([1., 2., 3., 4.])].size
        hwk_gpis = DF_Points['test_results_merged'].loc[DF_Points['test_results_merged'].isin([1., 3.])].size
        hfk_gpis = DF_Points['test_results_merged'].loc[DF_Points['test_results_merged'].isin([2., 3.])].size
        hboth_gpis = DF_Points['test_results_merged'].loc[DF_Points['test_results_merged'].isin([3.])].size

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
        img[DF_Points.index] = DF_Points['test_results_merged'].values

        # mask array where invalid values (nans) occur
        img_masked = np.ma.masked_invalid(img.reshape((180 * 4, 360 * 4)))

        f = plt.figure(num=None, figsize=(20, 10), dpi=200, facecolor='w', edgecolor='k')

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

        plotfile_name = '%s_%s' % ('COMBINED_test_results', breaktime_str)
        plt.savefig(os.path.join(workdir, 'global_images', plotfile_name + '.png'), dpi=f.dpi)
        plt.close()

def inhomo_plot_with_stats(in_file, out_path, name_base=None):
    # type: (str) -> dict
    '''
    :param workdir: path to directory containing nc files from HomogeneityTesting
    :return: None
    '''
    ncpath = in_file
    if not name_base:
        path = os.path.normpath(in_file)
        filename = path.split(os.sep)[-1]
        name_base = 'STATUS_%s' % in_file

    ccigrid = SMECV_Grid_v042()  # type: CellGrid()
    land_points = ccigrid.get_grid_points()[0]

    ncfile = xr.open_dataset(ncpath)
    for time in ncfile['time']:
        DF_Points = pd.DataFrame(index=range(0, 720 * 1440))
        DF_Points_from_file = ncfile.sel(time=time).to_dataframe().reset_index(inplace=False).set_index(['gpi'])
        #DF_Points_from_file = ncfile.to_dataframe().reset_index(inplace=False).set_index(['location_id'])
        DF_Points_from_file = DF_Points_from_file.loc[DF_Points_from_file.index.dropna()]
        for name in DF_Points_from_file.columns.values:
            DF_Points[name] = DF_Points_from_file[name]

        breaktime_str = str(np.datetime64(np.datetime_as_string(time.values,
                                                                timezone='local')[:10]))

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

        f = plt.figure(num=None, figsize=(20, 10), dpi=200, facecolor='w', edgecolor='k')

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

        title = 'HomogeneityTesting \n %s|%s|Breaktime:%s' % (ncfile.attrs['test_prod'],
                                                              ncfile.attrs['ref_prod'],
                                                              breaktime_str)
        plt.title(title, fontsize=20)

        plt.annotate(textbox, fontsize=15, xy=(0.025, 0.05), xycoords='axes fraction',
                     bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 3})

        plotfile_name = '%s_%s' % (name_base, breaktime_str)
        plt.savefig(os.path.join(out_path, plotfile_name + '.png'), dpi=f.dpi)
        plt.close()



def compare_RTM_RTG(RTM_dir, RTG_dir, cci_prod, model_prod,
                    filename='GLOBAL_test_results_before_adjustment.nc'):

    times_obj = CCITimes(cci_prod)
    times = times_obj.get_times(as_datetime=False)

    cols = times['breaktimes'].size
    fig, axs = plt.subplots(1, cols, figsize=(cols * 3, 5), facecolor='w', edgecolor='k', sharey=False)

    fig.suptitle('Comparison Break Detection RTM,RTG', fontsize=20, fontweight='bold')
    fig.subplots_adjust(top=0.85, wspace=0.2, hspace=0.1)
    axs = axs.ravel()
    comp_meas = []

    ncfile_rtm = xr.open_dataset(os.path.join(RTM_dir, filename))
    ncfile_rtg = xr.open_dataset(os.path.join(RTG_dir, filename))

    compare_points = {}

    for i, breaktime in enumerate(times['breaktimes']):
        #filename = filename_pattern % breaktime
        compare_points[breaktime] = {'RTM_found_break_RTG_not': [],
                                     'RTG_found_break_RTM_not': []}

        DF_Points_rtg_from_file = ncfile_rtg.sel(time=breaktime).to_dataframe()\
            .reset_index(inplace=False).set_index(['gpi'])


        #DF_Points_rtg_from_file = ncfile_rtg.to_dataframe().reset_index(inplace=False) \
        #    .set_index(['gpi'])

        DF_Points_ismn = DF_Points_rtg_from_file.loc[DF_Points_rtg_from_file.index.dropna()] \
            .rename(columns={'test_results': 'RTG'})
        DF_Points_ismn = DF_Points_ismn[np.isfinite(DF_Points_ismn['RTG'])][['RTG', 'lon', 'lat']]

        DF_Points_rtm_from_file = ncfile_rtm.sel(time=breaktime).to_dataframe()\
            .reset_index(inplace=False).set_index(['gpi'])

        #DF_Points_rtm_from_file = ncfile_rtm.to_dataframe().reset_index(inplace=False) \
        #    .set_index(['gpi'])[['test_results']]

        DF_Points_model = DF_Points_rtm_from_file.reset_index().dropna().set_index('gpi')\
            .rename(columns={'test_results': 'RTM'})

        DF_Points_model = DF_Points_model[np.isfinite(DF_Points_model['RTM'])][['RTM']]

        DF_Points_merged = pd.concat([DF_Points_model, DF_Points_ismn], axis=1) \
            .dropna(how='any')

        DF_Points_merged['diff'] = np.nan
        index = 0
        if not DF_Points_merged.empty:
            for gpi in DF_Points_merged.index.values:
                if DF_Points_merged['RTM'].loc[gpi] == np.nan or DF_Points_merged['RTG'].loc[gpi] == np.nan:
                    continue
                elif DF_Points_merged['RTM'].loc[gpi] == DF_Points_merged['RTG'].loc[gpi]:
                    DF_Points_merged = DF_Points_merged.set_value(gpi, 'diff', 0)
                elif DF_Points_merged['RTM'].loc[gpi] != 4. and DF_Points_merged['RTG'].loc[gpi] == 4.:
                    DF_Points_merged = DF_Points_merged.set_value(gpi, 'diff', 1)
                    compare_points[breaktime]['RTM_found_break_RTG_not'].append(gpi)
                elif DF_Points_merged['RTM'].loc[gpi] == 4. and DF_Points_merged['RTG'].loc[gpi] != 4.:
                    DF_Points_merged = DF_Points_merged.set_value(gpi, 'diff', -1)
                    compare_points[breaktime]['RTG_found_break_RTM_not'].append(gpi)
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
        axs[i].set_title('Break: %s' % breaktime)

    fig.savefig(os.path.join(RTG_dir, 'RTM_vs_RTG'))
    print('Number of points to compare:' + str(comp_meas))
    for point_bt, data_dict in compare_points.iteritems():
        ds = dict_csv_wrapper(data_dict)
        ds.write(os.path.join(RTG_dir, '%s_RTM_vs_RT_conflicting_points.csv' % str(point_bt)))
    return compare_points



def show_processed_gpis(in_file, out_path, plottype, name_base=None):
    '''
    Calculate spatial plots for the areas where Homogeneity Tests were (not)
    performed
    '''
    ncpath = in_file
    if not name_base:
        path = os.path.normpath(in_file)
        filename = path.split(os.sep)[-1]
        name_base = 'STATUS_%s' % in_file

    ccigrid = SMECV_Grid_v042()  # type: CellGrid()
    land_points = ccigrid.get_grid_points()[0]

    ncfile = xr.open_dataset(ncpath)
    for time in ncfile['time']:

        DF_Points = pd.DataFrame(index=range(0, 720 * 1440))
        DF_Points_from_file = ncfile.sel(time=time).to_dataframe().reset_index(inplace=False).set_index(['gpi'])
        # DF_Points_from_file = ncfile.to_dataframe().reset_index(inplace=False).set_index(['location_id'])
        DF_Points_from_file = DF_Points_from_file.loc[DF_Points_from_file.index.dropna()]
        for name in DF_Points_from_file.columns.values:
            DF_Points[name] = DF_Points_from_file[name]

        breaktime_str = str(np.datetime64(np.datetime_as_string(time.values,
                                                                timezone='local')[:10]))

        if plottype == 'test_status':
            statname = 'test_status'
            title = 'HomogeneityTesting Coverage \n %s|%s|%s' % (ncfile.attrs['test_prod'],
                                                                 ncfile.attrs['ref_prod'],
                                                                 breaktime_str)

        else:
            statname = 'adjustment_status'
            title = 'BreakAdjustment Coverage \n %s|%s|%s' % (ncfile.attrs['test_prod'],
                                                              ncfile.attrs['ref_prod'],
                                                              breaktime_str)

        meta = ncfile.variables[statname].attrs

        colors = ["#9b59b6", "#2ecc71", "#3498db", "#95a5a6", "#e74c3c",
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

        f = plt.figure(num=None, figsize=(20, 10), dpi=200, facecolor='w', edgecolor='k')

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
            cb.ax.text((float(j) / float(N)) + float(1. / float(N * 2)), .5, str(j), fontsize=15, ha='center', va='center')

        plt.title(title, fontsize=20)

        stats = DF_Points[statname].dropna().values
        all_gpis = float(land_points.size)

        # Grab the groups and their values from the meta data of the nc file

        textbox = 'Total ground points: %i\n' % all_gpis

        for val in sorted(meta.keys()):
            stat = float(stats[np.where(stats == float(val))].size / all_gpis) * 100.
            textbox += str(val) + ': ' + meta[val] + ': ' + '%.2f%%\n' % stat

        plt.annotate(textbox, fontsize=10, xy=(0.025, 0.05), xycoords='axes fraction',
                     bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 3})

        plotfile_name = '%s_%s' % (name_base,  breaktime_str)

        plt.savefig(os.path.join(out_path, plotfile_name+'.png'), dpi=f.dpi)
        plt.close()



def longest_homog_period_plots(workdir, startyear_plot=True, endyear_plot=True):
    # Create Plot of the longest homogeneouse Period of all nc files in the working dir
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

    log = LogFile(workdir)
    test_prod, ref_prod = log.get_products()

    lons = (np.arange(360 * 4) * 0.25) - 179.875
    lats = (np.arange(180 * 4) * 0.25) - 89.875
    lons, lats = np.meshgrid(lons, lats)

    cmap = plt.cm.jet

    levels = [8, 11.5, 15.5, 21, 30]
    cmaplist = [(0.5, 0, 0, 1), (0.7, 0, 0, 1), (1, 0.2, 0.2, 1), (1, 0.2, 0.2, 1),
                (0, 0.7, 0.7, 1), (0, 0.6, 0.7, 1), (0, 0.6, 0.7, 1), (0, 0.3, 1, 1), (0, 0.3, 1, 1),
                (0, 0, 0.5, 1), (0, 0, 0.5, 1), (0, 0, 0.5, 1)]

    cmap = cmap.from_list('LongestPeriodCMap', cmaplist, cmap.N)

    img = np.empty(DF_Points.index.size, dtype='float32')
    img.fill(np.nan)
    img[DF_Points.index] = DF_Points['max_Period'].values

    # mask array where invalid values (nans) occur
    img_masked = np.ma.masked_invalid(img.reshape((180 * 4, 360 * 4)))

    f = plt.figure(num=None, figsize=(20, 10), dpi=80, facecolor='w', edgecolor='k')

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
    plt.savefig(os.path.join(workdir, filename + '.png'), dpi=f.dpi)

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
            plt.savefig(os.path.join(workdir, filename + '.png'), dpi=f.dpi)
'''
if __name__ == '__main__':
    # df = pd_from_2Dnetcdf(r"H:\HomogeneityTesting_data\output\CCI33\HomogeneityTest_merra2_1991-08-01.nc", return_only_tested=False)
    # extreme_break(r'H:\HomogeneityTesting_data\output\CCI22EGU','gldas_v2','cci_22')
    compare_RTM_RTG(r'D:\users\wpreimes\datasets\HomogeneityTesting_data\CCI41_global_modelsfromdaily_adjustment',
                 r'D:\users\wpreimes\datasets\HomogeneityTesting_data\v37',
                 'CCI_41_COMBINED',
                 'merra2')

    '''
    out_path = r'D:\users\wpreimes\datasets\HomogeneityTesting_data\v34'
    show_processed_gpis(os.path.join(out_path,'GLOBAL_lin_model_params_before_adjustment.nc'),out_path,
                        'adjustment_status', 'ADJUSTMENT_STATUS_before_adjustment')

    inhomo_plot_with_stats(os.path.join(out_path,'GLOBAL_test_results_before_adjustment.nc'),out_path,
                         "TEST_RESULTS_before_adjustment.nc")
    inhomo_plot_with_stats(os.path.join(out_path,'GLOBAL_test_results_after_adjustment.nc'),out_path,
                         "TEST_RESULTS_after_adjustment.nc")

    show_processed_gpis(os.path.join(out_path,'GLOBAL_lin_model_params_after_adjustment.nc'),out_path,
                       'adjustment_status', 'ADJUSTMENT_STATUS_after_adjustment')
    show_processed_gpis(os.path.join(out_path,'GLOBAL_test_results_after_adjustment.nc'),out_path,
                        'test_status', 'TEST_STATUS_after_adjustment')
    show_processed_gpis(os.path.join(out_path,'GLOBAL_test_results_before_adjustment.nc'),out_path,
                        'test_status', 'TEST_STATUS_before_adjustment')
    '''
    # longest_homog_period_plots(r"D:\users\wpreimes\datasets\HomogeneityTesting_data\output\CCI41_ADJUST_global")
    #show_processed_gpis('/data-write/USERS/wpreimes/HomogeneityTesting_data/v1','test',"HomogeneityTestResult_2007-01-01_image.nc")
    #inhomo_plot_with_stats(r'D:\users\wpreimes\datasets\HomogeneityTesting_data\CCI41_global_modelsfromdaily_adjustment',
    #                     "GLOBAL_test_results_before_adjustment.nc")
    # compare_RTM_RTG(r'D:\users\wpreimes\datasets\HomogeneityTesting_data\output\CCI41_noAdjust',
    #                r'D:\users\wpreimes\datasets\HomogeneityTesting_data\output\CCI41_ISMN',
    #                'CCI_41_COMBINED',
    #                'merra2')
    # show_tested_gpis(r"H:\HomogeneityTesting_data\output\v15","HomogeneityTest_merra2_2011-10-01")
    # valid_months_plot(r'H:\HomogeneityTesting_data\output\CCI_available_data_plots', 'CCI_33_COMBINED','merra2', 'M', 0.33)
