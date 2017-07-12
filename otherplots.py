# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 17:06:28 2017

@author: wpreimes
"""
import numpy as np
import math
import os, glob
from datetime import datetime, timedelta
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from HomogeneityTesting.import_satellite_data import QDEGdata_M
from HomogeneityTesting.points_to_netcdf import globalCellgrid, points_to_netcdf
from HomogeneityTesting.otherfunctions import regress, calc_longest_homogeneous_period
from HomogeneityTesting.save_data import load_Log
from HomogeneityTesting.nlcmap import nlcmap

from scipy import stats

import pygeogrids.netcdf as nc

from points_to_netcdf import pd_from_2Dnetcdf
from netCDF4 import Dataset

def inhomo_plot_with_stats(workdir):
    # type: (str) -> None
    '''
    :param workdir: path to directory containing nc files from HomogeneityTesting
    :return: None
    '''
    fileslist = glob.glob(os.path.join(workdir, "HomogeneityTest*.nc"))

    landgrid = nc.load_grid(r"D:\users\wpreimes\datasets\grids\qdeg_land_grid.nc")
    land_points = landgrid.get_grid_points()[0]

    for filename in fileslist:
        DF_Points = pd_from_2Dnetcdf(filename, return_only_landpoints=False)

        filename = filename.replace(workdir, '')
        filename = filename.replace('.', '_')
        splitname = filename.split('_')

        log = load_Log(workdir)
        products = log.get_products()
        ref_prod = products['ref_prod']

        breaktime = datetime.strptime(splitname[2], '%Y-%m-%d')

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

        img = np.empty(DF_Points['lon'].values.size, dtype='float32')
        img.fill(np.nan)
        img[DF_Points.index.values] = DF_Points['test_results'].values

        # mask array where invalid values (nans) occur
        img_masked = np.ma.masked_invalid(img.reshape((180 * 4, 360 * 4)))

        f = plt.figure(num=None, figsize=(20, 10), dpi=90, facecolor='w', edgecolor='k')

        m = Basemap(projection='mill', llcrnrlat=-60, urcrnrlat=80,
                    llcrnrlon=-180, urcrnrlon=180, resolution='c', )

        m.drawcoastlines()
        m.drawcountries()
        lons, lats = np.meshgrid(np.unique(DF_Points['lon'].values),
                                 np.unique(DF_Points['lat'].values))

        im = m.pcolormesh(lons, np.flipud(lats), img_masked, cmap=cmap, latlon=True)

        im.set_clim(vmin=1, vmax=4)

        cb = m.colorbar(im, "bottom", size="5%", pad="8%")
        cb.outline.set_visible(False)
        cb.set_ticks([])

        for t in cb.ax.get_xticklabels():
            t.set_fontsize(20)

        cb.set_label('Break detected by Test:', fontsize=20, labelpad=-50)
        for j, lab in enumerate(['Wilkoxon', 'Fligner-Killeen', 'Both', 'None']):
            cb.ax.text((2 * j + 1) / 8.0, .5, lab, fontsize=15, ha='center', va='center')

        title = 'HomogeneityTesting \n %s|%s|Breaktime:%s' % ('CCI', ref_prod, breaktime.strftime("%Y-%m-%d"))
        plt.title(title, fontsize=20)

        plt.annotate(textbox, fontsize=15, xy=(0.025, 0.05), xycoords='axes fraction',
                     bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 3})

        filename = 'HomogeneityTest_%s_(breaktime-%s)' % (ref_prod, breaktime.strftime("%Y-%m-%d"))
        plt.savefig(workdir + '\\' + filename + '.png', dpi=f.dpi)


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


def show_tested_gpis(workdir):
    '''
    Calculate spatial plots for the areas where Homogeneity Tests were (not)
    performed
    '''
    fileslist = glob.glob(os.path.join(workdir, "HomogeneityTest*.nc"))

    landgrid = nc.load_grid(r"D:\users\wpreimes\datasets\grids\qdeg_land_grid.nc")
    land_points = landgrid.get_grid_points()[0]

    log = load_Log(workdir)
    products = log.get_products()
    ref_prod = products['ref_prod']
    test_prod = products['test_prod']

    for filename in fileslist:
        DF_Points = pd_from_2Dnetcdf(filename, return_only_landpoints=False)

        status_var = Dataset(filename).variables['status']

        filename = filename.replace(workdir, '')
        filename = filename.replace('.', '_')
        splitname = filename.split('_')
        breaktime = datetime.strptime(splitname[2], '%Y-%m-%d')

        colors = ["#2ecc71", "#9b59b6", "#3498db", "#95a5a6", "#e74c3c",
                  "#34495e", "#FFC200", "#FFC0CB", "#FF1493", "#ffffff"]

        N = len(colors)
        cmap = LinearSegmentedColormap.from_list('status',
                                                 colors,
                                                 N=N)

        img = np.empty(DF_Points['lon'].values.size, dtype='float32')
        img.fill(np.nan)
        img[DF_Points.index.values] = DF_Points['status'].values

        # mask array where invalid values (nans) occur
        img_masked = np.ma.masked_invalid(img.reshape((180 * 4, 360 * 4)))

        f = plt.figure(num=None, figsize=(20, 10), dpi=90, facecolor='w', edgecolor='k')

        m = Basemap(projection='mill', llcrnrlat=-60, urcrnrlat=80, \
                    llcrnrlon=-180, urcrnrlon=180, resolution='c', )

        m.drawcoastlines()
        m.drawcountries()

        lons, lats = np.meshgrid(np.unique(DF_Points['lon'].values),
                                 np.unique(DF_Points['lat'].values))

        im = m.pcolormesh(lons, np.flipud(lats), img_masked, cmap=cmap, latlon=True)

        im.set_clim(vmin=0, vmax=N)

        cb = m.colorbar(im, "bottom", size="5%", pad="8%")
        cb.outline.set_visible(True)
        cb.set_ticks([])

        for t in cb.ax.get_xticklabels():
            t.set_fontsize(20)

        cb.set_label('Break detected by Test:', fontsize=20, labelpad=-50)

        for j in range(N):
            cb.ax.text((float(j)/float(N)) + float(1./float(N*2)), .5, str(j), fontsize=15, ha='center', va='center')

        title = 'HomogeneityTesting Coverage \n %s|%s|%s' % (test_prod, ref_prod, breaktime.strftime("%Y-%m-%d"))
        plt.title(title, fontsize=20)

        stats = DF_Points['status'].dropna().values
        all_gpis = float(land_points.size)

        # Grab the groups and their values from the meta data of the nc file
        meta ={}
        groups_string = status_var.Values
        groups_strings = groups_string.split(',')
        for string in groups_strings:
            string = string.split('=')
            meta.update({int(string[0]): string[1]})

        textbox = 'Total ground points: %i\n' % all_gpis

        for val, s in meta.iteritems():
            stat = float(stats[np.where(stats == float(val))].size / all_gpis) * 100.
            textbox += str(val)+': ' + s + ': ' + '%.2f%%\n' % stat

        plt.annotate(textbox, fontsize=10, xy=(0.025, 0.05), xycoords='axes fraction',
                     bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 3})

        filename = 'RT_coverage_%s_%s' % (ref_prod, breaktime.strftime("%Y-%m-%d"))
        plt.savefig(workdir + '\\' + filename + '.png', dpi=f.dpi)



'''
def longest_homog_period_plots(workdir, create_start_year_plot=True, create_end_year_plot=True):
    
    #Create Plot of the longest homogeneouse Period of all nc files in the working dir
    #:param workdir: str
    #    Path to the directory containing the nc files returned by Homogeneity Test
    #:param same_colorbar: bool
    #    True to use fixed colobar (for comparsion of results for differenct cci versions
    #    False to use a different colorbar for different cci versions
    #:return: None
    
    DF_Period = calc_longest_homogeneous_period(workdir)

    if create_end_year_plot:
        points_to_netcdf(DF_Period[['endyear']], workdir, None, 'startendyears', None, None)
    if create_start_year_plot:
        points_to_netcdf(DF_Period[['startyear']], workdir, None, 'startendyears', None, None)



    lons = (np.arange(360 * 4) * 0.25) - 179.875
    lats = (np.arange(180 * 4) * 0.25) - 89.875
    lons, lats = np.meshgrid(lons, lats)

    cmap = plt.cm.jet

    levels = [8, 11.5, 15.5, 21]
    cmaplist = [(0.5, 0, 0, 1), (0.7, 0, 0, 1), (1, 0.2, 0.2, 1), (1, 0.2, 0.2, 1),
                (0, 0.7, 0.7, 1), (0, 0.6, 0.7, 1), (0, 0.6, 0.7, 1), (0, 0.3, 1, 1),
                (0, 0, 0.5, 1), (0, 0, 0.5, 1), (0, 0, 0.5, 1)]

    cmap = cmap.from_list('LongestPeriodCMap', cmaplist, cmap.N)
    cmap_nonlin = nlcmap(cmap, levels)
    # cmaplist = [cmap(i) for i in range(cmap.N)]


    img = np.empty(lons.size, dtype='float32')
    img.fill(np.nan)
    img[DF_Period.index.values] = DF_Period['max_Period'].values

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
    # cb.outline.set_visible(False)
    # cb.set_ticks(cmap_nonlin.transformed_levels)

    cb.set_ticklabels(["%.1f" % lev for lev in levels])
    for t in cb.ax.get_xticklabels():
        t.set_fontsize(20)

    cb.set_label('Period Length [years]:', fontsize=20, labelpad=-75)

    title = 'Length of Longest Homogeneous Period\n %s|%s' % ('CCI', ref_prod)
    plt.title(title, fontsize=20)

    # plt.annotate(textbox, fontsize=15,xy=(0.025, 0.05), xycoords='axes fraction',bbox={'facecolor':'grey', 'alpha':0.5, 'pad':3})

    filename = 'RTM_LongestHomogeneousePeriod_%s' % (ref_prod)
    plt.savefig(workdir + '\\' + filename + '.png', dpi=f.dpi)
'''

def calc_longest_homogeneous_period(workdir, test_prod, ref_prod, same_colorbar=False):
    # TODO: This function is bad...
    fileslist = glob.glob(os.path.join(workdir, "DF_Points_%s*.csv" % ref_prod))
    starttimes = []
    endtimes = []
    breaktimes = []
    periodsizes = []
    for filename in fileslist:
        filename = filename.replace(workdir, '')
        filename = filename.replace('.', '_')
        splitname = filename.split('_')
        starttimes.append(datetime.strptime(splitname[3], '%Y-%m-%d'))
        endtimes.append(datetime.strptime(splitname[5], '%Y-%m-%d'))
        breaktimes.append(datetime.strptime(splitname[4], '%Y-%m-%d'))

    DF_Period = pd.concat([pd.read_csv(f, index_col=0, usecols=[0, 6], header=0,
                                       names=['gpi_quarter', breaktimes[i]]) for i, f in enumerate(fileslist)], axis=1)

    i = 1
    for starttime, breaktime, endtime in zip(starttimes, breaktimes, endtimes):
        # DF_Period['%i.1'%i] = np.where(DF_Period[breaktime]==4.,(endtime-starttime).days, ((breaktime-starttime).days),(endtime-breaktime).days))
        DF_Period['%i' % i] = (breaktime - starttime).days
        periodsizes.append((breaktime - starttime).days)
        # DF_Period['%i.2'%i] = (endtime-breaktime).days
        i += 1
    DF_Period['%i' % i] = (endtimes[-1] - breaktimes[-1]).days
    periods = map(list, DF_Period[DF_Period.columns.values[-i:]].values)
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

    for periods in results:
        if max(periods) in periodsizes:
            max_period.append(np.nan)
            startdays.append(np.nan)

        else:
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
            startyears.append(startdate.year)
            endyears.append(enddate.year)

    startdates = np.asarray(startdates)
    enddates = np.asarray(enddates)
    max_period = np.asarray(max_period) / 365.
    DF_Period['max_Period'] = max_period
    DF_Period['startdate'] = startdates
    DF_Period['endate'] = enddates
    DF_Period['startyear'] = startyears
    DF_Period['endyear'] = endyears

    points_to_netcdf(DF_Period[['endyear']], workdir, None, 'startendyears', None, None)
    points_to_netcdf(DF_Period[['startyear']], workdir, None, 'startendyears', None, None)
    lons = (np.arange(360 * 4) * 0.25) - 179.875
    lats = (np.arange(180 * 4) * 0.25) - 89.875
    lons, lats = np.meshgrid(lons, lats)

    cmap = plt.cm.jet

    '''
    if test_prod=='cci_22':
        cbrange=(8,21)
    elif test_prod=='cci_31':
        cbrange=(6,20)
    '''

    if (test_prod[:6] == 'cci_31') or (same_colorbar == True):
        levels = [8, 11.5, 15.5, 21]
        cmaplist = [(0.5, 0, 0, 1), (0.7, 0, 0, 1), (1, 0.2, 0.2, 1), (1, 0.2, 0.2, 1),
                    # (0.3,0.3,0.3,1),(0.5,0.5,0.5,1),(0.5,0.5,0.5,1),
                    (0, 0.7, 0.7, 1),
                    (0, 0.6, 0.7, 1), (0, 0.6, 0.7, 1), (0, 0.3, 1, 1), (0, 0, 0.5, 1), (0, 0, 0.5, 1), (0, 0, 0.5, 1)]
        # cmaplist=[(0.5,0,0,1),(0.7,0,0,1),(1,0.2,0.2,1),(1,0.2,0.2,1),
        # (0,0.7,0.7,1),(0,0.7,0.7,1),
        # (0,0.6,0.7,1),(0,0.3,1,1),(0,0,0.8,1),(0,0,0.8,1),(0,0,0.5,1)]
    elif test_prod[:6] == 'cci_22':
        cmaplist = [(0.5, 0, 0, 1), (0.7, 0, 0, 1), (1, 0.2, 0.2, 1), (1, 0.2, 0.2, 1),
                    # (0.3,0.3,0.3,1),(0.5,0.5,0.5,1),(0.5,0.5,0.5,1),
                    (0, 0.7, 0.7, 1),
                    (0, 0.6, 0.7, 1), (0, 0.6, 0.7, 1), (0, 0.6, 0.7, 1), (0, 0.3, 1, 1), (0, 0, 0.8, 1),
                    (0, 0, 0.8, 1), (0, 0, 0.5, 1), (0, 0, 0.5, 1), (0, 0, 0.5, 1)]
        levels = [8, 10, 15.5, 21, 24]

    cmap = cmap.from_list('LongestPeriodCMap', cmaplist, cmap.N)
    # cmap_nonlin = nlcmap(cmap, levels)
    # cmaplist = [cmap(i) for i in range(cmap.N)]


    img = np.empty(lons.size, dtype='float32')
    img.fill(np.nan)
    img[DF_Period.index.values] = DF_Period['max_Period'].values

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
    # cb.outline.set_visible(False)
    # cb.set_ticks(cmap_nonlin.transformed_levels)

    cb.set_ticklabels(["%.1f" % lev for lev in levels])
    for t in cb.ax.get_xticklabels():
        t.set_fontsize(20)

    cb.set_label('Period Length [years]:', fontsize=20, labelpad=-75)

    title = 'Length of Longest Homogeneous Period\n %s|%s' % ('CCI', ref_prod)
    plt.title(title, fontsize=20)

    # plt.annotate(textbox, fontsize=15,xy=(0.025, 0.05), xycoords='axes fraction',bbox={'facecolor':'grey', 'alpha':0.5, 'pad':3})

    filename = 'RTM_LongestHomogeneousePeriod_%s' % (ref_prod)
    plt.savefig(workdir + '\\' + filename + '.png', dpi=f.dpi)

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

# df = pd_from_2Dnetcdf(r"H:\workspace\HomogeneityTesting\output\CCI33\HomogeneityTest_merra2_1991-08-01.nc", return_only_tested=False)
# extreme_break(r'H:\workspace\HomogeneityTesting\output\CCI22EGU','gldas_v2','cci_22')
# inhomo_plot_with_stats(r"H:\workspace\HomogeneityTesting\output\CCI33")
#calc_longest_homogeneous_period(r"H:\workspace\HomogeneityTesting\output\CCI31EGU",'cci_22','merra2')
#show_tested_gpis(r"H:\workspace\HomogeneityTesting\output\CCI33")
# inhomo_plot_with_stats(r'H:\workspace\HomogeneityTesting\output\v32')
# compare_RTM_RTG(r'H:\workspace\HomogeneityTesting\output\CCI31EGU','merra2')
# show_tested_gpis(r"H:\workspace\HomogeneityTesting\output\v30",'gldas-merged')
