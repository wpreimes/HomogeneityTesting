# -*- coding: utf-8 -*-
"""
Created on Jän 18 16:26 2018

@author: wpreimes
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

# original_breaks
def adjustment_ts_plot(dataframe, breaktimes, gpi, refdata_name, gpi_test_results, resample='M',
                       title=None, save_path=None, ax=None, legend=False, stats=False):
    '''
    Saves a plot of the passed data from the adjustment process
    :param dataframe: pd.DataFrame
        Data in DF, must contain: testdata (adjusted), testdata_original (original ccism),
                                  refdata_original (original reference)
    :param breaktimes: list
        List of breaktimes for the ccism product
    :param gpi: int
        GPI of the tested point
    :param ref_name: str
        Name of the reference data set
    :param gpi_test_result: dict
        Test results for the gpi, weather a break was detected or not
    :return:
        None
    '''

    if not ax:
        fig_full_ts, ax = plt.subplots(figsize=(20, 6))


    plot_df = dataframe[['refdata', 'testdata_before_adjustment', 'testdata_after_adjustment']].resample(resample).mean()


    not_adjused_values = plot_df.loc[plot_df['testdata_after_adjustment'] == plot_df['testdata_before_adjustment']].index
    plot_df.loc[not_adjused_values, 'testdata_after_adjustment'] = np.nan
    #plot_df = plot_df.rename(columns={'refdata': refdata_name,
    #                                  'testdata_after_adjustment': 'CCI SM adjusted',
    #                                  'testdata_before_adjustment': 'CCI SM unadjusted'})
    plot_df = plot_df.rename(columns={'refdata': refdata_name,
                                      'testdata_after_adjustment': 'CCI SM adjusted',
                                      'testdata_before_adjustment': 'CCI SM in'})
    plot_df.plot(title=title,
                 style=['k-', 'r-', 'b--'],
                 ax=ax,
                 legend=legend)
    ax.set_xlabel('Years')
    ax.set_ylabel('SM')
    for breaktime in breaktimes:
        if breaktime in gpi_test_results['homogeneous_adjusted']:
            color = 'green'
            linestyle = 'solid'
        elif breaktime in gpi_test_results['homogeneous_original']:
            color = 'green'
            linestyle = 'dashed'
        elif breaktime in gpi_test_results['inhomogeneous_adjusted']:
            color = 'red'
            linestyle = 'solid'
        else:
            color = 'grey'
            linestyle = 'dashed'

        ax.axvline(breaktime, linestyle=linestyle, color=color, lw=2)

        if stats:
            median_p1_adjusted = np.median(plot_df.loc[:breaktime, 'CCI SM adjusted'])
            median_p2_adjusted = np.median(plot_df.loc[breaktime + pd.DateOffset(1), 'CCI SM adjusted'])
            var_p1_adjusted = np.var(plot_df.loc[:breaktime, 'CCI SM adjusted'])
            var_p2_adjusted = np.var(plot_df.loc[breaktime + pd.DateOffset(1), 'CCI SM adjusted'])
            varratio_adjusted = var_p1_adjusted / var_p2_adjusted
            mediandiff_adjusted = median_p1_adjusted - median_p2_adjusted
            #TODO: Add the stats to the time frames plots

    if save_path:
        if not os.path.isdir(save_path):
            os.mkdir(save_path)
        fig_full_ts.savefig(os.path.join(save_path, '%i_adjusted_full.png' % gpi))



