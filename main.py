# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 17:05:55 2017

@author: wpreimes
"""

from typing import Union
import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime
from multiprocessing import Process, Queue
from pynetcf.time_series import GriddedNcIndexedRaggedTs
from otherfunctions import split, create_workfolder
from interface import BreakTestBase, get_test_meta
from cci_timeframes import CCITimes
from save_data import RegularGriddedCellData, LogFile
from grid_functions import split_cells
from adjustment import LinearAdjustment, get_adjustment_meta
from otherplots import show_processed_gpis, inhomo_plot_with_stats, longest_homog_period_plots, ts_adjustment_plot
from smecv_grid.grid import SMECV_Grid_v042
import warnings
import pandas as pd
from timeseries_plots import adjustment_ts_plot


def timeseries_stats(dataframe, breaktime, timeframe, return_ts=False):
    #TODO: Move this somewhere

    df = dataframe.copy(True)
    data_dict = {}
    breaktime_next = breaktime + pd.DateOffset(1)
    # Mean of adjusted testdata after breaktime
    mean_after_break_testdata = np.nanmean(df.loc[breaktime_next:timeframe[1], 'testdata'].values)
    # Mean of adjusted testdata before breaktime
    mean_before_break_testdata = np.nanmean(df.loc[timeframe[0]:breaktime, 'testdata'].values)
    # Mean of unadjusted testdata after break
    mean_after_break_original = np.nanmean(df.loc[breaktime_next:timeframe[1], 'testdata_original'].values)
    # Mean of testdata before break
    mean_before_break_original = np.nanmean(df.loc[timeframe[0]:breaktime, 'testdata_original'].values)
    # Mean of refdata after break
    mean_after_break_refdata_original = np.nanmean(df.loc[breaktime_next:timeframe[1], 'refdata_original'].values)
    # Mean of refdata before break
    mean_before_break_refdata_original = np.nanmean(df.loc[timeframe[0]:breaktime, 'refdata_original'].values)
    # Mean of refdata after break
    if 'refdata' in dataframe.columns.values:
        mean_after_break_refdata= np.nanmean(df.loc[breaktime_next:timeframe[1], 'refdata'].values)
        # Mean of refdata before break
        mean_before_break_refdata = np.nanmean(df.loc[timeframe[0]:breaktime, 'refdata'].values)


    # Variance of adjusted testdata
    var_after_break_testdata = np.nanvar(df.loc[breaktime_next:timeframe[1], 'testdata'].values)
    var_before_break_testdata = np.nanvar(df.loc[timeframe[0]:breaktime, 'testdata'].values)
    # Variance of unadjusted testdata
    var_after_break_original = np.nanvar(df.loc[breaktime_next:timeframe[1], 'testdata_original'].values)
    var_before_break_original = np.nanvar(df.loc[timeframe[0]:breaktime, 'testdata_original'].values)
    # Variance of original refdata after break
    var_after_break_refdata_original = np.nanvar(df.loc[breaktime_next:timeframe[1], 'refdata_original'].values)
    var_before_break_refdata_original = np.nanvar(df.loc[timeframe[0]:breaktime, 'refdata_original'].values)
    if 'refdata' in dataframe.columns.values:
        # Variance of refdata after break
        var_after_break_refdata = np.nanvar(df.loc[breaktime_next:timeframe[1], 'refdata'].values)
        var_before_break_refdata = np.nanvar(df.loc[timeframe[0]:breaktime, 'refdata'].values)




    data_dict['mean_after_break_testdata'] = mean_after_break_testdata
    data_dict['mean_before_break_testdata'] = mean_before_break_testdata
    data_dict['mean_after_break_original'] = mean_after_break_original
    data_dict['mean_before_break_original'] = mean_before_break_original
    data_dict['mean_after_break_refdata_original'] = mean_after_break_refdata_original
    data_dict['mean_before_break_refdata_original'] = mean_before_break_refdata_original
    if 'refdata' in dataframe.columns.values:
        data_dict['mean_after_break_refdata'] = mean_after_break_refdata
        data_dict['mean_before_break_refdata'] = mean_before_break_refdata

    data_dict['var_after_break_testdata'] = var_after_break_testdata
    data_dict['var_before_break_testdata'] = var_before_break_testdata
    data_dict['var_after_break_original'] = var_after_break_original
    data_dict['var_before_break_original'] = var_before_break_original
    data_dict['var_after_break_refdata_original'] = var_after_break_refdata_original
    data_dict['var_before_break_refdata_original'] = var_before_break_refdata_original
    if 'refdata' in dataframe.columns.values:
        data_dict['var_after_break_refdata'] = var_after_break_refdata
        data_dict['var_before_break_refdata'] = var_before_break_refdata


    if return_ts:
        Q_adj_ref_after_break = df.loc[breaktime_next:timeframe[1], 'testdata'].values - \
                                df.loc[breaktime_next:timeframe[1], 'refdata_original'].values
        Q_orig_ref_after_break = df.loc[breaktime_next:timeframe[1], 'testdata_original'].values - \
                                 df.loc[breaktime_next:timeframe[1], 'refdata_original'].values
    else:
        df = False

    return df, data_dict


def process_for_cells(q, workfolder, skip_times, times_obj, save_objects, cells, log_file, test_prod, ref_prod, anomaly,
                      process_no, save_adjusted_data=True):
    # Function for multicore processing

    print('Start process %i' % process_no)

    ########################################## TODO: Values to change
    if ref_prod == 'ISMN-Merge':
        min_data_for_temp_resampling = 0.1
        min_data_size = 3
        spearman_check = {'r': 0, 'p': 0.1}
    else:
        min_data_for_temp_resampling = 0.3  # If CCI data contains less values, the monthly value is resampled as nan
        min_data_size = 5  # Minimum number of monthly values before/after break to perform homogeneity testing # TODO: different value?
        spearman_check = {'r': 0, 'p': 0.05}

    alpha = 0.01
    adjusted_data_path = os.path.join(workfolder, test_prod + '_ADJUSTED')
    resample_method = 'M'  # or D
    refdata_correction_for = 'timeframe'  # 'timeframe', 'iteration' or 'gpi'
    tests = {'mean': 'wilkoxon', 'var': 'scipy_fligner_killeen'}
    backward_extended_timeframes = True
    model_plots = True
    ts_plots = (True, 'M')
    quantile_filtering_increase = 0.01  # Increase of quantile per iteration (0 ,0.03, 0.04,... and 1, 0.97, 0.94,....]
    first_iter_mode = 'both'  # 'd' (BAD!!!) or 'both' for adjustment if first iteration found wk break
    max_retries = 15  # max number of iterations if break is still found after adjustment
    adjust_always = False
    priority = False  # or ['mean', 'var']
    models_from = 'daily'
    remove_data_until_iter = 10  # Remove bad values for each iteration until this one, than stay at the amount
    filter_parts = 'both'
    adjust_all_values_before_break = True
    process_logs = False
    ##########################################

    if process_no == 0:
        log_file.add_line('-------------------------------')
        log_file.add_line(
            '%s: %s' % ('alpha - Significance level for rejecting the Null-Hypothesis (for finding a break)', alpha))
        log_file.add_line('%s: %s' % ('adjusted data saved to', adjusted_data_path))
        log_file.add_line('%s: %s' % ('min percent. of valid data for resampling', min_data_for_temp_resampling))
        log_file.add_line('%s: %s' % ('min number of point pairs for modelling', min_data_size))
        log_file.add_line('%s: %s' % ('parameters for checking correlaton between refdata and testdata',
                                      spearman_check))
        log_file.add_line('%s: %s' % ('temporal merging resolution', resample_method))
        log_file.add_line('%s: %s' % ('break test, which are performed to detect certain kinds of breaks', tests))
        log_file.add_line('%s: %s' % ('Fitting of reference data to test data', refdata_correction_for))
        log_file.add_line(
            '%s: %s' % ('Extending time frames over homogeneous, tested break times', backward_extended_timeframes))
        log_file.add_line('%s: %s' % (
            'Increase rate per iteration for quantile filtering of difference values', quantile_filtering_increase))
        log_file.add_line('%s: %s' % ('adjustment method for first iteration', first_iter_mode))
        log_file.add_line('%s: %s' % ('Maximum number of adjustment iterations', max_retries))
        log_file.add_line('%s: %s' % ('Adjust always, also if no break was found', adjust_always))
        log_file.add_line('%s: %s' % ('Prioritize kind of break to be removed (allows lower prioritzed to be added)',
                                      priority))
        log_file.add_line('%s: %s' % ('create linear models with values', models_from))
        log_file.add_line('%s: %s' % ('Perform quantile filtering until iteration', remove_data_until_iter))
        log_file.add_line('%s: %s' % ('Filter both parts for first iter and for other iters:', filter_parts))
        log_file.add_line('%s: %s' % ('adjust whole time series before the break time for each found break',
                                      adjust_all_values_before_break))

        log_file.add_line('-------------------------------')

    if process_logs:
        process_log = LogFile(workfolder, 'process_log_%i.txt' % process_no,
                              {'Process No.': process_no,
                               'Started at': datetime.now().strftime('%Y-%m-%d, %H:%M')})

    test_obj = BreakTestBase(test_prod,
                             ref_prod,
                             tests,
                             alpha,  # TODO: Choose higher alpha value eg 0.05 or 0.1
                             anomaly)

    grid = SMECV_Grid_v042()

    if save_adjusted_data:
        dataset = GriddedNcIndexedRaggedTs(path=adjusted_data_path,
                                           grid=grid, mode='w')
        max_retries = max_retries
        min_retries = 1 if save_adjusted_data == 'always' else 0
    else:
        dataset = None
        max_retries = 1
        min_retries = 0

    if model_plots:
        model_plots_dir = os.path.join(workfolder, 'model_plots')
        if not os.path.isdir(model_plots_dir):
            os.mkdir(model_plots_dir)

    bad_gpis = []  # GPIS where break was not removed
    for icell, cell in enumerate(cells):

        print('Processing QDEG Cell %i (iteration %i of %i)' % (cell, icell + 1, len(cells)))
        if process_logs: process_log.add_line('Start Testing Cell: %i at %s' % (cell, datetime.now().strftime('%Y-%m-%d %H:%M')))

        for itergpi, gpi in enumerate(grid.grid_points_for_cell(cell)[0]):
            if gpi not in [391440]: continue
            if process_logs:
                process_log.add_line('-----------------------------------------------------------------------------',1)
                process_log.add_line('processing gpi %i (%i of %i)'
                                     % (gpi, itergpi, grid.grid_points_for_cell(cell)[0].size), 1)

            try:  # load testdata and reference data
                df_time = test_obj.read_gpi(gpi,
                                            start=test_obj.range[0],
                                            end=test_obj.range[1])

                df_time['testdata'] = df_time['testdata_original'].copy(deep=True)
                df_time['testdata_RAW'] = df_time['testdata_original'].copy(deep=True)
            except:
                if process_logs: process_log.add_line('ERROR: loading data failed', 1)
                continue

            if process_logs: process_log.add_line('reference data bias correction', 1)
            corrected_refdata = test_obj.ref_data_correction(df_time.loc[:, ['testdata', 'refdata_original']],
                                                             'refdata_original',
                                                             ignore_exceptions=True)

            df_time.loc[:, 'refdata_original'] = corrected_refdata

            times = times_obj.get_times(gpi, as_datetime=True)
            timeframes, breaktimes = times['timeframes'], times['breaktimes']

            gpi_test_results = {'homogeneous_original':[],
                                'homogeneous_adjusted':[],
                                'inhomogeneous_adjusted':[]} # List of tested and approved break times for extending time frames
            if ts_plots[0]:
                ts_frames_fig = plt.figure(figsize=(20, 10))

            for i, (timeframe, breaktime) in enumerate(zip(timeframes, breaktimes)):
                #TODO: !!!!!!!!!!!!!!!!!!!!
                #if breaktime in [datetime(2012,7,1), datetime(2007,1,1)]: continue
                # if breaktime == datetime(1991,8,5):

                lin_model_params_before_adjustment = {'adjustment_status': 0}
                lin_model_params_after_adjustment = {'adjustment_status': 0}
                test_results_before_adjustment = {'test_status': 0}
                test_results_after_adjustment = {'test_status': 0}

                retries = 0
                # In case that after adjustment a break is still found, re-correct the reference data to
                # the (insufficiently) corrected test data and repeat adjustment (< max_retries times)
                while (retries <= min_retries) or (retries < max_retries):
                    log_level = retries + 2
                    if process_logs: process_log.add_line('process timeframe %i (%s to %s), retry %i'
                                                          %(i, str(timeframe[0].date()), str(timeframe[1].date()), retries), log_level)


                    if backward_extended_timeframes:  # Cut data to timeframe
                        while timeframe[1] in gpi_test_results['homogeneous_adjusted'] or \
                              timeframe[1] in gpi_test_results['homogeneous_original']:  # Use end of previous time frame
                            timeframe[1] = times_obj.get_adjacent(gpi, timeframe, -1)[1]
                            print('expand tf to %s' %str(timeframe[1].date()))

                    data_daily = df_time[timeframe[0]:timeframe[1]].copy()  # change with adjustment iterations

                    try:  # Homogeneity Testing
                        # TODO: Fit refdata again to testdata, but only for the timeframe that is tested
                        if (refdata_correction_for == 'timeframe' and retries == 0) or refdata_correction_for == 'timeframe_iteration':
                            if process_logs: process_log.add_line('reference data bias correction for timeframe', log_level)
                            corrected_refdata = test_obj.ref_data_correction(data_daily.loc[timeframe[0]:timeframe[1],
                                                                             ['testdata', 'refdata_original']],
                                                                             'refdata_original')
                            data_daily.loc[timeframe[0]:timeframe[1], 'refdata'] = corrected_refdata

                        elif refdata_correction_for == 'gpi' and retries == 0:
                            data_daily.loc[timeframe[0]:timeframe[1], 'refdata'] = data_daily.loc[timeframe[0]:timeframe[1],
                                                                                   'refdata_original']
                        elif retries > 0:
                            pass
                        else:
                            raise Exception("Choose 'gpi', 'timeframe', or 'timeframe_iteration' for refdata correction")

                        data_daily = data_daily.dropna(subset=['testdata', 'refdata'])
                        # Calculate difference TimeSeries

                        data_daily['Q'] = data_daily['testdata'] - data_daily['refdata']

                        if process_logs: process_log.add_line('monthly resampling', log_level)
                        data_resampled, resample_info = test_obj.temp_resample(data_daily,
                                                                               'M',
                                                                               min_data_for_temp_resampling)

                        # Checks if there is data left and calculates spearman correlation

                        corr, pval = test_obj.check_corr(data_resampled, spearman_check['r'], spearman_check['p'])

                        data_resampled, _, _ = test_obj.group_by_breaktime(data_resampled,
                                                                           breaktime,
                                                                           min_data_size=min_data_size,
                                                                           ignore_exception=False)

                        # Run Tests
                        if process_logs: process_log.add_line('running break tests', log_level)
                        test_result = test_obj.run_tests(data=data_resampled)

                        # TODO: Check this for a timeframe with a break
                        failed_tests, break_found_by, found_break = test_obj.check_testresult(test_result)
                        if failed_tests and process_logs: process_log.add_line('failed tests: %s' %str(failed_tests),
                                                                              log_level)
                        if found_break and process_logs: process_log.add_line('break found by: %s' %str(break_found_by),
                                                                             log_level)

                    except Exception as e:
                        # Testing Failed
                        adjust_obj = None
                        print '%s: Testing failed: %s' % (str(breaktime.date()), e)
                        test_result = {'test_status': int(str(e)[0])}
                        test_results_before_adjustment.update(test_result)
                        if process_logs: process_log.add_line('ERROR: %s' %str(e), log_level)
                        break

                    try:
                        # TODO: Filter only for first iteration, for next iterations use residuals to remove bad values

                        lq, uq = 0.05, 0.95
                        if retries <= remove_data_until_iter:  # Remove more values where ref and test is different for first 5 iterations
                            mult = retries
                        else:  # For later iterations dont remove more but remove as much as before
                            mult = remove_data_until_iter
                        lq = lq + mult * quantile_filtering_increase
                        uq = uq - mult * quantile_filtering_increase

                        if process_logs: process_log.add_line('Quantile Filtering: LQ (%f), UQ (%f)' %(lq, uq), log_level)

                        if models_from == 'daily':
                            adjustment_data = data_daily[['testdata','refdata','Q']]
                        else:
                            adjustment_data = data_resampled[['testdata','refdata','Q']]

                        adjust_obj = LinearAdjustment(data=adjustment_data,
                                                      breaktime=breaktime,
                                                      adjust_part='first',
                                                      adjust_param='both' if retries > 0 else first_iter_mode,
                                                      filter_parts='both' if (filter_parts and retries == 0) else filter_parts,  # if retries==0 else 'first',
                                                      filter_quantiles=(lq,uq),
                                                      model_plots=(retries, adjust_obj.model_plots) if (
                                                              retries > 0 and adjust_obj) else model_plots)

                        if process_logs: process_log.add_line('correction parameters: %f (slope ratio), %f (intercept difference)'
                                                             % (adjust_obj.slope_ratio, adjust_obj.intercept_diff), log_level)
                        lin_model_params = adjust_obj.get_lin_model_params()


                        if process_logs:
                            process_log.add_line(
                                '--var_ratio_refdata: %f' % lin_model_params['VarRatio_timeframes_refdata(before/after)'],
                                log_level)
                            process_log.add_line(
                                '--var_ratio_testdata: %f' % lin_model_params['VarRatio_timeframes_testdata(before/after)'],
                                log_level)
                            process_log.add_line(
                                '--mean_diff_refdata: %f' % lin_model_params['MeanDiff_timeframes_refdata(before-after)'],
                                log_level)
                            process_log.add_line(
                                '--mean_diff_testdata: %f' % lin_model_params['MeanDiff_timeframes_testdata(before-after)'],
                                log_level)


                        if save_adjusted_data:  # break was found, attempt adjustment
                            if found_break or adjust_always:
                                if process_logs:
                                    corr = data_daily.corr().loc['testdata','refdata']
                                    process_log.add_line(
                                        'Correlation between refdata and testdata before adjustment: %f' % corr, log_level)
                                if adjust_all_values_before_break:
                                    data_to_adjust = pd.concat([df_time.loc[:breaktime,:],
                                                                data_daily.loc[breaktime + pd.DateOffset(1):,:]], axis=0)
                                else:
                                    data_to_adjust = data_daily
                                data_adjusted = adjust_obj.adjust_data(data_to_adjust.loc[:,['testdata']],
                                                                       return_part='both' if retries == 0 else 'first')

                                adj_status = 1  # '1: Adjusted Data for time frame saved'

                                #Replace data_daily testdata with adjusted testdata
                                data_daily['testdata'] = data_adjusted.loc[timeframe[0]:, 'adjusted']
                                #Save adjusted data to gpi time series
                                df_time.loc[data_adjusted.index, 'testdata'] = data_adjusted['adjusted']
                                #Save the bias corrected reference data as well
                                df_time.loc[data_daily.index, 'refdata'] = data_daily['refdata']

                                corr = data_daily.corr().loc['testdata', 'refdata']
                                if process_logs: process_log.add_line(
                                    'Correlation between refdata and testdata after adjustment: %f' % corr, log_level)

                            else:
                                adj_status = 8  # '8: Adjusted Data is not being stored'

                            # Some statistical parameters that should be fit
                            data_daily, stats = timeseries_stats(data_daily, breaktime, timeframe, True)
                        else:
                            adj_status = 8  # '8: Adjusted Data is not being stored'


                    except Exception as e:
                        print '%s: Adjustment failed: %s' % (str(breaktime.date()), e)
                        adj_status = int(str(e)[0])
                        lin_model_params = {'adjustment_status': adj_status}
                        adjust_obj = None
                        if process_logs: process_log.add_line('ERROR: %s' % str(e), log_level)
                    finally:
                        lin_model_params.update({'adjustment_status': adj_status})

                        if retries == 0:
                            test_results_before_adjustment.update(test_result)
                            lin_model_params_before_adjustment.update(lin_model_params)
                        else:
                            test_results_after_adjustment.update(test_result)
                            lin_model_params_after_adjustment.update(lin_model_params)

                    if not found_break:  # no break was found
                        if retries == 0:
                            gpi_test_results['homogeneous_original'].append(breaktime)
                            if process_logs: process_log.add_line('No break detected', log_level)
                            print '%s: No break detected' % str(breaktime.date())
                            break
                        else:  # break was removed after some iteration
                            df_time['testdata_original'] = df_time['testdata']
                            gpi_test_results['homogeneous_adjusted'].append(breaktime)
                            if process_logs: process_log.add_line('Break removed after %i retries' % retries, log_level)
                            print '%s: Break removed after %i iteration(s)' % (str(breaktime.date()), retries)
                            break
                    elif adjust_obj:
                        retries += 1
                    else:
                        break

                if model_plots:
                    if adjust_obj and retries > 0:  # TODO add this and retries > 0: # show plots only for points that were adjusted
                        adjust_obj.save_plots(model_plots_dir, '%i_%s' % (gpi, str(breaktime.date())))
                    else:
                        plt.close('all')


                if 'test_results' in test_results_after_adjustment:
                    _, _, found_break = test_obj.check_testresult(test_results_after_adjustment)
                    if found_break:
                        if process_logs: process_log.add_line('Could not remove %s break after %i retries'
                                                              % (str(break_found_by), retries), 2)
                        print ('%s: Could not remove %s break after %i retries :('
                               % (str(breaktime.date()), str(break_found_by), retries))
                        bad_gpis.append(gpi)
                        gpi_test_results['inhomogeneous_adjusted'].append(breaktime)

                        lin_model_params_after_adjustment['adjustment_status'] = 4  # '4: max number of iterations reached'

                        # TODO: Test this, include this?
                        # If there is still a break replace values by original ones
                        if process_logs: process_log.add_line('save UNADJUSTED values as break was not removed', 1)
                        # TODO: This inserts false valus for adjusting all values before break
                        df_time.loc[data_adjusted.index, 'testdata'] = df_time.loc[
                            data_adjusted.index, 'testdata_original']
                        if priority:
                            '''
                            # Check if the prioritized break was removed and save if it was, else replace with original values
                            if test_obj.compare_testresults(test_results_before_adjustment,
                                                            test_results_after_adjustment,
                                                            priority=priority):
                                # If adjustment did not improve results
                                df_time[timeframe[0]:timeframe[1]] = original_values[timeframe[0]:timeframe[1]]
                                lin_model_params_after_adjustment[
                                    'adjustment_status'] = '5: max. iter. reached w.o improvements'
                            '''
                            raise NotImplementedError
                        else:
                            pass
                    else:
                        df_time.loc[data_adjusted.index, 'testdata_original'] =  df_time.loc[data_adjusted.index, 'testdata']


                save_objects['test_results_before_adjustment'].add_data(test_results_before_adjustment,
                                                                        gpi, time=breaktime)
                save_objects['lin_model_params_before_adjustment'].add_data(lin_model_params_before_adjustment,
                                                                            gpi, time=breaktime)
                if save_adjusted_data:
                    save_objects['lin_model_params_after_adjustment'].add_data(lin_model_params_after_adjustment,
                                                                               gpi, time=breaktime)
                    save_objects['test_results_after_adjustment'].add_data(test_results_after_adjustment,
                                                                           gpi, time=breaktime)

                if ts_plots[0] and 'refdata' in data_daily.columns:
                    cols = 3
                    ax = ts_frames_fig.add_subplot(cols, int(np.ceil(float(breaktimes.size) / float(cols))), i+1)
                    plot_df = data_daily[['testdata', 'testdata_original', 'refdata']].rename(
                        columns={'testdata':'testdata_after_adjustment', 'testdata_original': 'testdata_before_adjustment'})
                    adjustment_ts_plot(plot_df, [breaktime], gpi, ref_prod, gpi_test_results,
                                       title= 'Adjustment %s (GPI:%i)' % (str(breaktime.date()), gpi),
                                       resample=ts_plots[1], ax = ax, legend=False)



            if save_adjusted_data:
                dataset.write(gpi, df_time.rename(columns={'testdata': test_obj.test_prod + '_ADJUSTED'}))
            if ts_plots[0]:
                save_path = os.path.join(workfolder, 'ts_plots')
                plot_df = df_time[['testdata_original', 'refdata_original', 'testdata_RAW']].rename(
                    columns={'testdata_original': 'testdata_after_adjustment',
                             'testdata_RAW': 'testdata_before_adjustment',
                             'refdata_original': 'refdata'})

                adjustment_ts_plot(plot_df, breaktimes, gpi, ref_prod, gpi_test_results, resample=ts_plots[1],
                                   title='CCI SM Adjustment (GPI:%i)' % gpi, save_path=save_path, legend=True)
                ts_frames_fig.subplots_adjust(hspace=0.5)
                ts_frames_fig.savefig(os.path.join(save_path, '%i_adjusted_frames' % gpi))


        # Add Info to log file
        log_file.add_line('%s: Finished testing for cell %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), cell))
        log_file.add_line('Bad GPIS for cell %s: %s' % (cell, str(bad_gpis)))

    q.put(True)


def start(test_prod, ref_prod, path, cells_identifier='global', skip_times=None, anomaly=None, save_adjusted_data=None,
          parallel_processes=8):
    # type: (str, str, str, Union[list,str], Union[list,None], Union[str,None], bool, int) -> None
    '''

    :param test_prod: cci_version_product
    :param ref_prod: merra2
    :param path:
    :param cells_identifier: global, Australia, North_America or [cellslist]
    :param skip_times: [int] eg. [0,1,2] to skip certain breaktimes
    :param anomaly: None or 'ccirange' or 'timeframe'
    :param perform_adjustment: path to where adjusted data is saved
    :param parallel_processes: number of parallel processes
    :return:
    :return:
    '''
    if skip_times:
        warnings.warn('Ignoring Breaktimes with activated adjustment!')

    workfolder = create_workfolder(path)

    log_file = LogFile(workfolder, 'log.txt',
                       {'Test Product': test_prod,
                        'Reference Product': ref_prod,
                        'Anomaly Data': anomaly,
                        'Processed Cells': cells_identifier})

    times_obj = CCITimes(test_prod, ignore_position=True,
                         skip_breaktimes=skip_times)  # TODO: activate for location conditional Timeframes

    results = ['test_results_before_adjustment', 'lin_model_params_before_adjustment']
    if save_adjusted_data:
        results = results + ['test_results_after_adjustment', 'lin_model_params_after_adjustment']

    save_objects = dict.fromkeys(results)

    for name in save_objects.keys():
        save_objects[name] = RegularGriddedCellData(grid=SMECV_Grid_v042(),
                                                    path=os.path.join(workfolder, 'temp_homogtest_cellfiles', name),
                                                    times=times_obj.get_times(None, as_datetime=True)['breaktimes'],
                                                    resolution=(0.25, 0.25))

    cells = list(
        split(split_cells(cells_identifier, SMECV_Grid_v042()), parallel_processes))  # split cells for processes

    processes = []
    q = Queue()
    finished_processes = []

    for process_no in range(parallel_processes):
        cells_for_process = cells[process_no]
        p = Process(target=process_for_cells, args=(q,
                                                    workfolder,
                                                    skip_times,
                                                    times_obj,
                                                    save_objects,
                                                    cells_for_process,
                                                    log_file,
                                                    test_prod,
                                                    ref_prod,
                                                    anomaly,
                                                    process_no,
                                                    save_adjusted_data))
        processes.append(p)
        p.start()

    for process in processes:
        finished_processes.append(q.get(True))

    for process in processes:
        process.join()

    # _, saved_gpis = join_files(workfolder, csv_files)
    print('Finished Testing (and Adjustment)')
    log_file.add_line('=====================================')
    # Global files and images from testing

    global_images_folder_name = 'global_images_plots'

    for name, save_obj in save_objects.iteritems():
        var_meta = {'test_results': get_test_meta()[0],
                    'test_status': get_test_meta()[1],
                    'adjustment_status': get_adjustment_meta()}
        global_meta = {'test_prod': test_prod, 'ref_prod': ref_prod}

        filename = 'GLOBAL_' + name + '.nc'
        if 'before_adjustment' in filename:
            plotfile_name = '%s_before_adjustment'
        else:
            plotfile_name = '%s_after_adjustment'

        # TODO: Parallelise this
        save_obj.make_global_file(filepath=workfolder,
                                  filename=filename,
                                  fill_nan=False,
                                  mfdataset=False,
                                  keep_cell_files=True,
                                  drop_variables=None,
                                  global_meta_dict=global_meta,
                                  var_meta_dicts=var_meta)  # Merge test result cell files to global file

        log_file.add_line('Merged cell files to global netcdf file: %s' % filename)

        # TODO: Put this in the plotting functions
        if not os.path.isdir(os.path.join(workfolder, global_images_folder_name)):
            os.mkdir(os.path.join(workfolder, global_images_folder_name))

        image_plots_path = os.path.join(workfolder, global_images_folder_name)

        if 'test_results' in filename:
            inhomo_plot_with_stats(os.path.join(workfolder, filename), image_plots_path,
                                   plotfile_name % 'TEST_RESULT')

            log_file.add_line('Nice Test Results Plots for %s as %s with meta data: %s' % (filename,
                                                                                           plotfile_name,
                                                                                           str(var_meta[
                                                                                                   'test_results'])))

            show_processed_gpis(os.path.join(workfolder, filename), image_plots_path,
                                'test_status', plotfile_name % 'TEST_STATUS')

            log_file.add_line('Nice Test Status Plots for %s as %s with meta data: %s' % (filename,
                                                                                          plotfile_name,
                                                                                          str(var_meta['test_status'])))
        elif 'lin_model_params' in filename:
            show_processed_gpis(os.path.join(workfolder, filename), image_plots_path,
                                'adjustment_status', plotfile_name % 'ADJUSTMENT_STATUS')

            log_file.add_line('Nice Adjustment Status Plots for %s as %s with meta data: %s' % (filename,
                                                                                                plotfile_name,
                                                                                                str(var_meta[
                                                                                                        'adjustment_status'])))

    return


if __name__ == '__main__':
    # Refproduct must be one of gldas-merged,gldas-merged-from-file,merra2,ISMN-merge
    # Testproduct of form cci_*version*_*product*

    start('CCI_41_COMBINED',
          'merra2',
          r'D:\users\wpreimes\datasets\HomogeneityTesting_data',
          # ''/data-write/USERS/wpreimes/HomogeneityTesting_data',
          cells_identifier=[2173],  # 'global', a continent or a list of cells
          skip_times=[1,4],  # list of breaktimes to skip
          anomaly=False,  # False, timeframe or ccirange
          save_adjusted_data=True,
          parallel_processes=1)
