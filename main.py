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
from otherfunctions import csv_read_write, join_files, split, create_workfolder
from interface import BreakTestBase
from cci_timeframes import CCITimes
from save_data import RegularGriddedCellData, LogFile
from grid_functions import cells_for_continent
from adjustment import LinearAdjustment
from otherplots import show_processed_gpis, inhomo_plot_with_stats, longest_homog_period_plots
from smecv_grid.grid import SMECV_Grid_v042
import warnings
import pandas as pd


def split_cells(cells_identifier, grid):
    if cells_identifier == 'global':
        cells = grid.get_cells()
    elif isinstance(cells_identifier, str):
        cells = cells_for_continent(cells_identifier)
    else:
        cells = cells_identifier
    return cells

def process_for_cells(q, workfolder, skip_times, times_obj, save_objects, cells, log_file, test_prod, ref_prod, anomaly,
                      process_no, save_adjusted_data=True):

    # Function for multicore processing

    print('Start process %i' % process_no)

    ########################################## TODO: Values to change
    adjusted_data_path = os.path.join(workfolder, test_prod + '_ADJUSTED')
    min_data_for_temp_resampling = 0.33  # If CCI data contains less values, the monthly value is resampled as nan
    min_data_size = 3  # Minimum number of monthly values before/after break to perform homogeneity testing # TODO: different value?
    resample_method = 'M'  # or D
    refdata_correction_for = 'timeframe'  # iteration, full or timeframe
    tests = {'mean': 'wilkoxon', 'var': 'scipy_fligner_killeen'}
    backward_extended_timeframes = False
    model_plots = True
    quantile_filtering_increase = 0.01  # Increase of quantile per iteration (0 ,0.03, 0.04,... and 1, 0.97, 0.94,....]
    first_iter_mode = 'both'  # 'd' (BAD!!!) or 'both' for adjustment if first iteration found wk break
    max_retries = 10 # max number of iterations if break is still found after adjustment
    adjust_always = False
    priority = ['mean','var']
    ##########################################

    if process_no == 0:
        log_file.add_line('-------------------------------')
        log_file.add_line('%s: %s' % ('adjusted data saved to', adjusted_data_path))
        log_file.add_line('%s: %s' % ('min percent. of valid data for resampling', min_data_for_temp_resampling))
        log_file.add_line('%s: %s' % ('min number of point pairs for modelling', min_data_size))
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
        log_file.add_line('-------------------------------')

    test_obj = BreakTestBase(test_prod,
                             ref_prod,
                             tests,
                             0.01,  # TODO: Choose higher alpha value eg 0.05 or 0.1
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


    for icell, cell in enumerate(cells):

        print 'Processing QDEG Cell %i (iteration %i of %i)' % (cell, icell + 1, len(cells))
        log_file.add_line('%s: Start Testing Cell %i' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), cell))

        for iteration, gpi in enumerate(grid.grid_points_for_cell(cell)[0]):
            print 'processing gpi %i (%i of %i)' % (gpi, iteration, grid.grid_points_for_cell(cell)[0].size)

            homogeneous_breaktimes = []  # List of tested and approved break times for extending time frames

            try: # load testdata and reference data
                df_time = test_obj.read_gpi(gpi,
                                            start=test_obj.range[0],
                                            end=test_obj.range[1])
                # TODO: Do refdata correction for every iteration?
                corrected_refdata = test_obj.ref_data_correction(df_time[['testdata', 'refdata']])
                df_time['refdata'] = corrected_refdata
            except:
                continue

            original_values = df_time.copy()

            times = times_obj.get_times(gpi, as_datetime=True)
            timeframes, breaktimes = times['timeframes'], times['breaktimes']

            for i, (timeframe, breaktime) in enumerate(zip(timeframes, breaktimes)):

                lin_model_params_before_adjustment = {}
                lin_model_params_after_adjustment = {}
                test_results_before_adjustment = {}
                test_results_after_adjustment = {}

                retries = 0
                # In case that after adjustment a break is still found, re-correct the reference data to
                # the (insufficiently) corrected test data and repeat adjustment (< max_retries times)
                while (retries <= min_retries) or (retries < max_retries):

                    if backward_extended_timeframes:  # Cut data to timeframe
                        while timeframe[1] in homogeneous_breaktimes:  # Use end of previous time frame
                            timeframe[1] = times_obj.get_adjacent(gpi, timeframe, -1)[1]

                    data_daily = df_time[timeframe[0]:timeframe[1]].copy()  # change with adjustment iterations

                    try:  # Homogeneity Testing
                        data_daily = data_daily.dropna(subset=['testdata', 'refdata'])

                        data_daily['Q'] = data_daily['testdata'] - data_daily['refdata'] # Calculate difference TimeSeries

                        data_resampled = test_obj.temp_resample(data_daily,
                                                                'M',
                                                                min_data_for_temp_resampling)


                        # Checks if there is data left and calculates spearman correlation
                        corr, pval = test_obj.check_corr(data_resampled)
                        data_resampled, _, _ = test_obj.group_by_breaktime(data_resampled,
                                                                           breaktime,
                                                                           min_data_size=min_data_size,
                                                                           ignore_exception=False)

                        # Correct bias in reference data, again (this time for monthly values)
                        #TODO: Do this again, changes the model of part 2 for eacht iteration?
                        #if refdata_correction_for == 'iteration' or retries == 0:
                        #    data_resampled.loc[:,'refdata'] = test_obj.ref_data_correction(data_resampled.loc[:,('testdata', 'refdata')])
                        #    data_resampled.loc[:, 'Q'] = data_resampled.loc[:, 'testdata'] - data_resampled.loc[:, 'refdata']

                        # Run Tests
                        test_result = test_obj.run_tests(data=data_resampled)

                        # TODO: Check this for a timeframe with a break
                        failed_tests, break_found_by, found_break = test_obj.check_testresult(test_result)

                    except Exception as e:
                        # Testing Failed
                        adjust_obj = None
                        print '%s: Testing failed: %s' % (str(breaktime.date()), e)
                        test_result = {'test_status': int(str(e)[0])}
                        test_results_before_adjustment.update(test_result)
                        break

                    try:
                        #TODO: Filter only for first iteration, for next iterations use residuals to remove bad values
                        #if retries ==0:
                        lq = 0.05
                        uq = 0.95
                        data_daily_filtered = test_obj.quantile_filtering(data_daily,# if retries == 0 else data_daily_filtered,
                                                                          breaktime,
                                                                          'both', #if retries==0 else 'first',
                                                                          lq if retries == 0 else lq + retries * quantile_filtering_increase,
                                                                          uq if retries == 0 else uq - retries * quantile_filtering_increase)

                        adjust_obj = LinearAdjustment(data_daily_filtered, #TODO: DataDaily, DataDailyFiltered oder MonthlyData
                                                      breaktime,
                                                      'first',
                                                      'both',
                                                      (retries, adjust_obj.model_plots) if (retries > 0 and adjust_obj) else model_plots)

                        lin_model_params = adjust_obj.get_lin_model_params()

                        if save_adjusted_data:  # break was found, attempt adjustment
                            if found_break or adjust_always:
                                data_adjusted = adjust_obj.adjust_data(data_daily, 'both' if retries == 0 else 'first')
                                adj_status = 0 #'0: Adjusted Data for time frame saved'

                                df_time.loc[data_adjusted.index, 'testdata'] = data_adjusted
                            else:
                                adj_status = 8 #'8: Adjusted Data is not being stored'
                        else:
                            adj_status = 8 #'8: Adjusted Data is not being stored'
                    except Exception as e:
                        print '%s: Adjustment failed: %s' % (str(breaktime.date()), e)
                        adj_status = int(str(e)[0])
                        lin_model_params = {'adjustment_status': adj_status}
                        adjust_obj = None
                    finally:
                        lin_model_params.update({'adjustment_status': adj_status})

                        if retries == 0:
                            test_results_before_adjustment.update(test_result)
                            lin_model_params_before_adjustment.update(lin_model_params)
                        else:
                            test_results_after_adjustment.update(test_result)
                            lin_model_params_after_adjustment.update(lin_model_params)

                    if not found_break:  # no break was found
                        homogeneous_breaktimes.append(breaktime)
                        if retries == 0:
                            print '%s: No break detected' % str(breaktime.date())
                            break
                        else:  # break was removed after some iteration
                            print '%s: Break removed after %i iteration(s)' % (str(breaktime.date()), retries)
                            break
                    elif adjust_obj:
                        retries += 1
                    else:
                        break


                if model_plots:
                    if adjust_obj and retries > 0:# TODO add this and retries > 0: # show plots only for points that were adjusted
                        adjust_obj.save_plots(model_plots_dir, '%i_%s' %(gpi, str(breaktime.date())))
                    else:
                        plt.close('all')

                if test_results_after_adjustment:
                    _, _, found_break = test_obj.check_testresult(test_results_after_adjustment)
                    if found_break:
                        print ('%s: Could not remove break after %i retries :(' % (str(breaktime.date()), max_retries))
                        lin_model_params_after_adjustment['adjustment_status'] = '4: max number of iterations reached'

                    #TODO: Test this
                    if test_obj.compare_testresults(test_results_before_adjustment,
                                                    test_results_after_adjustment,
                                                    priority=priority):
                        # If adjustment did not improve results
                        df_time[timeframe[0]:timeframe[1]] = original_values[timeframe[0]:timeframe[1]]
                        lin_model_params_after_adjustment['adjustment_status'] = '5: max. iter. reached w.o improvements'

                #elif test_obj.check_testresult(test_results_before_adjustment)

                '''
                if 'test_status' in test_results_before_adjustment.keys():
                    test_results_before_adjustment['test_status'] = np.nan
                if 'test_status' in test_results_after_adjustment.keys():
                    test_results_after_adjustment['test_status'] = np.nan
                if 'adjustment_status' in lin_model_params_after_adjustment.keys():
                    lin_model_params_after_adjustment['adjustment_status'] = np.nan
                if 'adjustment_status' in lin_model_params_before_adjustment.keys():
                    lin_model_params_before_adjustment['adjustment_status'] = np.nan
                '''

                save_objects['test_results_before_adjustment'].add_data(test_results_before_adjustment,
                                                                        gpi, time=breaktime)
                save_objects['test_results_after_adjustment'].add_data(test_results_after_adjustment,
                                                                       gpi, time=breaktime)
                save_objects['lin_model_params_before_adjustment'].add_data(lin_model_params_before_adjustment,
                                                                            gpi, time=breaktime)
                save_objects['lin_model_params_after_adjustment'].add_data(lin_model_params_after_adjustment,
                                                                           gpi, time=breaktime)

            if save_adjusted_data:
                dataset.write(gpi, df_time[['testdata']].rename(columns={'testdata': test_obj.test_prod + '_ADJUSTED'}))
        # Add Info to log file
        log_file.add_line('%s: Finished testing for cell %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), cell))

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
    '''
    if skip_times:
        warnings.warn('Ignoring Breaktimes with activated adjustment!')

    workfolder = create_workfolder(path)

    log_file = LogFile(workfolder, {'Test Product': test_prod,
                                    'Reference Product': ref_prod,
                                    'Anomaly Data': anomaly,
                                    'Processed Cells': cells_identifier})

    times_obj = CCITimes(test_prod, ignore_position=True,
                         skip_times=skip_times)  # TODO: activate for location conditional Timeframes

    results = ['test_results_before_adjustment', 'lin_model_params_before_adjustment']
    if save_adjusted_data:
        results = results + ['test_results_after_adjustment', 'lin_model_params_after_adjustment']

    save_objects = dict.fromkeys(results)

    for name in save_objects.keys():
        save_objects[name] = RegularGriddedCellData(grid=SMECV_Grid_v042(),
                                                    path=os.path.join(workfolder, 'temp_homogtest_cellfiles', name),
                                                    times=times_obj.get_times(None, as_datetime=True)['breaktimes'],
                                                    resolution=(0.25,0.25))


    cells =  list(split(split_cells(cells_identifier, SMECV_Grid_v042()), parallel_processes)) # split cells for processes

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

    #_, saved_gpis = join_files(workfolder, csv_files)
    print('Finished Testing (and Adjustment)')
    log_file.add_line('=====================================')
    # Global files and images from testing


    for name, save_obj in save_objects.iteritems():
        #TODO: Parallelise this
        filename = 'GLOBAL_' + name + '.nc'
        global_file_name = save_obj.make_global_file(workfolder,filename, False, False, keep_cell_files=True)  # Merge test result cell files to global file

    #post_process_grid = save_obj.save_subgrid(saved_gpis,
    #                                          'breaktest_grid.nc',
    #                                          'Breaktest')  # Save test gpis subset to gridfile

    log_file.add_line('Merged files to global file: %s' % global_file_name)
    '''
    for i, image_file in enumerate(image_files_names, start=1):  # Create spatial plots from test results and coverage
        log_file.add_line('  NC Image File %i : %s' % (i, image_file))
        meta = show_processed_gpis(workfolder, 'test', image_file)
        log_file.add_line('  Test Results Plot %i : %s' % (i, str(meta)))
        stats = inhomo_plot_with_stats(workfolder, image_file)
        log_file.add_line('  Test Results Plot %i : %s' % (i, str(meta)))
    fn_long_per_plot = longest_homog_period_plots(workfolder)  # Create plot of longest homogeneous period
    log_file.add_line('  Plot of Longest Homogeneous Period : %s' % fn_long_per_plot)

    if perform_adjustment:
        for i, image_file in enumerate(image_files_names, start=1):
            log_file.add_line('  Adjusment coverage Plot %i : %s' % (i, str(stats)))
            meta = show_processed_gpis(workfolder, 'adjustment', image_file)
    '''

if __name__ == '__main__':
    # Refproduct must be one of gldas-merged,gldas-merged-from-file,merra2,ISMN-merge
    # Testproduct of form cci_*version*_*product*

    start('CCI_41_COMBINED',
          'merra2',
          r'D:\users\wpreimes\datasets\HomogeneityTesting_data',
          cells_identifier=[2244],  # 'global', a continent or a list of cells
          skip_times=None,  # list of breaktimes to skip
          anomaly=False,  # False, timeframe or ccirange
          save_adjusted_data=True,
          parallel_processes=1)
