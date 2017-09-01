# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 17:05:55 2017

@author: wpreimes
"""
from typing import Union
import os
from datetime import datetime
from multiprocessing import Process, Queue
from pynetcf.time_series import GriddedNcIndexedRaggedTs
from otherfunctions import resample_to_monthly, csv_read_write, join_files, split, create_workfolder
from interface import HomogTest
from cci_timeframes import CCITimes
from save_data import Results2D, extract_adjust_infos, extract_test_infos, save_Log, GlobalResults
from grid_functions import cells_for_continent
from adjustment import adjustment, adjustment_params
from otherplots import show_processed_gpis, inhomo_plot_with_stats, longest_homog_period_plots
from smecv_grid.grid import SMECV_Grid_v042

def process_for_cells(q, workfolder, save_obj, times_obj, cells, log_file, test_prod, ref_prod, anomaly, backward_extended_timeframes,
                      process_no, perform_adjustment=True):

    # Function for multicore processing

    print('Start process %i' % process_no)
    process_csv_file = 'saved_points_%s.csv' % process_no
    ########################################## TODO: Values to change
    adjusted_data_path = os.path.join(workfolder,test_prod+'_adjusted')
    get_all_B = False # TODO: Try this
    min_monthly_values = 10 # If CCI data contains less values, the monthly value is resampled as nan
    min_data_size = 3 # Minimum number of monthly values before/after break to perform homogeneity testing # TODO: different value?
    ##########################################
    tests = ['wilkoxon', 'fligner_killeen']

    test_obj = HomogTest(test_prod,
                         ref_prod,
                         tests,
                         0.01,  # TODO: Choose higher alpha value eg 0.05 or 0.1
                         anomaly,
                         perform_adjustment)

    if perform_adjustment:
        dataset = GriddedNcIndexedRaggedTs(path=adjusted_data_path, grid=(save_obj.pre_process_grid), mode='w')
        max_retries = 5
    else:
        max_retries = 1
        dataset = None

    for icell, cell in enumerate(cells):

        print 'Processing QDEG Cell %i (iteration %i of %i)' % (cell, icell + 1, len(cells))
        log_file.add_line('%s: Start Testing Cell %i' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), cell))
        grid_points = (save_obj.pre_process_grid).grid_points_for_cell(cell)[0]

        for iteration, gpi in enumerate(grid_points):
            #if gpi not in [369909, 365585, 375665, 320970, 322412]: continue

            if iteration%10 == 0:
                print 'processing gpi %i (%i of %i)' % (gpi, iteration, grid_points.size)
            '''
            if test_obj.ref_prod == 'ISMN-merge':
                valid_insitu_gpis = test_obj.ismndata.gpis_with_netsta

                if gpi not in valid_insitu_gpis.keys():
                    continue
            '''
            # test_obj.save_as_mat(gpi=gpi)
            try:
                df_time = test_obj.read_gpi(gpi, start=test_obj.range[0], end=test_obj.range[1])
            except:
                continue

            times = times_obj.get_times(gpi, as_datetime=True)
            timeframes = times['timeframes']
            breaktimes = times['breaktimes']

            homogeneous_breaktimes = []
            for i, (timeframe, breaktime) in enumerate(zip(timeframes, breaktimes)):
                retries = 0
                adjresult = False
                testresult = False
                # In case that after adjustment a break is still found, re-correct the reference data to
                # the (insufficiently) corrected test data and repeat adjustment (< max_retries times)
                while retries < max_retries and (not testresult or
                                                     not test_obj.check_testresult(testresult)[1]):

                    if backward_extended_timeframes: # Cut data to timeframe
                        while timeframe[1] in homogeneous_breaktimes:
                            timeframe[1] = times_obj.get_adjacent(gpi, timeframe, -1)[1] #Use previous timeframe

                    data_daily = df_time[timeframe[0]:timeframe[1]] # keep this for adjustment #TODO: Allow Extended Timeframes?
                    try:  # Homogeneity Testing
                        data_daily = data_daily.dropna()
                        # Correct bias in reference data
                        data_daily['refdata'] = test_obj.ref_data_correction(data_daily[['testdata', 'refdata']])
                        # resample to monthly value
                        data = resample_to_monthly(data_daily, min_monthly_values).dropna()
                        # Checks if there is data left and calculates spearman correlation
                        corr, pval = test_obj.check_corr(data)
                        data, _, _= test_obj.group_by_breaktime(data,
                                                                breaktime,
                                                                min_data_size=min_data_size,
                                                                ignore_exception=False)

                        # Calculate difference TimeSeries
                        data['Q'] = data['testdata'] - data['refdata']
                        # Run Tests
                        _, testresult = test_obj.run_tests(data=data)
                        testresult.update({'test_status': '0: Testing successful'})
                    except Exception as e:
                        testresult = {'test_status': str(e)}
                        break

                    break_found_by, check_result = test_obj.check_testresult(testresult)
                    if check_result == True: # no break was found
                        homogeneous_breaktimes.append(breaktime)
                        if retries == 0:
                            #print '%s: Adjustment not necessary' % str(breaktime.date())
                            break  # No breaks found, escape loop
                        else:  # break was removed, escape loop
                            print '%s: Break removed after %i retries' % (str(breaktime.date()), retries)
                            break
                    elif get_all_B or perform_adjustment:  # break was found, attempt adjustment
                        '''
                        # Adjustment 
                        1) Perform adjustment on the dataframe
                            ---Use equal amouunt of data after the break than before and vice versa (never more than to next/previous break point +margin
                            ----Use maximal amount of adjusted data (from breaktime until end of dataset)
                            -- Adjust model param 1, model param 2 or both
                        2) replacce data in df_time with adjusted data
                        3) Iterate over all breaktimes
                        4) Save df time to file
                        5) Iterate over all gpis
                        '''
                        try:  # Time Series Adjustment
                            B, corr = adjustment_params(data=data,   #todo: or monthly??
                                                        breaktime=breaktime)

                            #TODO: point dependent adjustment parameter?
                            '''
                            if break_found_by == ['wilkoxon']:
                                adjust_param = 'add'
                            elif break_found_by == ['fligner_killeen']:
                                adjust_param = 'mult'
                            else:
                                adjust_param = 'both'
                            '''
                            adjusted, adjresult = adjustment(data_daily=data_daily,
                                                             B=B,
                                                             breaktime=breaktime,
                                                             adjust_param='both',
                                                             adjust_part='first',
                                                             return_part='all' if i == 0 else 'first',
                                                             check_adjustment=True)


                            adjresult.update({'adj_status': '0: Adjustment performed'})
                            df_time.loc[adjusted.index, 'testdata'] = adjusted
                            # Replace df_time (all values) in timeframe with adjusted values
                            retries += 1
                        except Exception as e:
                            adjresult = {'adj_status': str(e)}
                            print '%i: Adjustment failed due to exception: %s' % (gpi, e)
                            break
                    else:
                        # There is a break, but adjustment is not selected.
                        break

                # Merge test results and adjustresults for gpi and add to save object
                if testresult['test_status'][0] == '0' and not test_obj.check_testresult(testresult)[1]: #todo: ugly
                    print ('%s: Could not remove break after %i retries :(' % (str(breaktime.date()), max_retries))
                testresult = extract_test_infos(testresult)
                if not adjresult:
                    adjresult = {'adj_status': '9: Not adjusted'}
                adjresult = extract_adjust_infos(adjresult)
                merged_results = testresult.copy()
                merged_results.update(adjresult)
                save_obj.add_data(gpi, str(breaktime.date()), merged_results)

            if dataset:
                dataset.write(gpi, df_time[['testdata']].rename(columns={'testdata': test_obj.test_prod + '_adjusted'}))

        saved_points = save_obj.save_to_gridded_netcdf()  # Save data for the cell to netcdf file

        csv_read_write(os.path.join(workfolder, process_csv_file), 'write', saved_points)

        # Add Info to log file
        log_file.add_line('%s: Finished testing for cell %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), cell))

    q.put(process_csv_file)


def start(test_prod, ref_prod, path, cells='global', skip_times=None, anomaly=None, perform_adjustment=None,
          backward_extended_timeframes=True, parallel_processes=8):
    # type: (str, str, str, Union[list,str], Union[list,None], str, Union[str,None], bool, int) -> None
    '''

    :param test_prod: cci_version_product
    :param ref_prod: merra2
    :param path:
    :param cells: global, Australia, North_America or [cellslist]
    :param skip_times: [int] eg. [0,1,2] to skip certain breaktimes
    :param anomaly: None or 'ccirange' or 'timeframe'
    :param perform_adjustment: path to where adjusted data is saved
    :param backward_extended_timeframes: expand timeframe over tested breaktimes, if no break was found
    :param parallel_processes: number of parallel processes
    :return:
    '''
    if skip_times and perform_adjustment:
        raise Warning('Ignoring Breaktimes with activated adjustment!')

    times_obj = CCITimes(test_prod,ignore_position=True, skip_times=skip_times) #TODO: activate for location conditional Timeframes

    grid = SMECV_Grid_v042()

    workfolder = create_workfolder(path)

    log_file = save_Log(workfolder, test_prod, ref_prod, anomaly, cells)

    if cells == 'global':
        cells = grid.get_cells()
    elif isinstance(cells,str):
        cells = cells_for_continent(cells)

    save_obj = Results2D(grid, workfolder, times_obj.get_times(None, as_datetime=False)['breaktimes'],
                         cell_files_path=os.path.join(workfolder,'temp_homogtest_cellfiles' ))

    processes = []
    q = Queue()
    csv_files = []
    # Split cells in equally sized packs for multiprocessing
    for process_no, cell_pack in enumerate(list(split(cells, parallel_processes))):
        p = Process(target=process_for_cells, args=(q, workfolder, save_obj, times_obj, cell_pack,
                                                    log_file, test_prod, ref_prod, anomaly, backward_extended_timeframes,
                                                    process_no, perform_adjustment))
        processes.append(p)
        p.start()

    for process in processes:
        csv_files.append(q.get(True))

    for process in processes:
        process.join()

    _, saved_gpis = join_files(workfolder, csv_files)
    print('Finished Testing (and Adjustment)')
    log_file.add_line('=====================================')
    post_process_grid = save_obj.save_subgrid(saved_gpis)  # Save test gpis subset to gridfile

    # Global files and images from testing
    save_obj_glob = GlobalResults(save_obj, post_process_grid)
    global_file_name = save_obj_glob.save_global_file(keep_cell_files=True)  # Merge test result cell files to global file
    log_file.add_line('Merged files to global file: %s' % global_file_name)
    image_files_names = save_obj_glob.create_image_files()  # Create 2D image files from test results
    log_file.add_line('Create global image files:')
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
    log_file.add_line('=====================================')
    if perform_adjustment:
        log_file.add_line('%s: Start TS Adjustment' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        pass
        #Do TS adjutement and save resutls to path
    '''


if __name__ == '__main__':
    # Refproduct must be one of gldas-merged,gldas-merged-from-file,merra2,ISMN-merge
    # Testproduct of form cci_*version*_*product*
    start('CCI_41_COMBINED',
          'merra2',
          r'D:\users\wpreimes\datasets\HomogeneityTesting_data\output',
          cells=[777,742],
          skip_times=None,
          anomaly=False,
          backward_extended_timeframes = True,
          perform_adjustment=True,
          parallel_processes=2)
