# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 17:05:55 2017

@author: wpreimes
"""
from typing import Union
import os
from datetime import datetime
from pygeogrids.netcdf import load_grid
from multiprocessing import Process, Queue
from pynetcf.time_series import GriddedNcIndexedRaggedTs
import pandas as pd
from interface import HomogTest
from cci_timeframes import get_timeframes
from save_data import Results2D, save_Log, GlobalResults
from grid_functions import cells_for_continent
from adjustment import regression_adjustment
import numpy as np
from otherplots import show_tested_gpis, inhomo_plot_with_stats, longest_homog_period_plots
import csv
'''
Breaktimes sind Punkte an denen Inhomogenitäten vermutet werden
Laut Processing overview und Research letter:
    Combined Data: (ignore September 1987) , August 1991, January 1998, July 2002, Januar 2007, 
                    October 2011, July 2012, (additional May 2015)

used alpha in research letter:0.01
'''
def split(el, n):
    '''
    Split list of cells in n approx. equal parts for multiprocessing
    :param el: list of elements to split
    :param n: number of lists to split input up into
    :return: list
    '''
    k, m = divmod(len(el), n)
    return (el[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))

def create_workfolder(path):
    # type: (str) -> str
    i = 1
    while os.path.exists(os.path.join(path, 'v'+str(i))):
        i += 1
    else:
        os.makedirs(os.path.join(path, 'v'+str(i)))
    
    workfolder = os.path.join(path, 'v'+str(i))
    print('Create workfolder:%s' % str(workfolder))
        
    return workfolder

def join_files(filefolder, filelist):
    merged_file_name = 'saved_points.csv'
    merged = []

    for file in filelist:
        filepath = os.path.join(filefolder, file)
        data = csv_read_write(filepath, 'read')
        for row in data:
            merged.append(row)
        os.remove(filepath)

    merged_int = map(int, [item for sublist in merged for item in sublist])
    path = csv_read_write(os.path.join(filefolder, 'saved_points.csv'), 'write', merged_int)

    return path, merged_int


def csv_read_write(csv_path, mode, data=None):
    if mode == 'write':
        if not os.path.isfile(csv_path):
            with open(csv_path, 'wb') as file:
                wr = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
                wr.writerow(data)
        else:
            if os.path.isfile(csv_path):
                with open(csv_path, 'ab') as file:
                    wr = csv.writer(file,delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
                    wr.writerow(data)
        return csv_path
    if mode == 'read':
        return_data = []
        with open(csv_path, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                return_data.append(row)
        return return_data

def process_for_cells(q, workfolder, save_obj, cells, log_file, test_prod, ref_prod, anomaly, process_no, adjusted_ts_path=None):
    #Function for multicore processing

    print('Start process %i' %process_no)
    process_csv_file = 'saved_points_%s.csv' %process_no

    tests = ['wilkoxon', 'fligner_killeen']

    test_obj = HomogTest(test_prod,
                         ref_prod,
                         tests,
                         0.01,
                         anomaly,
                         adjusted_ts_path)
    if adjusted_ts_path:
        dataset = GriddedNcIndexedRaggedTs(path=adjusted_ts_path, grid=(save_obj.pre_process_grid), mode='w')

    for icell, cell in enumerate(cells):

        print 'Processing QDEG Cell %i (iteration %i of %i)' % (cell, icell+1, len(cells))

        log_file.add_line('%s: Start Testing Cell %i' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), cell))
        grid_points = (save_obj.pre_process_grid).grid_points_for_cell(cell)[0]


        for iteration, gpi in enumerate(grid_points):
            if iteration%10 == 0:
                print 'processing gpi %i of %i' %(iteration, grid_points.size)
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

            adjusted_dataframes = {str(breaktime.date()):[] for breaktime in test_obj.breaktimes}
            adjustment_results = {str(breaktime.date()):[] for breaktime in test_obj.breaktimes}

            for i, (timeframe, breaktime) in enumerate(zip(test_obj.timeframes, test_obj.breaktimes)):
                success = False
                if adjusted_ts_path:
                    max_retries = 3
                else:
                    max_retries=0
                retries = 0
                adjresult = {'slope': np.nan, 'intercept': np.nan, 'part1_B1': np.nan, 'part1_B2': np.nan,
                                     'part2_B1': np.nan, 'part2_B2': np.nan,
                                     'B1_aft_adjust': np.nan, 'B2_aft_adjust': np.nan,
                                     'status_adj': '5: Point not tested'}
                while not success:
                    # Cut data to timeframe
                    data = df_time[timeframe[0]:timeframe[1]]
                    try: # Homogeneity Testing
                        data = data.dropna()
                        # Checks if there is data left and calculates spearman correlation
                        corr, pval = test_obj.check_corr(data)
                        data = test_obj.group_by_breaktime(data, breaktime, min_data_size=3) #TODO: different value?
                        # Correct bias in reference data
                        data['bias_corr_refdata'] = test_obj.ref_data_correction(data[['testdata','refdata']])
                        # Calculate difference TimeSeries
                        data['Q'] = data['testdata'] - data['bias_corr_refdata']
                        # Run Tests
                        _, testresult = test_obj.run_tests(data=data)
                        testresult.update({'status_test': '0: Testing successful'})
                    except Exception as e:
                        testresult = {'status_test': str(e)}
                        success = True # Finish this breaktime
                        continue

                    if all(h == 0 for h in [testresult[test]['h'] for test in test_obj.tests]):
                        success = True # Finish this breaktime
                        #if retries == 0:
                        #    print '%i: No break found' % gpi
                        #else:
                        #    print '%i: Break removed' % gpi
                    elif retries == max_retries:
                        success = True # After max nbr of repetitions, give up
                        #print('%i: Max retries reached' % gpi)
                    else:
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
                        if adjusted_ts_path:
                            try: # Time Series Adjustment
                                df_adjusted, adjresult = regression_adjustment(data = data,
                                                                               breaktime = breaktime,
                                                                               adjust_param = 'both',
                                                                               adjust_part = 'first',
                                                                               return_part = 'all' if i==1 else 'first')

                                adjresult.update({'status_adj': '0: Adjustment performed'})

                                df_time.loc[df_adjusted.index]['testdata']= df_adjusted   # Replace df_time (all values) in timeframe with adjusted values
                                retries += 1
                            except Exception as e:
                                print str(e)
                                df_adjusted = data # is actually un-adjusted
                                adjresult = {'slope': np.nan, 'intercept': np.nan, 'part1_B1': np.nan, 'part1_B2': np.nan,
                                            'part2_B1': np.nan, 'part2_B2': np.nan,
                                            'B1_aft_adjust': np.nan, 'B2_aft_adjust': np.nan, 'status_adj': str(e)}

                                print '%i: Adjustment failed' %gpi
                                success = True # Not a success, but continue

                    if adjusted_ts_path and not adjresult: # If no break was found
                        adjresult = {'slope': np.nan, 'intercept': np.nan, 'part1_B1': np.nan, 'part1_B2': np.nan,
                                     'part2_B1': np.nan, 'part2_B2': np.nan,
                                     'B1_aft_adjust': np.nan, 'B2_aft_adjust': np.nan,
                                     'status_adj': '1: No adjustment necessary'}


                gpi_result = {'testresult' : testresult}
                if adjusted_ts_path:
                    gpi_result['adjresult'] = adjresult

                #Add Test results and adjustment results to data saving buffer
                save_obj.add_data(gpi, str(breaktime.date()), gpi_result)

            if adjusted_ts_path:
                dataset.write(gpi, df_time[['testdata']].rename(columns={'testdata':test_obj.test_prod + '_adjusted'}))


        saved_points = save_obj.save_to_gridded_netcdf() # Save data for the cell to netcdf file

        csv_read_write(os.path.join(workfolder, process_csv_file), 'write', saved_points)

        # Add Info to log file
        log_file.add_line('%s: Finished testing for cell %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), cell))

    q.put(process_csv_file)



def start(test_prod, ref_prod, path, cells='global',skip_times=None, anomaly=None, adjusted_ts_path=None,
          parallel_processes=8):
    # type: (str, str, str, Union[list,str], Union[list,None], str, Union[str,None]) -> None
    '''

    :param test_prod: cci_version_product
    :param ref_prod: merra2
    :param path:
    :param cells: global, Australia, North_America or [cellslist]
    :param skip_times: [int] eg. [0,1,2] to skip certain breaktimes
    :param anomaly: None or 'ccirange' or 'timeframe'
    :param adjusted_ts_path: path to where adjusted data is saved
    :return:
    '''
    if skip_times and adjusted_ts_path:
        raise Warning('Ignoring Breaktimes with activated adjustment!')

    testtimes = get_timeframes(test_prod, skip_times=skip_times)

    timeframes = testtimes['timeframes']
    breaktimes = testtimes['breaktimes']

    
    grid = load_grid(r"D:\users\wpreimes\datasets\grids\qdeg_land_grid.nc")

    workfolder = create_workfolder(path)
       
    log_file = save_Log(workfolder, test_prod, ref_prod, anomaly, cells)

    if isinstance(cells, str):
        if cells == 'global':
            cells = grid.get_cells()
        else:
            cells = cells_for_continent(cells)

    save_obj = Results2D(grid, workfolder, breaktimes)


    processes = []
    q = Queue()
    csv_files = []
    # Split cells in equally sized packs for multiprocessing
    for process_no, cell_pack in enumerate(list(split(cells, parallel_processes))):
        p = Process(target=process_for_cells, args=(q, workfolder, save_obj, cell_pack,
                                              log_file, test_prod, ref_prod, anomaly, process_no, adjusted_ts_path))
        processes.append(p)
        p.start()

    for process in processes:
        csv_files.append(q.get(True))

    for process in processes:
        process.join()

    _, saved_gpis = join_files(workfolder, csv_files)
    print('Finished Testing (and Adjustment)')
    log_file.add_line('=====================================')
    gridfile = save_obj.save_subgrid(saved_gpis)  # Save test gpis subset to gridfile
    log_file.add_line('Saved Grid of Tested Points: %s' % gridfile)

    # Global files and images from testing
    save_obj = GlobalResults(workfolder, 'gridded_files', breaktimes)
    global_file_name = save_obj.save_global_file(keep_cell_files=True) # Merge test result cell files to global file
    log_file.add_line('Merged files to global file: %s' % global_file_name)
    image_files_names = save_obj.create_image_files() # Create 2D image files from test results
    log_file.add_line('Create global image files:')
    for i, image_file in enumerate(image_files_names, start=1): # Create spatial plots from test results and coverage
        log_file.add_line('  NC Image File %i : %s' % (i, image_file))
        meta = show_tested_gpis(workfolder, image_file)
        log_file.add_line('  Test Results Plot %i : %s' % (i, str(meta)))
        stats = inhomo_plot_with_stats(workfolder, image_file)
        log_file.add_line('  Coverage Plot %i : %s' % (i, str(stats)))
    fn_long_per_plot = longest_homog_period_plots(workfolder)   # Create plot of longest homogeneous period
    log_file.add_line('  Plot of Longest Homogeneous Period : %s' % fn_long_per_plot)


    '''
    log_file.add_line('=====================================')
    if adjusted_ts_path:
        log_file.add_line('%s: Start TS Adjustment' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        pass
        #Do TS adjutement and save resutls to path
    '''

if __name__ == '__main__':
    # Refproduct must be one of gldas-merged,gldas-merged-from-file,merra2,ISMN-merge
    # Testproduct of form cci_*version*_*product*
    start('cci_31_combined',
          'merra2',
          r'H:\HomogeneityTesting_data\output',
          cells=[777],
          skip_times=None,
          anomaly=False,
          adjusted_ts_path=None,
          parallel_processes=1)
