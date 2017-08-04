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
import pandas as pd
from interface import HomogTest
from cci_timeframes import get_timeframes
from save_data import Results_2D, save_Log
from grid_functions import cells_for_continent
from adjustment import regression_adjustment
import numpy as np

#from otherplots import show_tested_gpis, inhomo_plot_with_stats
#from longest_homogeneous_period import calc_longest_homogeneous_period
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

def process_for_cells(q, workfolder, save_obj, cells, log_file, test_prod, ref_prod, anomaly, adjusted_ts_path=None):
    #Function for multicore processing


    test_obj = HomogTest(test_prod,
                         ref_prod,
                         0.01,
                         anomaly,
                         adjusted_ts_path)


    # Object for saving spatial information to netcdf file


    processed_cells = []

    for icell, cell in enumerate(cells):
        print 'Processing QDEG Cell %i (iteration %i of %i)' % (cell, icell, len(cells))

        log_file.add_line('%s: Start Testing Cell %i' % (str(datetime.now()), cell))
        grid_points = (save_obj.grid).grid_points_for_cell(cell)[0]

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
                    _, testresult = test_obj.run_tests(data=data,
                                                          tests=['wk', 'fk'])
                    testresult.update({'status': '0: Testing successful'})
                except Exception as e:
                    testresult = {'status': str(e)}

                #Add Test results to data saving buffer
                save_obj.add_data(gpi, str(breaktime.date()), testresult)

                '''
                if any(h == 1 for h in [testresult['wilkoxon']['h'], testresult['fligner_killeen']['h']]):
                    try: # Time Series Adjustment
                        df_adjusted, adjresult = regression_adjustment(data = data,
                                                                       breaktime = breaktime,
                                                                       adjust_part = 'first',
                                                                       return_part = 'all' if i==1 else 'first',
                                                                       test_adjusted = True)

                        adjresult.update({'status': '0: Adjustment performed'})
                    except Exception as e:
                        df_adjusted = data
                        adjresult = {'slope': np.nan, 'intercept': np.nan, 'part1_B1': np.nan, 'part1_B2': np.nan,
                                       'part2_B1': np.nan, 'part2_B2': np.nan}
                        adjresult = {'status': str(e)}
                else:
                    df_adjusted = data
                    adjresult = {'status': '1: No break detected'}

                adjusted_dataframes[str(breaktime.date())] = df_adjusted
                adjustment_results[str(breaktime.date())] = adjresult

                df_concat = pd.concat(adjusted_dataframes.values(), axis=1)







                if adjusted_ts_path:
                    adjust df_results

                # wenn testresults break finden adjustment über timeframe
                # Replace df_time[timeframe] with adjusted data
                '''



        # Add Info to log file
        save_obj.save_to_netcdf() # save rest to file
        log_file.add_line('%s: Finished testing for cell %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), cell))
        processed_cells.append(cell)
    if q:
        q.put(processed_cells)
    else:
        return processed_cells


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

    save_obj = Results_2D(workfolder, grid, breaktimes, buffer_size=100)

    if parallel_processes == 1:
        processed_cells = process_for_cells(None, workfolder, save_obj, cells, log_file, test_prod, ref_prod, anomaly)
    else:
        q = Queue()
        for cell_pack in list(split(cells, parallel_processes)): # Split cells in equally sized packs for multiprocessing
            p = Process(target=process_for_cells, args=(q, workfolder, save_obj, cell_pack,
                                                     log_file, test_prod, ref_prod, anomaly))
            p.start()

        processed_cells = []
        for i in range(len(breaktimes)):
            processed_cells.append(q.get(True))
    log_file.add_line('=====================================')
    log_file.add_line('Processed Cells: %s' %str(processed_cells))
    global_file_name = save_obj.save_global_file()
    log_file.add_line('Merged files to global file: %s' % global_file_name)
    image_files_names = save_obj.create_image_files()
    log_file.add_line('Created global image files: %s' % image_files_names)
    '''
    # Create nc file with longest period data
    ncfilename = calc_longest_homogeneous_period(workfolder, create_netcdf=True)
    log_file.add_line('%s: Saved longest homogeneous Period, startdate, enddate to %s' % (datetime.now().strftime('%Y-%m-%d_%H:%M:%S'),
                                                                                          ncfilename))
    for breaktime in breaktimes:
        #Create coverage plots
        meta = show_tested_gpis(workfolder, filename)
        log_file.add_line('%s: Create coverage plots with groups %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                                        str(meta)))
        #Create plots of test results with statistics
        stats = inhomo_plot_with_stats(workfolder, filename)
        log_file.add_line('%s: Create nice Test Results plot with stats for breaktime %s: %s' % (datetime.now().strftime('%Y-%m-%d_%H:%M:%S'),
                                                                                                 str(i), str(stats)))
    '''
    log_file.add_line('=====================================')
    if adjusted_ts_path:
        log_file.add_line('%s: Start TS Adjustment' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        pass
        #Do TS adjutement and save resutls to path


if __name__ == '__main__':
    # Refproduct must be one of gldas-merged,gldas-merged-from-file,merra2,ISMN-merge
    # Testproduct of form cci_*version*_*product*
    start('cci_31_combined',
          'merra2',
          r'H:\HomogeneityTesting_data\output',
          cells=[855], skip_times=None,
          anomaly=False,
          adjusted_ts_path=r'D:\users\wpreimes\datasets\CCI_adjusted',
          parallel_processes=1)
