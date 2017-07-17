# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 17:05:55 2017

@author: wpreimes
"""

import sys
if r"H:\workspace" not in sys.path:
    sys.path.append(r"H:\workspace")

from typing import Union
import os
from datetime import datetime
from pygeogrids.netcdf import load_grid
from multiprocessing import Process, Queue

from HomogeneityTesting.interface import HomogTest
from HomogeneityTesting.otherfunctions import cci_timeframes
from HomogeneityTesting.save_data import SaveResults, save_Log
from HomogeneityTesting.grid_functions import grid_points_for_cells
from HomogeneityTesting.otherplots import show_tested_gpis, inhomo_plot_with_stats
from HomogeneityTesting.longest_homogeneous_period import calc_longest_homogeneous_period
'''
Breaktimes sind Punkte an denen Inhomogenitäten vermutet werden
Laut Processing overview und Research letter:
    Combined Data: (ignore September 1987) , August 1991, January 1998, July 2002, Januar 2007, 
                    October 2011, July 2012, (additional May 2015)

used alpha in research letter:0.01
'''

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

def single_test_for_breaktime(q, workfolder, grid, grid_points, log_file, test_prod, ref_prod, timeframe, breaktime, anomaly):
    #Function for multicore processing

    test_obj = HomogTest(test_prod,
                         ref_prod,
                         timeframe,
                         breaktime,
                         0.01,
                         anomaly)

    filename = 'HomogeneityTest_%s_%s' % (test_obj.ref_prod, test_obj.breaktime.strftime("%Y-%m-%d"))


    save_obj = SaveResults(workfolder, grid, filename, buffer_size=300)

    log_file.add_line('%s: Start Testing Timeframe and Breaktime: %s and %s'
                      % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), timeframe, breaktime))

    for iteration, gpi in enumerate(grid_points):
        if iteration % 1000 == 0:
            print 'Processing QDEG Point %i (iteration %i of %i)' % (gpi, iteration, len(grid_points))

        if test_obj.ref_prod == 'ISMN-merge':
            valid_insitu_gpis = test_obj.ismndata.gpis_with_netsta

            if gpi not in valid_insitu_gpis.keys():
                continue

        try:
            # test_obj.save_as_mat(gpi=gpi)
            df_time, testresult = test_obj.run_tests(gpi=gpi,
                                                     tests=['wk', 'fk'])
            testresult.update({'status': '0: Testing successful'})
        except Exception as e:
            df_time = None
            testresult = {'status': str(e)}

        save_obj.fill_buffer(gpi, testresult)
        if iteration == len(grid_points) - 1:
            # The last group of grid points is saved to file, also if the buffer size is not yet reached
            save_obj.save_to_netcdf()

    # Add Info to log file
    log_file.add_line('%s: Finished testing for timeframe %s' % (datetime.now().strftime('%Y-%m-%d_%H:%M:%S'),
                                                                 timeframe))
    log_file.add_line('Saved results to: HomogeneityTest_%s.csv' % breaktime)
    #Return filename
    q.put(filename)

def start(test_prod, ref_prod, path, cells='global',skip_times=None, anomaly=False, adjusted_ts_path=None):
    # type: (str, str, str, Union[list,str], Union[list,None], bool) -> None

    testtimes = cci_timeframes(test_prod, skip_times=skip_times)

    timeframes = testtimes['timeframes']
    breaktimes = testtimes['breaktimes']

    
    grid = load_grid(r"D:\users\wpreimes\datasets\grids\qdeg_land_grid.nc")
    
    if cells == 'global' or not cells:
        grid_points = grid.get_grid_points()[0]
    else:
        grid_points = grid_points_for_cells(grid, cells)
    
    workfolder = create_workfolder(path)
       
    log_file = save_Log(workfolder, test_prod, ref_prod, anomaly, cells)

    q = Queue()
    for breaktime, timeframe in zip(breaktimes, timeframes):
        p = Process(target=single_test_for_breaktime, args=(q, workfolder, grid, grid_points,
                                                            log_file, test_prod, ref_prod, timeframe,
                                                            breaktime, anomaly))
        p.start()

    filenames = []
    for i in range(len(breaktimes)):
        filenames.append(q.get(True))

    log_file.add_line('=====================================')
    #Plotting and netcdf files
    for i, filename in enumerate(filenames):
        # Create nc file with longest period data
        ncfilename = calc_longest_homogeneous_period(workfolder, create_netcdf=True)
        log_file.add_line('%s: Saved longest homogeneous Period, startdate, enddate to %s' % (datetime.now().strftime('%Y-%m-%d_%H:%M:%S'),
                                                                                              ncfilename))
        #Create coverage plots
        meta = show_tested_gpis(workfolder, filename)
        log_file.add_line('%s: Create coverage plots with groups %s' % (datetime.now().strftime('%Y-%m-%d_%H:%M:%S'),
                                                                        str(meta)))
        #Create plots of test results with statistics
        stats = inhomo_plot_with_stats(workfolder, filename)
        log_file.add_line('%s: Create nice Test Results plot with stats for breaktime %s: %s' % (datetime.now().strftime('%Y-%m-%d_%H:%M:%S'),
                                                                                                 str(i), str(stats)))
    log_file.add_line('=====================================')
    if adjusted_ts_path:
        log_file.add_line('%s: Start TS Adjustment' % (datetime.now().strftime('%Y-%m-%d%H:%M:%S')))
        pass
        #Do TS adjutement and save resutls to path


if __name__ == '__main__':
    # Refproduct must be one of gldas-merged,gldas-merged-from-file,merra2,ISMN-merge
    # Testproduct of form cci_*version*_*product*
    start('cci_31_passive',
          'merra2',
          r'H:\HomogeneityTesting_data\output',
          cells='global', skip_times=None,
          anomaly=False,
          adjusted_ts_path=None)
