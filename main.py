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

from interface import HomogTest
from cci_timeframes import get_timeframes
from save_data import SaveResults, save_Log
from grid_functions import grid_points_for_cells
from otherplots import show_tested_gpis, inhomo_plot_with_stats
from longest_homogeneous_period import calc_longest_homogeneous_period
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

def test_for_cells(q, workfolder, grid, cells, log_file, test_prod, ref_prod, anomaly):
    #Function for multicore processing

    test_obj = HomogTest(test_prod,
                         ref_prod,
                         0.01,
                         anomaly)

    filename = 'HomogeneityTest_%s_%s' % (test_obj.test_prod, test_obj.ref_prod)


    save_obj = SaveResults(workfolder, grid, filename, buffer_size=300)

    for cell in cells:
        log_file.add_line('%s: Start Testing Cell %i' % (str(datetime.now()), cell))
        grid_points = grid.grid_points_for_cell(cell)

        for iteration, gpi in enumerate(grid_points):
            if iteration % 1000 == 0:
                print 'Processing QDEG Point %i (iteration %i of %i)' % (gpi, iteration, len(grid_points))
            '''
            if test_obj.ref_prod == 'ISMN-merge':
                valid_insitu_gpis = test_obj.ismndata.gpis_with_netsta

                if gpi not in valid_insitu_gpis.keys():
                    continue
            '''
            # test_obj.save_as_mat(gpi=gpi)
            try:
                df_time = test_obj.read_gpi(gpi, start=test_obj.range[0], end=test_obj.range[1])
                for timeframe, breaktime in zip(test_obj.timeframes, test_obj.breaktimes):
                    data = df_time[timeframe[0]:timeframe[1]] #TODO: dropna here??

                    df_results, testresult = test_obj.run_tests(data=data,
                                                                breaktime=breaktime,
                                                                min_data_size=3,    #TODO: differnt value?
                                                                tests=['wk', 'fk'])
                    # wenn testresults break finden adjustment über timeframe
                    # Replace df_time[timeframe] with adjusted data

                #Save adjusted TS afterwards!
                testresult.update({'status': '0: Testing successful'})
            except Exception as e:
                df_results = None
                testresult = {'status': str(e)}

            save_obj.fill_buffer(gpi, testresult)
            if iteration == len(grid_points) - 1:
                # The last group of grid points is saved to file, also if the buffer size is not yet reached
                save_obj.save_to_netcdf()

        # Add Info to log file
        log_file.add_line('%s: Finished testing for cell %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), cell))
    log_file.add_line('Saved results to: HomogeneityTest_%s.nc' % breaktime)
    #Return filename
    q.put(filename)


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
    
    if cells == 'global':
        cells = grid.get_cells()
    
    workfolder = create_workfolder(path)
       
    log_file = save_Log(workfolder, test_prod, ref_prod, anomaly, cells)

    q = Queue()
    for cell_pack in list(split(cells, parallel_processes)): # Split cells in equally sized packs for multiprocessing
        p = Process(target=test_for_cells, args=(q, workfolder, grid, cell_pack,
                                                 log_file, test_prod, ref_prod, anomaly))
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
        log_file.add_line('%s: Create coverage plots with groups %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                                        str(meta)))
        #Create plots of test results with statistics
        stats = inhomo_plot_with_stats(workfolder, filename)
        log_file.add_line('%s: Create nice Test Results plot with stats for breaktime %s: %s' % (datetime.now().strftime('%Y-%m-%d_%H:%M:%S'),
                                                                                                 str(i), str(stats)))
    log_file.add_line('=====================================')
    if adjusted_ts_path:
        log_file.add_line('%s: Start TS Adjustment' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        pass
        #Do TS adjutement and save resutls to path


if __name__ == '__main__':
    # Refproduct must be one of gldas-merged,gldas-merged-from-file,merra2,ISMN-merge
    # Testproduct of form cci_*version*_*product*
    start('adjusted_cci',
          'merra2',
          r'H:\HomogeneityTesting_data\output',
          cells='Australia', skip_times=None,
          anomaly=None,
          adjusted_ts_path=None,
          parallel_processes=8)
