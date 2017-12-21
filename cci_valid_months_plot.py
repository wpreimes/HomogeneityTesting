# -*- coding: utf-8 -*-
"""
Created on Dez 04 11:15 2017

@author: wpreimes
"""
from cci_timeframes import CCITimes
from smecv_grid.grid import SMECV_Grid_v042
from interface import BreakTestData
import pandas as pd
import numpy as np
from datetime import datetime
import os
from main import split_cells
from multiprocessing import Process, Queue
from otherfunctions import split
from save_data import RegularGriddedCellData


def process_for_cells(q, grid, time_obj, test_obj, save_obj, cells_subset, resample_method, min_monthly_values,
                      process_no):

    for cell in grid.subgrid_from_cells(cells_subset).get_cells():
        grid = grid.subgrid_from_cells(cells_subset)
        gpis = grid.grid_points_for_cell(cells_subset)[0]
        lons = grid.grid_points_for_cell(cells_subset)[1]
        lats = grid.grid_points_for_cell(cells_subset)[2]

        DF_Points = pd.DataFrame(index=gpis, data={'lon': lons, 'lat': lats})
        for breaktime in time_obj.get_times()['breaktimes']:
            DF_Points['before %s' % breaktime] = np.nan
            DF_Points['after %s' % breaktime] = np.nan

        for iteration, gpi in enumerate(grid.grid_points_for_cell(cell)[0]):
            times = time_obj.get_times(as_datetime=True)
            starts = np.flipud(np.array([timeframe[0] for timeframe in times['timeframes']]))
            ends = np.flipud(np.array([timeframe[1] for timeframe in times['timeframes']]))
            breaktimes = np.flipud(times['breaktimes'])

            print 'processing gpi %i (iteration %i of %i)' % (gpi, iteration, len(grid.grid_points_for_cell(cell)[0]))
            try:
                df_time = test_obj.read_gpi(gpi, starts[0], ends[-1])

                for start, breaktime, end in zip(starts, breaktimes, ends):

                    df_subtime = df_time[['testdata']][start:end].dropna()
                    df_resample, _ = test_obj.temp_resample(df_subtime, resample_method, min_monthly_values)
                    df_group, len_bef, len_aft = \
                        test_obj.group_by_breaktime(df_resample, breaktime, 3, ignore_exception=True)

                    data_dict = {'before_%s' % str(breaktime.date()):len_bef,
                                 'after_%s' % str(breaktime.date()):len_aft}

                    save_obj.add_data(data_dict, gpi=gpi, time=datetime(2000, 1, 1)) #Placeholder time

            except:
                continue

    q.put(True)

def valid_months_plot(workdir, cells_identifier, testproduct, refproduct, resample_method, min_monthly_values, parallel_processes):
    # NOT IN PROCESS

    # Count the number of valid months i.e, month that contain >=
    # the minimum amount of values
    # Count values BEFORE breaktime
    # Count value AFTER breaktime

    time_obj = CCITimes(testproduct, ignore_position=True)
    grid = SMECV_Grid_v042()

    test_obj = BreakTestData(testproduct, refproduct, False)

    save_obj = RegularGriddedCellData(grid=grid,
                                      path=os.path.join(workdir, 'temp_cellfiles'),
                                      times=[datetime(2000, 1, 1)], #Placeholder time
                                      resolution=(0.25, 0.25))


    cells = list(
        split(split_cells(cells_identifier, grid), parallel_processes))  # split cells for processes

    processes = []
    q = Queue()
    finished_processes = []

    for process_no in range(parallel_processes):
        cells_for_process = cells[process_no]
        p = Process(target=process_for_cells, args=(q,
                                                    grid,
                                                    time_obj,
                                                    test_obj,
                                                    save_obj,
                                                    cells_for_process,
                                                    resample_method,
                                                    min_monthly_values,
                                                    process_no))
        processes.append(p)
        p.start()

    for process in processes:
        finished_processes.append(q.get(True))

    for process in processes:
        process.join()

    save_obj.make_global_file(workdir,
                              'temp_coverage.nc',
                              True,
                              False,
                              True,
                              None)




if __name__ == '__main__':
    valid_months_plot(r'D:\users\wpreimes\datasets\HomogeneityTesting_data\cci41_validmonths',
                      'global',
                      'CCI_41_COMBINED',
                      'merra2',
                      'M',
                      0.33,
                      8)