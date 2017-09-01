# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 13:28:09 2017

@author: wpreimes
"""
import pandas as pd
import numpy as np
import pygeogrids as grids
import pygeogrids.netcdf as nc
import xarray as xr
from pynetcf.point_data import GriddedPointData
import os
import re
import shutil
from smecv_grid.grid import SMECV_Grid_v042


class save_Log(object):
    def __init__(self, workfolder, test_prod, ref_prod, anomaly, cells):
        # type: (str, str, str, str, list) -> None
        self.workfolder = workfolder

        with open(os.path.join(workfolder, 'log.txt'), 'w') as f:
            f.write('Log file for HomogeneityTesting \n')
            f.write('=====================================\n')
            f.write('Test Product: %s \n' % test_prod)
            f.write('Reference Product: %s \n' % ref_prod)
            f.write('Anomaly data: %s \n' % str(anomaly))
            f.write('Input Cells: %s \n' % cells)
            f.write('=====================================\n')
            f.write('\n')

    def add_line(self, string):
        # type: (str) -> None
        # Initialize and add line to log file for current process
        with open(os.path.join(self.workfolder, 'log.txt'), 'a') as f:
            f.write(string + '\n')

class load_Log(object):
    def __init__(self, workfolder):
        self.workfolder = workfolder

    def get_products(self):
        with open(os.path.join(self.workfolder,'log.txt'), mode='r') as f:
            lines = [line.rstrip('\n') for line in f]

        test_prod = lines[2].split(' ')[2]
        ref_prod = lines[3].split(' ')[2]

        return {'test_prod': test_prod, 'ref_prod': ref_prod}

class Results2D(object):
    def __init__(self, grid, path, breaktimes, cell_files_path=None):
        if not all(isinstance(breaktime, str) for breaktime in breaktimes):
            raise Exception("Breaktimes must be passed as list of strings")

        self.path = path
        if not cell_files_path:
            cell_files_path = os.path.join(self.path, 'gridded_files')
            if not os.path.isdir(cell_files_path):
                os.mkdir(cell_files_path)
        self.cell_files_path = cell_files_path
        self.breaktimes = breaktimes
        self.pre_process_grid = grid  # type: grids.CellGrid
        self.global_gpis = []
        self.data = {breaktime: {'gpi': []} for breaktime in self.breaktimes}

    def empty_temp_files(self):
        files = os.listdir(self.cell_files_path)
        for file in files:
            os.remove(os.path.join(self.cell_files_path, file))

        return self.cell_files_path

    def save_subgrid(self, gpis):
        filename = 'breaktest_grid.nc'
        subgrid = self.pre_process_grid.subgrid_from_gpis(gpis)

        nc.save_grid(os.path.join(self.path, filename),
                     grid=subgrid, subset_name='breaktest_flag',
                     subset_meaning='Tested Points')

        return subgrid

    def reset_data(self):
        self.data = {breaktime: {'gpi': []} for breaktime in self.breaktimes}


    def add_data(self, gpi, breaktime, data_dict):
        '''
        Fill buffer by adding data for gpi
        :param gpi:
        :param breaktime: str
        :param data_dict:
        :return:
        '''
        if gpi in self.data[breaktime]['gpi']:
            raise Exception("GPI already stored")
        else:
            self.data[breaktime]['gpi'].append(gpi)
            for name, value in data_dict.iteritems():
                if name not in self.data[breaktime].keys():
                    self.data[breaktime][name] = []
                self.data[breaktime][name].append(value)

    def save_to_gridded_netcdf(self):

        dataframes = []

        for breaktime, data_dict in self.data.iteritems():
            df = pd.DataFrame.from_dict(data_dict)
            df = df.set_index('gpi')
            df = df.rename(columns={column_name : column_name + '_' + breaktime for column_name in df.columns.values})
            dataframes.append(df)

        df = pd.concat(dataframes, axis=1) #type: pd.DataFrame

        if not df.index.is_unique:
            raise Exception('df gpis not unique')

        gpi_vals = df.to_dict('index')
        grid = self.pre_process_grid.subgrid_from_gpis(gpi_vals.keys())

        with GriddedPointData(self.cell_files_path, mode='a', grid=grid,
                              fn_format='{:04d}.nc') as nc:

            for loc_id, data_dict in gpi_vals.iteritems():
                nc.write(loc_id, data_dict)
                if loc_id not in self.global_gpis:
                    self.global_gpis.append(loc_id)
                else:
                    raise Exception('Trying to set existing GPI %i' %loc_id)
        return_gpis = self.global_gpis
        self.reset_data()
        self.global_gpis=[]
        return return_gpis


def extract_test_infos(data_dict):
    # type: (dict) -> dict
    status = int(data_dict['test_status'][0])
    if status not in [0, 6, 7]:
        return {'h_wk': np.nan, 'h_fk': np.nan, 'test_results': np.nan, 'test_status': status}
    else:
        wk = data_dict['wilkoxon']['h']
        fk = data_dict['fligner_killeen']['h']

        if type(wk) is str:
            # wk test failed and contains error message
            wk = np.nan
            status = 6 # This is just an identification number, no math meaning
        if type(fk) is str:
            # fk failed and contains error message
            fk = np.nan
            status = 7 # This is just an identification number, no math meaning

        all = 4.0
        if wk == 1:
            if fk == 1:
                all = 3.0
            else:
                all = 1.0
        elif fk == 1:
            all = 2.0

        return {'h_wk': wk, 'h_fk': fk, 'test_results': all, 'test_status': status}


def extract_adjust_infos(data_dict):
    # type: (dict) -> dict
    status = int(data_dict['adj_status'][0])
    if status not in [0]:
        return {'slope': np.nan, 'intercept': np.nan, 'part1_B1': np.nan, 'part1_B2': np.nan,
                'part2_B1': np.nan, 'part2_B2': np.nan, 'adj_status': status}
    else:
        data_dict['adj_status'] = status

        return data_dict  #DataDict from Adjustment incl. status



class GlobalResults(Results2D):

    def __init__(self, baseObject, post_process_grid):
        self.__class__ = type(baseObject.__class__.__name__,
                              (self.__class__, baseObject.__class__),
                              {})

        self.path = baseObject.path
        self.cell_files_path = baseObject.cell_files_path
        self.breaktimes = baseObject.breaktimes
        self.fn_global = 'homogtest_global.nc'
        self.fn_images_base = 'HomogeneityTestResult_%s_image.nc'

        self.post_process_grid = post_process_grid


    def save_global_file(self, keep_cell_files=False):

        with GriddedPointData(self.cell_files_path, grid=self.post_process_grid,
                              fn_format='{:04d}.nc') as nc:

            nc.to_point_data(os.path.join(self.path, self.fn_global))

        if not keep_cell_files:
            self.empty_temp_files()
        return self.fn_global

    def create_image_files(self):

        file_meta_dict = {
            'test_results': 'Break detection classes by Hopothesis tests.'
                            '1 = WK only, 2 = FK only, 3 = WK and FK, 4 = None',
            'test_status': '0 = Testing performed,'
                      '1 = No coinciding data for timeframe,'
                      '2 = Test TS and Ref TS do not match,'
                      '3 = Spearman correlation too low,'
                      '4 = Min. Dataseries len. not reached,'
                      '5 = neg/nan correl. aft. bias corr.,'
                      '6 = Error during WK testing,'
                      '7 = Error during FK testing,'
                      '8 = WK test and FK test failed,'
                      '9 = Error reading gpi',
            'adj_status': '0: Adjustment performed,'
                     '1 = No adjustment necessary,'
                     '2 = negative correlation,'
                     '3 = positive correlation insignificant,'
                     '4 = Model param 1 not sufficiently adapted,'
                     '5 = Model param 2 not sufficiently adapted,'
                     '6 = B tolerance after adj. not reached,'
                     '9 = Not adjusted' }

        global_file = xr.open_dataset(os.path.join(self.path, self.fn_global))
        df = global_file.to_dataframe()

        filenames =[]
        for breaktime_str in self.breaktimes:
            df_breaktime = df[['lat','lon']]
            for col_name in df.columns.values:
                if breaktime_str in col_name:
                    [var, breaktime] = re.split(r'[_](?=[0-9])', col_name)
                    df_breaktime[var] = df[col_name]
            df_breaktime['location_id'] = df['location_id'].astype(int)
            df_breaktime = df_breaktime.sort_values(['lat','lon']) \
                                       .set_index(['lat','lon'])

            global_image = df_breaktime.to_xarray()

            for name, val in file_meta_dict.iteritems():
                global_image.variables[name].attrs['Values']=val

            filename = self.fn_images_base % breaktime_str
            global_image.to_netcdf(os.path.join(self.path, filename))

            filenames.append(filename)
        return filenames


if __name__ == '__main__':
    gpis = range(392888,392901)+range(391448,391461)
    grid = SMECV_Grid_v042()
    save_obj = Results2D(grid, r'C:\Temp', ['2000-01-01'])
    for gpi in gpis:
        save_obj.add_data(gpi, '2000-01-01', {'testresult':{'status_test':'1'},
                                              'adjresult': {'status_adj': '1'}})
    global_gpis = save_obj.save_to_gridded_netcdf()

    glob_obj = GlobalResults(r'C:\Temp', r'C:\Temp\gridded_files', ['2000-01-01'])
    global_file_name = glob_obj.save_global_file(False)
    image_files_names = glob_obj.create_image_files()
