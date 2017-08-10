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
    '''
    Baiscally a big DataFrame over all passed GPIs with function to flush and save data
    '''
    def __init__(self, grid, path, breaktimes):
        # type: (grids.CellGrid,str, list, int) -> None

        if not all(isinstance(breaktime, str) for breaktime in breaktimes):
            raise Exception("Breaktimes must be passed as list of strings")

        self.path = path
        cell_files_path = os.path.join(self.path, 'gridded_files')
        if not os.path.isdir(cell_files_path):
            os.mkdir(cell_files_path)
        self.cell_files_path = cell_files_path
        self.breaktimes = breaktimes
        self.pre_process_grid = grid # type: grids.CellGrid
        self.global_gpis=[]
        self.data = {breaktime: {'gpi': []} for breaktime in self.breaktimes}

    def save_subgrid(self, gpis):
        filename = 'breaktest_grid.nc'
        subgrid = self.pre_process_grid.subgrid_from_gpis(gpis)

        nc.save_grid(os.path.join(self.cell_files_path, filename),
                     grid=subgrid, subset_name='breaktest_flag',
                     subset_meaning='Tested Points')
        return filename


    def reset_data(self):
        self.data = {breaktime: {'gpi': []} for breaktime in self.breaktimes}

    def add_data(self, gpi, breaktime, gpi_result):
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
            data_dict = {}
            for name, resdict in gpi_result.iteritems():
                if name == 'testresult':
                    data_dict.update(self.extract_infos_test(resdict))
                if name == 'adjresult':
                    data_dict.update(self.extract_infos_adj(resdict))

            self.data[breaktime]['gpi'].append(gpi)

            for name, value in data_dict.iteritems():
                if name not in self.data[breaktime].keys():
                    self.data[breaktime][name] = []
                self.data[breaktime][name].append(value)


    def extract_infos_test(self, data_dict):
        # type: (dict) -> dict
        status = int(data_dict['status_test'][0])
        if status not in [0, 6, 7]:
            return {'h_wk': np.nan, 'h_fk': np.nan, 'test_results': np.nan, 'status_test': status}
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

            return {'h_wk': wk, 'h_fk': fk, 'test_results': all, 'status_test': status}

    def extract_infos_adj(self, data_dict):
        # type: (dict) -> dict
        status = int(data_dict['status_adj'][0])
        if status not in [0]:
            return {'slope': np.nan, 'intercept': np.nan, 'part1_B1': np.nan, 'part1_B2': np.nan,
                   'part2_B1': np.nan, 'part2_B2': np.nan,
                   'B1_aft_adjust': np.nan, 'B2_aft_adjust': np.nan, 'status_adj': status}
        else:
            data_dict['status_adj'] = status
        return data_dict


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



class GlobalResults(object):

    def __init__(self, path, cellfiles_folder, breaktimes):
        self.path = path
        self.cellfile_dir = os.path.join(path, cellfiles_folder)
        self.breaktimes = breaktimes
        self.fn_global = 'homogtest_global.nc'
        self.fn_images_base = 'HomogeneityTestResult_%s_image.nc'
        self.gpis = nc.load_grid(os.path.join(self.cellfile_dir,'breaktest_grid.nc')).get_grid_points()[0]

        landgrid = nc.load_grid(r"D:\users\wpreimes\datasets\grids\qdeg_land_grid.nc")
        self.post_process_grid = landgrid.subgrid_from_gpis(self.gpis)

    def save_global_file(self, keep_cell_files=False):

        with GriddedPointData(self.cellfile_dir, grid=self.post_process_grid,
                              fn_format='{:04d}.nc') as nc:

            nc.to_point_data(os.path.join(self.path, self.fn_global))

        if not keep_cell_files:
            for cell in np.unique(self.post_process_grid.get_grid_points()[3]):
                filename = str(cell)
                while len(filename) < 4:
                    filename = str(0) + filename
                os.remove(os.path.join(self.cellfile_dir, filename+'.nc'))
        return self.fn_global

    def create_image_files(self):

        file_meta_dict = {
            'test_results': 'Break detection classes by Hopothesis tests.'
                            '1 = WK only, 2 = FK only, 3 = WK and FK, 4 = None',
            'status_test': '0 = Processing OK,'
                      '1 = No coinciding data for timeframe,'
                      '2 = Test TS and Ref TS do not match,'
                      '3 = Spearman correlation too low,'
                      '4 = Min. Dataseries len. not reached,'
                      '5 = neg/nan correl. aft. bias corr.,'
                      '6 = Error during WK testing,'
                      '7 = Error during FK testing,'
                      '8 = WK test and FK test failed,'
                      '9 = Error reading gpi'}

           # 'status_adj': '0: Adjustment performed,'
           #         '1 = No adjustment necessary,'
           #         '2 = negative correlation,'
           #         '3 = positive correlation insignificant'}

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
    grid = nc.load_grid(r"D:\users\wpreimes\datasets\grids\qdeg_land_grid.nc")
    save_obj = Results2D(grid, r'C:\Temp', ['2000-01-01'])
    for gpi in gpis:
        save_obj.add_data(gpi, '2000-01-01', {'testresult':{'status_test':'1'},
                                              'adjresult': {'status_adj': '1'}})
    global_gpis = save_obj.save_to_gridded_netcdf()

    glob_obj = GlobalResults(r'C:\Temp', r'C:\Temp\gridded_files', ['2000-01-01'])
    global_file_name = glob_obj.save_global_file(False)
    image_files_names = glob_obj.create_image_files()
