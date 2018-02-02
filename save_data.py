# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 13:28:09 2017

@author: wpreimes
"""
import pandas as pd
import numpy as np
import pygeogrids as grids
import xarray as xr
import os
from datetime import datetime
import rsdata.root_path as root
from smecv_grid.grid import SMECV_Grid_v042

class LogFile(object):
    def __init__(self, workfolder, logfile_name=None, initial_parameters=None):
        # type: (str, dict) -> None
        self.workfolder = workfolder
        if not logfile_name:
            logfile_name = 'log.txt'
        self.log_path = os.path.join(self.workfolder, logfile_name)
        if not os.path.isfile(self.log_path):
            with open(self.log_path, 'w') as f:
                if initial_parameters:
                    f.write('Log file for HomogeneityTesting \n')
                    f.write('=====================================\n')
                    for name, val in initial_parameters.iteritems():
                        f.write('%s: %s \n' % (name, str(val)))
                    f.write('=====================================\n')
                    f.write('\n')


    def add_line(self, string, level=0):
        # type: (str) -> None
        string = '  ' * level + string
        with open(self.log_path, 'a') as f:
            f.write(string + '\n')

    def get_products(self):
        '''
        Reads the name of the reference product and the test product from log file
        :return: str, str
        '''
        test_prod, ref_prod = None, None
        with open(self.log_path, mode='r') as f:
            lines = [line.rstrip('\n') for line in f]

        for line in lines:
            if not ':' in line: continue
            if 'Test Product' in line:
                test_prod = line.split(':')[1].strip()
            elif 'Reference Product' in line:
                ref_prod = line.split(':')[1].strip()

        if not test_prod or not ref_prod:
            raise Exception('cannot find reference product or test product in log file')

        return test_prod, ref_prod


class RegularGriddedCellData(object):
    #For saving data in cellfiles (for multiprocessing)
    def __init__(self, path, grid=None, times=None,
                 reference_time=datetime(1900,1,1,0,0),
                 resolution=(0.25,0.25)):
        '''
        Cell Data Class uses xarray to write dictionaries of data to specified GPIS.

        :param path: string
            Path to which the cell files are saved
        :param grid: regular grid
            Input regular grid on which the passed GPIs are, if not passed a regular grid created
            with the passed resolution
        :param times: list
            Datetimes for which images are saved
        :param reference_time: datetime
            Reference for creating floats from datetimes
        :param resolution: tuple
            Resolution of the regular grid (lat,lon)
        '''

        self.cell_files_path = path
        if not os.path.exists(self.cell_files_path):
            os.makedirs(self.cell_files_path)

        self.times = times

        self.reference_date = reference_time

        self.global_grid = grids.genreg_grid(*resolution).to_cell_grid(5)
        if not grid:
            self.grid = self.global_grid
        else:
            self.grid = grid

        self.reset_buffer()
        #FixMe: write data to buffer before saving to file instead of writing every point directly, BufferSize Problem

    def reset_buffer(self):
        self.cell_data_cell = None
        self.cell_data_file = None
        self.cell_data_frame = None
        self.cell_data_buffersize = None

    def change_cell(self, cell):
        self.cell_data_cell = cell
        self.load_cell_data(cell)

    def load_cell_data(self, cell, type='frame'):
        #if self.cell_data_frame:
        #    print('Cell changed, write buffer to file')
        #    self.write_netcdf()

        self.cell_data_file = self.create_cellfile_name(cell=cell)
        print self.cell_data_file
        if os.path.isfile(self.cell_data_file):
            dataset = xr.open_dataset(self.cell_data_file)
            if type == 'xr':
                return dataset
            else:
                self.cell_data_frame = dataset.to_dataframe()

                self.cell_data_frame = self.cell_data_frame.reorder_levels(['time', 'lat', 'lon'])
                self.cell_data_frame = self.cell_data_frame.sort_index()
                return self.cell_data_frame

        else:
            self.cell_data_frame = self.create_empty_cell_data_frame()
            if type == 'xr':
                return self.cell_data_frame.to_xarray()
            else:
                return self.cell_data_frame

    def create_empty_cell_data_frame(self):
        gpis, lons, lats = self.grid.grid_points_for_cell(self.cell_data_cell)

        df_concat = []
        for time in self.times:
            df = pd.DataFrame(data={'time': np.tile(time, lons.size),
                                    'lon': lons, 'lat': lats})
            df_concat.append(df)
        df = pd.concat(df_concat)
        df.set_index(['time', 'lat', 'lon'], inplace=True)
        return df


    def create_empty_global_data_frame(self):
        #TODO: Delete me
        gpis, lons, lats, cells = self.global_grid.get_grid_points()

        df_concat = []
        for time in self.times:
            df = pd.DataFrame(data={'time': np.tile(time, lons.size),
                                    'lon': lons, 'lat': lats})
            df_concat.append(df)
        df = pd.concat(df_concat)
        df.set_index(['time', 'lat', 'lon'], inplace=True)
        return df


    def create_cellfile_name(self, gpi=None, cell=None):
        # Returns filename (form:cellnumber.nc) and cell for passed gpi in passed grid
        if gpi is not None:
            grid_points = self.grid.get_grid_points()
            gpi_index = np.where(grid_points[0] == gpi)[0][0]
            cell = grid_points[3][gpi_index]
            # Create filename from cell
            file_pattern = str(cell)
            while len(file_pattern) < 4:
                file_pattern = str(0) + file_pattern

            return os.path.join(self.cell_files_path, file_pattern + '.nc')

        if cell is not None:
            file_pattern = str(cell)
            while len(file_pattern) < 4:
                file_pattern = str(0) + file_pattern

            return os.path.join(self.cell_files_path, file_pattern + '.nc')

    def write_xr(self):
        return self.cell_data_frame.to_xarray()

    def write_netcdf(self):
        dataset = self.write_xr()
        dataset.to_netcdf(self.cell_data_file)
        #self.reset_buffer()


    def add_data(self, data, gpi=None, lat=None, lon=None, cell=None, time=None, write_now=True):

        if not write_now: #TODO: Implement
            raise NotImplementedError

        if not data:
            return

        if not lat or not lon:
            lon, lat = self.grid.gpi2lonlat(gpi)
        elif not gpi:
            gpi = self.grid.find_nearest_gpi(lon, lat)
        else:
            raise Exception('Select gpi and/or lon/lat')

        if not cell:
            cell = self.grid.gpi2cell(gpi)

        if time not in self.times:
            raise Exception('Time not in object definition')

        if not self.cell_data_cell or cell != self.cell_data_cell:
            self.change_cell(cell)

        data['cell'] = cell
        data['gpi'] = gpi

        index = pd.MultiIndex.from_tuples([(time, lat, lon)], names=['time', 'lat','lon'])
        gpi_data_frame = pd.DataFrame(index = index, data=data)

        if not all(np.in1d(gpi_data_frame.columns.values, self.cell_data_frame.columns.values)):
            for col in gpi_data_frame.columns.values:
                if col not in self.cell_data_frame.columns.values:
                    self.cell_data_frame[col] = np.nan
        #TODO: This raises performance warning
        self.cell_data_frame.set_value(gpi_data_frame.index,
                                       gpi_data_frame.columns.values,
                                       gpi_data_frame.loc[gpi_data_frame.index].values)


        if write_now: #TODO: Buffer writing, last cell problem,
            self.write_netcdf()

    def check_saved_cells(self):
        cells = []
        for filename in os.listdir(self.cell_files_path):
            try:
                cell = int(filename[0:4])
            except:
                continue
            if cell in self.global_grid.get_cells():
                cells.append(cell)
        return cells

    def empty_temp_files(self):
        files = os.listdir(self.cell_files_path)
        for file in files:
            os.remove(os.path.join(self.cell_files_path, file))
        os.rmdir(os.path.join(self.cell_files_path))


    def make_global_file(self, filepath=None, filename= 'global.nc', fill_nan=True, mfdataset=False,
                         keep_cell_files=False, drop_variables=None, global_meta_dict = None, var_meta_dicts=None):
        '''
        Merge all cell files in the cell files path to a global netcdf image file.
        :param filepath: string
            Path in which the global image is saved
        :param fill_nan: boolean
            Select True to read all cell files directly and write to global file without filling missing cells.
            This may lead to anomalies in the merged rendered image if the AOI is incoherent but speed up process.
        :return: None
        '''
        # TODO: Add metadata from input
        if not filepath:
            filepath = os.path.join(self.cell_files_path)

        glob_file = os.path.join(filepath, filename)

        if not fill_nan and mfdataset:
            cell_data = xr.open_mfdataset(os.path.join(self.cell_files_path,'*.nc'))
            cell_data.to_netcdf(glob_file)
        else:
            #cell_data = cell_data.to_dataframe()
            #cell_data.to_xarray().to_netcdf(filepath)

            firstfile = os.listdir(self.cell_files_path)[0]
            variables = xr.open_dataset(os.path.join(self.cell_files_path, firstfile),
                                        drop_variables=drop_variables).variables.keys()
            if fill_nan:
                self.grid = self.global_grid
            else:
                self.grid = (self.global_grid).subgrid_from_cells(self.check_saved_cells())
            #TODO: implement multiprocessing and cell buffer, avoiding reopening large global file too often

            for i, cell in enumerate(self.grid.get_cells()):
                cellfile_name = self.create_cellfile_name(cell=cell)
                if os.path.isfile(cellfile_name):
                    cell_data = xr.open_dataset(cellfile_name,
                                                drop_variables=drop_variables)
                else:
                    self.cell_data_cell = cell
                    empty_cell_data = self.load_cell_data(cell, 'frame')
                    for var in variables:
                        if var not in empty_cell_data.index.names:
                            empty_cell_data[var] = np.nan
                        cell_data = empty_cell_data.to_xarray()

                if not os.path.isfile(glob_file):
                    cell_data.to_netcdf(glob_file)
                else:
                    try:
                        with xr.open_dataset(glob_file, autoclose=True) as global_data:
                            global_data_new = xr.merge([global_data, cell_data])
                            global_data_save = global_data.copy(deep=True)
                        global_data_new.to_netcdf(glob_file, mode='w')
                    except:
                        global_data_save.to_netcdf(glob_file, mode='w')
                        print('Could not merge file for cell %s to global file' % cell)
                cell_data.close()


        if global_meta_dict:
            with xr.open_dataset(glob_file) as global_data:
                global_data_save = global_data.copy(deep=True)
                global_data_save.attrs = global_meta_dict
            global_data_save.to_netcdf(glob_file, mode='w')

        if var_meta_dicts:
            with xr.open_dataset(glob_file) as global_data:
                global_data_save = global_data.copy(deep=True)
            for varname, var_meta in var_meta_dicts.iteritems():
                if varname in global_data.variables:
                    global_data_save[varname].attrs = var_meta
            global_data_save.to_netcdf(glob_file, mode='w')


        if not keep_cell_files:
            self.empty_temp_files()

        return

def test_simple_adding():

    from smecv_grid.grid import SMECV_Grid_v042
	
    grid = SMECV_Grid_v042()
    global_grid = grids.genreg_grid(0.25, 0.25).to_cell_grid(5)

    cells = [2137, 2173, 2209]
    path = os.path.join(root.u, 'test')

    times = pd.DatetimeIndex(start='2000-01-01', end='2000-01-3', freq='D')

    # dimensions
    celldata = RegularGriddedCellData(path, grid, [t for t in times])

    for cell in cells:
        gpis, lons, lats = grid.grid_points_for_cell(cell)
        for i, gpi in enumerate(gpis):
            print gpi
            for time in times:
                data_dict = {'var1': np.random.rand(1)[0],
                             'var2': np.random.rand(1)[0],
                             'var3': np.random.rand(1)[0]}

                celldata.add_data(data_dict, gpi=gpi, time=time)

    celldata.make_global_file(meta_dicts={'var1':{'m1':1, 'm2':2}})

if __name__ == '__main__':
    from cci_timeframes import CCITimes
    from smecv_grid.grid import SMECV_Grid_v042
    from interface import get_test_meta
    from adjustment import get_adjustment_meta


    workfolder =r'D:\users\wpreimes\datasets\HomogeneityTesting_data\v38'
    test_prod = 'CCI_41_COMBINED'
    skip_times = None
    save_adjusted_data = True
    times_obj = CCITimes(test_prod, ignore_position=True,
                         skip_times=skip_times)
    results = ['test_results_before_adjustment',
               'lin_model_params_before_adjustment',
               'test_results_after_adjustment',
               'lin_model_params_after_adjustment']
    save_objects = dict.fromkeys(results)
    for name in save_objects.keys():
        save_objects[name] = RegularGriddedCellData(grid=SMECV_Grid_v042(),
                                                    path=os.path.join(workfolder, 'temp_homogtest_cellfiles', name),
                                                    times=times_obj.get_times(None, as_datetime=True)['breaktimes'],
                                                    resolution=(0.25, 0.25))
    for name, save_obj in save_objects.iteritems():
        #if name != 'lin_model_params_after_adjustment':continue
        filename = 'GLOBAL_' + name + '.nc'
        var_meta = {'test_results': get_test_meta()[0],
                     'test_status': get_test_meta()[1],
                     'adjustment_status': get_adjustment_meta()}
        global_meta = {'test_prod': 'CCI_41_COMBINED',
                       'ref_prod': 'merra2'}

        save_obj.make_global_file(workfolder, filename, False, False, True, None,global_meta, var_meta)




