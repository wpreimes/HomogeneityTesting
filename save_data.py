# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 13:28:09 2017

@author: wpreimes
"""
import pandas as pd
import numpy as np
from pygeogrids.grids import CellGrid
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
            f.write('Processed Cells: %s \n' % cells)
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



class Results_2D(object):
    '''
    Baiscally a big DataFrame over all passed GPIs with function to flush and save data
    '''
    def __init__(self, path, grid, breaktimes, buffer_size=400):
        # type: (str, CellGrid, list, int) -> None

        if not all(isinstance(breaktime, str) for breaktime in breaktimes):
            raise Exception("Breaktimes must be passed as list of strings")

        self.total_number = [0, ]
        self.path = path
        self.grid = grid
        self.buffer_size = buffer_size
        self.breaktimes = breaktimes
        self.gpis = np.array([])
        self.data = {breaktime: {'gpi':[]} for breaktime in self.breaktimes}
        self.fn_global = 'global_points_homogtest.nc'
        self.fn_images_base = 'HomogeneityTest_%s'
    def reset_data(self):
        self.data = {breaktime: {'gpi':[]} for breaktime in self.breaktimes}

    def add_data(self, gpi, breaktime, data_dict):
        '''
        Fill buffer by adding data for gpi
        :param gpi:
        :param breaktime: str
        :param data_dict:
        :return:
        '''

        if breaktime not in self.data.keys():
            raise Exception("Breaktime not in Buffer Object")
        if gpi in self.data[breaktime]['gpi']:
            raise Exception("GPI already stored")
        else:
            data_dict = self.extract_infos(data_dict)
            self.data[breaktime]['gpi'].append(gpi)

            for name, value in data_dict.iteritems():
                if name not in self.data[breaktime].keys():
                    self.data[breaktime][name] = []
                self.data[breaktime][name].append(value)

            #Save the buffer to file as soon as ANY breaktime reaches the max size
            if len(self.data[breaktime]['gpi']) == self.buffer_size:
                self.save_to_netcdf()

    def extract_infos(self, data_dict):
        # type: (dict) -> dict
        status = int(data_dict['status'][0])
        if status not in [0, 6, 7]:
            return {'h_wk': np.nan, 'h_fk': np.nan, 'test_results': np.nan, 'status': status}
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

            return {'h_wk': wk, 'h_fk': fk, 'test_results': all, 'status': status}


    def save_to_netcdf(self):

        file_meta_dict = {
            'test_results': 'Break detection classes by Hopothesis tests.'
                            '1 = WK only, 2 = FK only, 3 = WK and FK, 4 = None',
            'status': 'Processing status and errors:'
                      '0 = Processing OK,'
                      '1 = No coinciding data for timeframe,'
                      '2 = Test TS and Ref TS do not match,'
                      '3 = Spearman correlation too low,'
                      '4 = Minimum Dataseries Length not reached,'
                      '5 = neg/nan correl. aft. bias corr.,'
                      '6 = Error during WK testing,'
                      '7 = Error during FK testing,'
                      '8 = WK test and FK test failed,'
                      '9 = Error reading gpi'}

        dataframes = []
        for breaktime, data_dict in self.data.iteritems():
            df = pd.DataFrame.from_dict(data_dict)
            df = df.set_index('gpi')
            df = df.rename(columns={column_name : column_name + '_' + breaktime for column_name in df.columns.values})
            dataframes.append(df)

        df = pd.concat(dataframes, axis=1)
        gpi_vals = df.to_dict('index')

        with GriddedPointData(self.path, mode='a', grid=self.grid,
                              fn_format='{:04d}.nc') as nc:

            for loc_id, data_dict in gpi_vals.iteritems():
                nc.write(loc_id, data_dict)
                if loc_id not in self.gpis:
                    self.gpis = np.append(self.gpis, loc_id)


        self.gpis = np.append(self.gpis, df.index.values)
        self.reset_data()

    def save_global_file(self):

        grid = self.grid.subgrid_from_gpis(np.unique(self.gpis))
        with GriddedPointData(self.path, grid=grid,
                              fn_format='{:04d}.nc') as nc:

            nc.to_point_data(os.path.join(self.path, self.fn_global))

        return self.fn_global

    def create_image_files(self):
        global_file = xr.open_dataset(os.path.join(self.path, self.fn_global))
        df = global_file.to_dataframe()

        filenames =[]
        for breaktime_str in self.breaktimes:
            df_breaktime = df[['lat','lon']]
            for col_name in df.columns.values:
                if breaktime_str in col_name:
                    [var, breaktime] = re.split(r'[_](?=[0-9])', col_name)
                    df_breaktime[var] = df[col_name]

            df_breaktime = df_breaktime.sort_values(['lat','lon']) \
                                       .set_index(['lat','lon'])
            df_breaktime.index
            df_breaktime.index=df_breaktime.index.drop_duplicates()
            print df_breaktime
            global_image = df_breaktime.to_xarray()
            filename = self.fn_images_base % breaktime_str
            global_image.to_netcdf(os.path.join(self.path, filename))
            filenames.append(filename)
        return filenames



'''     
test_obj.DF_Points['test_results']=np.nan
test_obj.DF_Points['message']='Not processed'


test_obj.save_test_results(gpi,'h_fk',testresult['FlignerKilleen']['h'])
test_obj.save_test_results.set_value(gpi,'h_wk',testresult['Wilkoxon']['h'])

test_obj.save_test_results.set_value(gpi,'message','Processing OK')  
        
#In case something went wrong
test_obj.DF_Points.set_value(gpi,'h_fk',np.nan)
test_obj.DF_Points.set_value(gpi,'h_wk',np.nan)
test_obj.DF_Points.set_value(gpi,'message',str(e))
'''
