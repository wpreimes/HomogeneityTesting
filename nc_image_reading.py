# -*- coding: utf-8 -*-
"""
Created on Jän 16 13:46 2018

@author: wpreimes
"""
import xarray as xr
import pandas as pd


class RegularGriddedCellData_Reader(object):
    #For saving data in cellfiles (for multiprocessing)
    def __init__(self, nc_image_file):
        '''
        Wrapper Class for reading 2D images at a certain time stamp from an nc file
        :param nc_image_file: path to the netcdf image file
        '''
        self.ncfile = xr.open_dataset(nc_image_file)

    def read(self, time, var_names=None):
        DF_Points_from_file = self.ncfile.sel(time=time).to_dataframe() \
            .reset_index(inplace=False).set_index(['gpi'])

        DF_Points_from_file = DF_Points_from_file.loc[DF_Points_from_file.index.dropna(),:]
        if var_names:
            DF_Points_from_file = DF_Points_from_file.loc[:, var_names]
        return DF_Points_from_file

    def get_variables(self):
        return self.ncfile.variables.keys()

