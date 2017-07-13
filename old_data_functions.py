# -*- coding: utf-8 -*-
"""
Created on Jul 13 12:15 2017

@author: wpreimes
"""

import glob
import os
import pandas as pd
import numpy as np
from points_to_netcdf import points_to_netcdf, globalCellgrid
from netCDF4 import Dataset


def csv_to_netcdf(workdir):
    # type: (str) -> None
    '''
    Turns csv file from old version of homogneity testing into netcdf file
    :param workdir:
    :return: None
    '''
    fileslist = glob.glob(os.path.join(workdir, "DF_Points_merra2*"))


    for filename in fileslist:
        DF_Points = pd.read_csv(filename, index_col=0)

        filename = filename.replace(workdir, '')
        filename = filename.replace('.', '_')
        splitname = filename.split('_')
        breaktime = splitname[4]

        if 'h_all' in DF_Points.columns.values:
            DF_Points = DF_Points.rename(columns={'h_all': 'test_results'})
        if 'h_WK' in DF_Points.columns.values:
            DF_Points = DF_Points.rename(columns={'h_WK': 'h_wk'})
        if 'h_FK' in DF_Points.columns.values:
            DF_Points = DF_Points.rename(columns={'h_FK': 'h_fk'})

        points_to_netcdf(dataframe=DF_Points[['test_results', 'h_wk', 'h_fk']], path=workdir, filename="HomogeneityTest_merra2_%s" % breaktime)


def correct_grid(workdir):
    '''
    Will not be needed for newer versions, for old version last lon was missing due to use of false globgrid
    :param filepath:
    :return:
    '''
    # TODO: Delete this function when everything runs
    fileslist = glob.glob(os.path.join(workdir, "oldHomogeneityTest*.nc"))
    filenames = [afile.replace(workdir + '\\', '') for afile in fileslist]


    grid = globalCellgrid()
    grid_points = grid.get_grid_points()
    glob_lon = np.unique(grid_points[1])
    glob_lat = np.unique(grid_points[2])

    for filename in filenames:

        ncfile = Dataset(os.path.join(workdir,filename), 'r')

        file_lat = np.flipud(ncfile.variables['lat'][:])
        file_lon = ncfile.variables['lon'][:]

        var_names = []
        for var_name in ncfile.variables.keys():
            if ncfile.variables[var_name].dimensions == ('lat','lon'):
                var_names.append(var_name)

        data = {name: np.flipud(ncfile.variables[name][:]) for name in var_names}

        if 'status' in var_names and data['status'].dtype == 'int64':
            data['status'] = data['status'].astype(float)
            data['status'][np.where(data['status']<0)] = np.nan


        if not (glob_lon.size == file_lon.size and all(np.equal(glob_lon, file_lon))):
            for i, lon in enumerate(glob_lon):
                if lon not in file_lon:
                    for var_name in data.keys():
                        data[var_name] = np.insert(data[var_name], i, np.nan, axis=1)

        if not (glob_lat.size == file_lat.size and all(np.equal(glob_lat, file_lat))):
            for i, lat in enumerate(glob_lat):
                if lat not in file_lat:
                    for var_name in data.keys():
                        data[var_name] = np.insert(data[var_name], i, np.nan, axis=0)  #np.full(lons.size, np.nan)


        data_flat ={}
        for name, data in data.iteritems():
            data_flat[name] = data.flatten()

        dataframe = pd.DataFrame(data=data_flat)
        var_meta_dicts ={'status': {}}
        for name in ncfile.variables['status'].ncattrs():
            var_meta_dicts['status'][name] = getattr(ncfile.variables['status'],name)

        breaktime = filename.replace('.','_').split('_')[2]

        points_to_netcdf(dataframe[var_names],workdir,None,'HomogeneityTest_merra2_%s' % breaktime,None,var_meta_dicts)


#correct_grid(r"H:\workspace\HomogeneityTesting\output\CCI33")
#csv_to_netcdf(r"H:\workspace\HomogeneityTesting\output\CCI31EGU")