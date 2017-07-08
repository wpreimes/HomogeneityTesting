# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 15:11:28 2017

@author: wpreimes
"""

from __future__ import print_function
# Import landgrid

import pygeogrids.netcdf as nc
from pygeogrids.grids import BasicGrid
import os
import rsdata.root_path as root
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from netCDF4 import Dataset, date2num


def get2Dpos(gpis, globalgrid, landgrid):
    grid_points = globalgrid.get_grid_points()
    lats = np.unique(grid_points[2])[::-1]
    lons = np.unique(grid_points[1])
    lon, lat = landgrid.gpi2lonlat(gpis)
    y = datamask(lons, lon)
    x = datamask(lats, lat)
    return x, y


def globalCellgrid():
    # symettrical grid
    lon_res = 0.25
    lat_res = 0.25
    offs_h = lon_res / 2.
    offs_v = lat_res / 2.
    # create meshgrid
    lon, lat = np.meshgrid(
        np.arange(-180 + offs_h, 180 - offs_h, lon_res),
        np.arange(-90 + offs_v, 90 - offs_v, lat_res)
    )

    return BasicGrid(lon.flatten(), lat.flatten()).to_cell_grid(cellsize=5.)


def datamask(x, y):
    index = np.argsort(x)
    sorted_x = x[index]
    sorted_index = np.searchsorted(sorted_x, y)

    yindex = np.take(index, sorted_index, mode="clip")
    mask = x[yindex] != y

    return np.ma.array(yindex, mask=mask)


def create_cellfile_name(gpi, grid):
    # Returns filename (form:cellnumber.nc) and cell for passed gpi in passed grid
    grid_points = grid.get_grid_points()
    gpi_index = np.where(grid_points[0] == gpi)[0][0]
    cell = grid_points[3][gpi_index]
    # Create filename from cell
    file_pattern = str(cell)
    while len(file_pattern) < 4:
        file_pattern = str(0) + file_pattern

    return cell, file_pattern


def update_loc_var(ncfile, data, name, grid, idx):
    if name in ncfile.variables.keys():
        contt = ncfile.variables[name][:]
        if idx != None:
            if contt.ndim == 2:
                x, y = get2Dpos(idx, grid[0], grid[1])
                contt[x, y] = data
            else:
                contt[idx] = data
        else:
            contt = data
        ncfile.variables[name][:] = contt
    else:
        if name in ncfile.dimensions.keys():
            dimension = [ncfile.dimensions[name]]
            dimension_size = dimension[0].size
        elif name == 'location_id':
            dimension = [ncfile.dimensions['locations']]
            dimension_size = dimension[0].size
        elif u'lat' and u'lon' in ncfile.dimensions:
            dimension = [ncfile.dimensions[u'lat'], ncfile.dimensions[u'lon']]
            dimension_size = dimension[0].size * dimension[1].size
        else:
            dimension = [ncfile.dimensions['locations']]
            dimension_size = dimension[0].size
        # If variable does not exist, create it with correct size and retry
        try:
            if isinstance(data, str):
                dtype = str
                contt = np.array([''] * dimension_size, dtype=object)
            else:
                contt = np.full(dimension_size, np.nan)
                if np.asarray(data).dtype == int:
                    dtype = float
                else:
                    dtype = np.asarray(data).dtype
            ncvar = ncfile.createVariable(varname=name,
                                          datatype=dtype,
                                          dimensions=tuple([dim.name for dim in dimension]),
                                          zlib=False)
            ncvar[:] = contt
        except Exception:
            print('Cannot save data for %s to file' % name)

        return update_loc_var(ncfile, data, name, grid, idx)


def points_to_netcdf(dataframe,
                     path,
                     index_col_name=None,
                     filename=None,
                     file_meta_dict=None,
                     var_meta_dicts=None):
    '''
    Write spatial data (data series, data frame) to file:
        -pandas object must contain GPIs as index or in the selected 
        column (index_col_name)
    Parameters
    ----------
    dataframe (mandatory): pandas data frame or data series
        pandas object with data for writing to file
        for time series data: date time as index
        for spatial data: gpi as index
    path (mandatory): string
        path where netcdf file is saved to
    index_col_name (optional): string
        name of the column with time/location data in the pandas object
    filename (optional): string
        for time series data: filename is automatically "*cell*.nc"
        for spatial data: select file name
    file_meta_dict (optional): dictionary
        additional meta information on the netcdf file
    var_meta_dict (optional): dictionary of dictionaries
        additional meta information on the written variables
        for each column in the dataframe, there is 1 dictionary in this list
    overwrite (optional): boolean
        If a (cell)file already exits at the chosen location, existing ground 
        point data is overwritten
    '''

    grid = nc.load_grid(os.path.join(root.r, 'Datapool_processed', 'GLDAS', 'GLDAS_NOAH025_3H.020',
                                     'ancillary', 'GLDASv2_025_land_grid.nc'))

    if not filename:
        filename = 'global'

    # Create or open netcdf cell file
    if os.path.isfile(os.path.join(path, filename + '.nc')):
        ncfile = Dataset(os.path.join(path, filename + '.nc'), "a", format="NETCDF4")
    else:
        ncfile = Dataset(os.path.join(path, filename + '.nc'), "w", format="NETCDF4")
    try:
        globgrid = globalCellgrid()
        grid_points = grid.get_grid_points()
        global_grid_points = globgrid.get_grid_points()

        # TODO: Why -1
        latitudes, longitudes = np.unique(global_grid_points[2])[::-1], np.unique(global_grid_points[1])
        locations = grid_points[0]

        if index_col_name:
            locs = dataframe[index_col_name]
        else:
            locs = dataframe.index
        # glob_pos contains the indices of points to process in the overall grid
        pos = datamask(np.array(locations), np.array(locs))

        n_gpis = locations.size

        # Create data dimensions for Time series and global image
        if not ncfile.dimensions:
            ncfile.createDimension(dimname='locations', size=n_gpis)
            ncfile.createDimension(dimname='lat', size=latitudes.size)
            ncfile.createDimension(dimname='lon', size=longitudes.size)
        # TODO: Add Metadata for netcdf file to dict
        if not ncfile.ncattrs():
            meta_dict = {'geospatial_lon_min': longitudes[0],
                         'geosptial_lon_max': longitudes[-1],
                         'geospatial_lat_min': latitudes[-1],
                         'geospatial_lat_max': latitudes[0],
                         'id': 'global',
                         'date_created': datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
            if file_meta_dict:
                meta_dict.update(file_meta_dict)
            ncfile.setncatts(meta_dict)


            # Create variable for locations and add value
        # GPI, LAT, LON werden beim erstellen immer gefüllt je nach grid unabhängig vom GPI
        # Statt None: gpi_index: Nur für den prozessierten gpi werden idx,lat,lon ins file gespeichert
        meta = {'long_name': 'Location Index', 'standard_name': 'GPI', 'valid_range': '[0 Grid Dependant'}
        update_loc_var(ncfile, locations, u'location_id', grid, pos)
        meta = {'units': 'degrees_east', 'long_name': 'location longitude', 'standard_name': 'longitude',
                'valid_range': '[-180. 180.]'}
        update_loc_var(ncfile, longitudes, u'lon', grid, None)
        ncfile.variables[u'lon'].setncatts(meta)
        meta = {'units': 'degrees_north', 'long_name': 'location latitude', 'standard_name': 'latitude',
                'valid_range': '[-90. 90.]'}
        update_loc_var(ncfile, latitudes, u'lat', grid, None)
        ncfile.variables[u'lat'].setncatts(meta)

        for i, var in enumerate(dataframe.columns.values):
            glob_pos = datamask(global_grid_points[0], locs.values)
            update_loc_var(ncfile, dataframe[var].values, var, [globgrid, grid], glob_pos)
            try:
                ncfile.variables[var].setncatts(var_meta_dicts[var])
            except KeyError:
                ##TODO: Make more useful auto meta data
                var_meta_auto = {'name': var, 'info': 'Automatically generated meta data'}
                ncfile.variables[var].setncatts(var_meta_auto)

    except Exception:
        # TODO: handle the case that no metadata was passed
        # print('Error during filling file %s'%filename)
        pass

    ncfile.close()
