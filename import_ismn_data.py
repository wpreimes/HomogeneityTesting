# -*- coding: utf-8 -*-
"""
Created on Jun 27 15:33 2017

@author: wpreimes
"""


import numpy as np
import pandas as pd
from rsroot import root_path
import os
import pytesmo.io.ismn.interface as ismn
import pickle
from otherfunctions import merge_ts, regress
from smecv_grid.grid import SMECV_Grid_v042
from grid_functions import cells_for_continent
from rsdata import root_path
from import_satellite_data import QDEGdata_D


class ISMNdataUSA(object):
    def __init__(self, scaleprod, max_depth=0.1,
                 files_path = None):
        # Create a list of gpis nearest to the stations of the dataset
        # If a gpi is nearest for multiple stations,
        # create a list of stations for these gpis that have to be merged
        # when importing data for the gpi
        self.path_ismn_usa = os.path.join(root_path.u, 'datasets', 'ISMN', 'insituUSA',
                                          'Data_seperate_files_19500101_20170321_2365493_xzeO_20170321')

        self.scaleprod = scaleprod
        self.scaledata = QDEGdata_D(products=[scaleprod])
        self.max_depth = max_depth
        self.ISMN_reader = ismn.ISMN_Interface(self.path_ismn_usa)
        self.networks = self.ISMN_reader.list_networks()
        self.grid = SMECV_Grid_v042()
        self.gpis_with_netsta = self._init_assign_gpis_stations()


    def _init_assign_gpis_stations(self):
        # For each station search the nearest gpi
        cells_USA = cells_for_continent('United_States')
        subgrid = self.grid.subgrid_from_cells(cells_USA)
        gpis_with_netsta = {}
        for i, network in enumerate(self.networks):
            print(network, '%i of %i' % (i, len(self.networks) - 1))
            stations = self.ISMN_reader.list_stations(network=network)
            for station in stations:
                station_obj = self.ISMN_reader.get_station(stationname=station, network=network)
                gpi, dist = subgrid.find_nearest_gpi(station_obj.longitude,
                                                     station_obj.latitude)

                variables = station_obj.get_variables()
                if 'soil moisture' in variables:
                    depths_from, depths_to = station_obj.get_depths('soil moisture')

                    # Check if any sensor measured in the correct depth
                    if any(np.around(depths_to, decimals=2) <= self.max_depth):

                        #station_timeframe = station_obj.get_min_max_obs_timestamp()
                        ## Check if station measured during the timeframe

                        #if (station_timeframe[0] < self.timeframe[1]) and \
                        #        (station_timeframe[1] > self.timeframe[0]):

                        if gpi in gpis_with_netsta.keys():
                            gpis_with_netsta[gpi].append((network, station))
                        else:
                            gpis_with_netsta.update({gpi: [(network, station)]})
        return gpis_with_netsta



    def read_all_for_stations_near_gpi(self, gpi):
        # Only bother with reading the station near the gpi if the gpi is in the objects list
        df_stations = pd.DataFrame()
        stations_dict = {}
        for i, (network, station) in enumerate(self.gpis_with_netsta[gpi]):
            stations_dict.update({station: {}})
            station_obj = self.ISMN_reader.get_station(stationname=station,
                                                       network=network)

            depths_from, depths_to = station_obj.get_depths('soil moisture')
            depths_from = np.unique(depths_from)
            depths_to = np.unique(depths_to)
            # Exclude depths below max depth
            depths_to = depths_to[np.where(depths_to < self.max_depth)]
            depths_from = depths_from[np.where(depths_to < self.max_depth)]

            # Iterate over all valid depths for all sensors in this depths
            for j, (depth_from, depth_to) in enumerate(zip(depths_from, depths_to)):
                sensors = station_obj.get_sensors('soil moisture', depth_from, depth_to)
                stations_dict[station].update({(round(depth_from, 3), round(depth_to, 3)): sensors.tolist()})
                for k, sensor in enumerate(sensors):
                    try:
                        data = station_obj.read_variable('soil moisture',
                                                         depth_from=depths_from[0],
                                                         depth_to=depths_to[0],
                                                         sensor=sensor).data
                        df_stations['%s_%s_%s_%s' % (station,
                                                     str(round(depth_from, 3)),
                                                     str(round(depth_to, 3)), sensor)] = \
                            data[data['soil moisture_flag'] == 'G']['soil moisture']
                    except:
                        continue

        return stations_dict, df_stations


    def get_scale_ts(self,gpi,starttime, endtime):
        return self.scaledata.read_gpi(gpi, starttime, endtime)

    def merge_stations_around_gpi(self, gpi, scale_ts):
        '''
        Reads the ISMN data for stations, for which the selected gp is the nearest one,
        merges multiple depths, multiple sensors and multiple stations via linear regression
        fitting (by interpolating missing values in base series (highest correlation)from
        :param gpi:
        :return:
        '''
        stations_dict, df_stations = self.read_all_for_stations_near_gpi(gpi)
        # Include only meansurements +-6H around 0:00H
        # TODO: can measurements from different sensors but same depth be averaged?
        # TODO: can measurements from different depths be averaged?
        df_sensor = pd.DataFrame()
        for s, station_name in enumerate(stations_dict.keys()):
            for d, depth in enumerate(stations_dict[station_name]):
                # For each depth merge TS from available sensors
                for i, sensor in enumerate(stations_dict[station_name][depth]):
                    ts_sensor = df_stations['%s_%s_%s_%s' % (station_name, depth[0], depth[1], sensor)]
                    df_sensor['%s_%i_%i' % (s, d, i)] = ts_sensor
                if df_sensor.empty:
                    continue

                # If a sensore measures daily, just use the daily values as represtative
                if np.unique(df_sensor.index.date).size == df_sensor.index.date.size:
                    df_sensor = df_sensor.resample('24H', base=18, label='right').mean()
                else:
                    df_sensor = df_sensor.between_time('18:00', '06:00', include_start=True, include_end=True)
                df_sensor = df_sensor.resample('12H', base=18, label='right').mean()
                df_sensor = df_sensor.at_time('06:00').shift(-6, freq='H')

        if not scale_ts:
            scale_ts = self.get_scale_ts(gpi,
                                         pd.to_datetime(df_sensor.index.values[0]).to_pydatetime(),
                                         pd.to_datetime(df_sensor.index.values[-1]).to_pydatetime())

        df_merged = pd.DataFrame()
        for column in df_sensor:
            merged = pd.concat([df_sensor[[column]].rename(columns={column:'refdata'}),
                                   scale_ts.rename(columns={self.scaleprod:'testdata'})],axis=1)
            ismn_bias_corr, R, pval, ress = regress(merged.dropna())
            df_merged[column] = ismn_bias_corr

        return df_merged.mean(axis=1, skipna=False)

        '''
                df_depth['Depth%i' % d] = merge_ts(df_sensor)
            if df_depth.empty:
                continue
            df_station[station_name] = merge_ts(df_depth)
        if df_station.empty:
            raise Exception('No valid insitu data around gpi')
        return df_station
        '''
        '''
        for station in df_station.columns.values:
            df_merged = pd.concat([df_cci, df_station[station]], axis=1)
            df_merged = df_merged.dropna()
            df_merged = df_merged.rename(columns={station: 'refdata', df_merged.columns.values[0]: 'testdata'})
            # bivariate linear correlation --> b,c to remove add/mult biases
            df_merged[station], R, pval, ress = regress(df_merged)
            df_station[station + '_bias_corr'] = df_merged[station]

        df_station = df_station[df_station.columns[df_station.columns.to_series().str.contains('_bias_corr')]]
        df_station['testdata'] = df_cci
        weights = df_station.corr()['testdata'] ** 2

        del weights['testdata']
        del df_station['testdata']
        # Normalization: Pos lin correlation betw CCI and each station-->weighting coeff
        sumweights = weights.sum()
        df_station['insitu'] = (df_station * weights).sum(axis=1, skipna=False) / sumweights
        df_station['testdata'] = df_cci

        df_station = df_station.resample('M').mean()
        
        return df_station['insitu'][self.timeframe[0]:self.timeframe[1]]
        '''

    def read_gpi(self, gpi, startdate, enddate, scalets=None):
        data_group = pd.DataFrame(index=pd.date_range(startdate, enddate, freq='D'))
        try:
            ts_ismnmerge = self.merge_stations_around_gpi(gpi, scalets)
            ts_ismnmerge.resample('D').mean()
        except:
            ts_ismnmerge = pd.Series(index=pd.date_range(start=startdate, end=enddate))

        if ts_ismnmerge.isnull().all():
            pass
        ts_ismnmerge.index = ts_ismnmerge.index.to_datetime().date
        ts_ismnmerge.index = pd.to_datetime(ts_ismnmerge.index)
        data_group['ISMN-Merge'] = ts_ismnmerge

        return data_group[startdate:enddate]

if __name__ == '__main__':
    from datetime import datetime
    gpi = 772117
    testprod = 'CCI_41_COMBINED'
    scaleprod = 'merra2'
    timeframe = [datetime(1978, 10, 26), datetime(2016, 12, 31)]
    scale_data = QDEGdata_D(products=[scaleprod])
    data = ISMNdataUSA(scaleprod, 0.1)
    data2 = QDEGdata_D(products=['merra2'])
    df = data.read_gpi(gpi, timeframe[0], timeframe[1])
    df2 = data2.read_gpi(gpi, timeframe[0], timeframe[1])
    print pd.concat([df, df2],axis=1)
