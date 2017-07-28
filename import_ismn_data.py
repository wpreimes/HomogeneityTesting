# -*- coding: utf-8 -*-
"""
Created on Jun 27 15:33 2017

@author: wpreimes
"""


import numpy as np
import pandas as pd

import os
from pygeogrids.netcdf import load_grid
import pytesmo.io.ismn.interface as ismn
import pickle
from otherfunctions import merge_ts, regress


class ISMNdataUSA(object):
    def __init__(self, timeframe, breaktime, max_depth=0.1):
        # Create a list of gpis nearest to the stations of the dataset
        # If a gpi is nearest for multiple stations,
        # create a list of stations for these gpis that have to be merged
        # when importing data for the gpi
        path_ismn_usa = os.path.join('U:\\', 'datasets', 'ISMN', 'insituUSA',
                                     'Data_seperate_files_19500101_20170321_2365493_xzeO_20170321')

        self.breaktime = breaktime
        self.timeframe = timeframe
        self.max_depth = max_depth
        self.path_ismn = path_ismn_usa
        self.ISMN_reader = ismn.ISMN_Interface(self.path_ismn)
        networks = self.ISMN_reader.list_networks()

        defaultfile = r'H:\HomogeneityTesting_data\ismn_files\USA_gpinetsta_%s_%s_%s.pkl' % (
            timeframe[0].strftime('%Y-%m-%d'),
            breaktime.strftime('%Y-%m-%d'),
            timeframe[1].strftime('%Y-%m-%d'))

        land_grid = load_grid(r"R:\Datapool_processed\GLDAS\GLDAS_NOAH025SUBP_3H\ancillary\GLDAS_025_grid.nc")

        if os.path.isfile(defaultfile):
            with open(defaultfile, 'rb') as f:
                self.gpis_with_netsta = pickle.load(f)
        else:
            print('File for stations near GPI not found. Creating...')
            self.gpis_with_netsta = {}
            # IDS of measurements of valid variable and depth

            for i, network in enumerate(networks):
                print(network, '%i of %i' % (i, len(networks) - 1))
                stations = self.ISMN_reader.list_stations(network=network)
                for station in stations:
                    station_obj = self.ISMN_reader.get_station(stationname=station, network=network)
                    gpi, dist = land_grid.find_nearest_gpi(station_obj.longitude,
                                                           station_obj.latitude)

                    variables = station_obj.get_variables()
                    if 'soil moisture' in variables:
                        depths_from, depths_to = station_obj.get_depths('soil moisture')
                        depths_from = np.unique(depths_from)
                        depths_to = np.unique(depths_to)

                        # Check if any sensor measured in the correct depth
                        if any(np.around(depths_to, decimals=2) <= self.max_depth):
                            station_timeframe = station_obj.get_min_max_obs_timestamp()
                            # Check if station measured during the timeframe

                            if (station_timeframe[0] < self.timeframe[1]) and \
                                    (station_timeframe[1] > self.timeframe[0]):

                                if gpi in self.gpis_with_netsta.keys():
                                    self.gpis_with_netsta[gpi].append((network, station))
                                else:
                                    self.gpis_with_netsta.update({gpi: [(network, station)]})

            with open(defaultfile, 'wb') as f:
                pickle.dump(self.gpis_with_netsta, f, pickle.HIGHEST_PROTOCOL)

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

    def merge_stations_around_gpi(self, gpi, df_cci):
        stations_dict, df_stations = self.read_all_for_stations_near_gpi(gpi)
        # Include only meansurements +-6H around 0:00H
        # TODO: can measurements from different sensors but same depth be averaged?
        # TODO: can measurements from different depths be averaged?
        df_station = pd.DataFrame()
        for s, station in enumerate(stations_dict.keys()):
            df_depth = pd.DataFrame()
            for d, depth in enumerate(stations_dict[station].keys()):
                # For each depth merge TS from available sensors
                df_sensor = pd.DataFrame()
                for i, sensor in enumerate(stations_dict[station][depth]):
                    ts_sensor = df_stations['%s_%s_%s_%s' % (station, depth[0], depth[1], sensor)]
                    df_sensor['Sensor%i' % i] = ts_sensor
                if df_sensor.empty:
                    continue

                # If a sensore measures daily, just use the daily values as represtative
                if np.unique(df_sensor.index.date).size == df_sensor.index.date.size:
                    df_sensor = df_sensor.resample('24H', base=18, label='right').mean()
                else:
                    df_sensor = df_sensor.between_time('18:00', '06:00', include_start=True, include_end=True)
                df_sensor = df_sensor.resample('12H', base=18, label='right').mean()
                df_sensor = df_sensor.at_time('06:00')

                df_depth['Depth%i' % d] = merge_ts(df_sensor)
            if df_depth.empty:
                continue
            df_station['Station%i' % s] = merge_ts(df_depth)
        if df_station.empty:
            raise Exception('No valid insitu data around gpi')
        df_station = df_station.shift(-6, freq='H')

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
