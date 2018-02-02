# -*- coding: utf-8 -*-
"""
Created on Jun 27 15:33 2017

@author: wpreimes
"""


import numpy as np
import pandas as pd
import os
import pytesmo.io.ismn.interface as ismn
from smecv_grid.grid import SMECV_Grid_v042
from rsdata import root_path
from import_satellite_data import QDEGdata_D
from csv_writer import dict_csv_wrapper
from pytesmo.scaling import add_scaled
from scipy import stats


def wavg(dataframe, weights):
    '''
    Calculate the weighted average of a numpy array
    :param dataframe: Dataframe
        dataframe with multiple time series to merge
    :param weights: np.array
        weights for combining the columns of the dataframe, must be same size as number of columns
    :return:
    '''
    df_averaged = pd.DataFrame(index=dataframe.index)
    for index, line in dataframe.iterrows():
        data = line.values
        if any(np.isnan(data)):
            data = np.ma.masked_invalid(data)
            weights = np.ma.masked_array(weights, data.mask)
        if (data.size == 0) or (weights.size == 0):
            return np.nan

        weighted_mean = np.average(data, weights=weights)
        df_averaged.set_value(index, 'ismn_weighted_average', weighted_mean)

    return df_averaged.iloc[:, 0]


def fill_holes_lress(fill_series, other_series):
    '''
    Perform regression of column refdata on column testdata as used for combining
    multiple insitu measurements
    '''
    # Find holes in max_series where there is data in other_series
    fill_series = fill_series.copy(True)
    other_series = other_series.copy(True)

    df = pd.concat([fill_series, other_series], axis=1).rename(columns={fill_series.name: 'fill_series',
                                                                        other_series.name: 'other_series'})

    # Drop values where both sets are nan, no interpolation possible here
    df = df.dropna(how='all')
    # Group blocks where contiguous values are Nan in the series to fill

    df.loc[:, 'nanmask'] = df['fill_series'].isnull().astype(bool)
    df.loc[:, 'nanmask_shift'] = df['nanmask'].shift(1)
    df.loc[:, 'nan_change'] = df['nanmask'] != df['nanmask_shift']
    df.loc[:, 'nangroup'] = df['nan_change'].cumsum()

    holes = df[['fill_series', 'other_series', 'nangroup', 'nanmask']].groupby(['nangroup', 'nanmask'])
    # calculate regression coefficients for holes from other_series
    for count, hole in holes:
        if hole.nanmask.all():
            slope, intercept, rvalue, pvalue, stderr = stats.linregress(range(len(hole.fill_series)),
                                                                        hole.other_series)

            df.ix[hole.index, 'fill_series'] = intercept + slope * hole.other_series

    return df['fill_series']

    # fill holes in max_series from regression coefficients found




def merge_ts(dataframe_in):
    '''
    Merge temporally coinciding timeseries via interpolation from
    bivariate linear regression.
    Input: Dataframe of temporally coinciding timeseries.
    -Find TS with most values. Use as reference to enhance.
    -Enhance missing values in reference via linear regression from other TS.
        Only use correlating (R>0.8) TS for interpolation.
    -Return gap filled reference TS
    '''
    dataframe = dataframe_in.copy(True)
    max_series_name = np.argmax(dataframe.notnull().sum().to_frame()[0])
    quick_corr = dataframe.corr().sort_values(max_series_name, ascending=False).index.values[1:]
    # Order other TSs via correlation to max_series and return list corr_ordered
    # Improve max_series first with best fitting other,than second best fitting other,...
    DF_other_series = pd.DataFrame()
    weights = np.array([])
    for other_series_name in dataframe[quick_corr]:
        subframe = dataframe[[max_series_name, other_series_name]].copy(True)
        corr, pval = stats.pearsonr(subframe.dropna()[max_series_name],
                                    subframe.dropna()[other_series_name])
        #print corr, pval
        if corr >= 0.8 and pval < 0.01:
            ismn_bias_corr = add_scaled(subframe.dropna(), 'linreg', label_in=other_series_name, label_scale=max_series_name)
            DF_other_series[other_series_name] = ismn_bias_corr[other_series_name + '_scaled_linreg']
            weights = np.append(weights, corr)

    max_series = dataframe.loc[:,max_series_name]
    #other_series_merged = DF_other_series.mean(axis=1)
    if max_series.isnull().all():
        return None
    elif all(DF_other_series.isnull().all()):
        # If there is no second TS, return the unchanged reference series
        return max_series
    else:
        other_series_merged = wavg(DF_other_series, weights)
        # Else perform interpolation from linear regression
        merged_ts = fill_holes_lress(max_series, other_series_merged)
        return merged_ts.rename('ismn_weighted_average')




class ISMNdata(object):
    def __init__(self, scaleprod=None, output_files_path=None, path=None, max_depth=0.1):
        # Create a list of gpis nearest to the stations of the dataset
        # If a gpi is nearest for multiple stations,
        # create a list of stations for these gpis that have to be merged
        # when importing data for the gpi
        if not path:
            self.path_ismn_data = os.path.join(root_path.u, 'datasets', 'ISMN', 'insituUSA',
                                          'Data_seperate_files_19500101_20170321_2365493_xzeO_20170321')
        else:
            self.path_ismn_data = path
        if scaleprod:
            self.scaleprod = scaleprod
            self.scaledata = QDEGdata_D(products=[scaleprod])
        else:
            self.scaleprod = False
            self.scaledata = False

        self.max_depth = max_depth
        self.ISMN_reader = ismn.ISMN_Interface(self.path_ismn_data)
        self.networks = self.ISMN_reader.list_networks()
        self.grid = SMECV_Grid_v042()
        self.output_path = output_files_path
        self.df_gpis = pd.DataFrame.from_dict(self.assign_netsta_to_gpis()).set_index('gpi')



    def assign_netsta_to_gpis(self):
        # For each station search the nearest gpi
        gpis_with_netsta = {'cell':[], 'gpi':[], 'network':[], 'station':[], 'dist':[],
                            'station_lat':[], 'station_lon':[]}
        for i, network in enumerate(self.networks):
            print(network, '%i of %i' % (i, len(self.networks) - 1))
            stations = self.ISMN_reader.list_stations(network=network)
            for station in stations:
                station_obj = self.ISMN_reader.get_station(stationname=station, network=network)
                lon = station_obj.longitude
                lat = station_obj.latitude
                gpi, dist = self.grid.find_nearest_gpi(lon,lat)

                variables = station_obj.get_variables()
                if 'soil moisture' in variables:
                    depths_from, depths_to = station_obj.get_depths('soil moisture')
                    # Check if any sensor measured in the correct depth
                    if any(np.around(depths_to, decimals=2) <= self.max_depth):
                        gpis_with_netsta['gpi'].append(gpi)
                        gpis_with_netsta['cell'].append(self.grid.gpi2cell(gpi))
                        gpis_with_netsta['network'].append(network)
                        gpis_with_netsta['station'].append(station)
                        gpis_with_netsta['dist'].append(round(dist,3))
                        gpis_with_netsta['station_lat'].append(round(lat,3))
                        gpis_with_netsta['station_lon'].append(round(lon,3))

        writer = dict_csv_wrapper(gpis_with_netsta)
        if self.output_path:
            writer.write(os.path.join(self.output_path, 'gpis_netsta.csv'))

        return writer.content

    def read_all_for_stations_near_gpi(self, df_for_gpi):
        # Only bother with reading the station near the gpi if the gpi is in the objects list
        df_stations = pd.DataFrame()
        stations_dict = {}
        for i, (network, station) in enumerate(zip(df_for_gpi['network'].values,
                                                   df_for_gpi['station'].values)):

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
        if self.scaledata:
            return self.scaledata.read_gpi(gpi, starttime, endtime) / 100
        else:
            return None


    def merge_stations_around_gpi(self, gpi, scale_ts):
        '''
        Reads the ISMN data for stations, for which the selected gp is the nearest one,
        merges multiple depths, multiple sensors and multiple stations via linear regression
        fitting (by interpolating missing values in base series (highest correlation)from
        :param gpi:
        :return:
        '''
        df_for_gpi = self.df_gpis.loc[gpi]
        stations_dict, df_stations = self.read_all_for_stations_near_gpi(df_for_gpi)
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

        if not scale_ts and self.scaleprod:
            scale_ts = self.get_scale_ts(gpi,
                                         pd.to_datetime(df_sensor.index.values[0]).to_pydatetime(),
                                         pd.to_datetime(df_sensor.index.values[-1]).to_pydatetime())

        if self.scaleprod:
            print('Merge with scale product %s' %self.scaleprod)
            #Fit all single ISMN Series to the common scale ts and calculate weighted mean from correlations
            df_merged = pd.DataFrame()
            df_weights = pd.DataFrame(index = ['weights'])
            for column in df_sensor:
                merged = pd.concat([df_sensor[[column]].rename(columns={column:'ismn'}),
                                       scale_ts.rename(columns={self.scaleprod:'reference'})],axis=1)

                #ismn_bias_corr, R, pval, ress = regress(merged.dropna())
                ismn_bias_corr = add_scaled(merged.dropna(), 'linreg', label_in='ismn', label_scale='reference')
                df_merged[column] = ismn_bias_corr['ismn_scaled_linreg']
                corr, p = stats.pearsonr(ismn_bias_corr['ismn_scaled_linreg'].values,
                                         ismn_bias_corr['reference'].values)
                corr = ismn_bias_corr.corr().loc['ismn_scaled_linreg','reference']
                if corr >0.8 and p < 0.01: #TODO: Change to 0.8 and 0.01
                    df_weights.set_value('weights', column, corr**2)
                else:
                    df_weights.set_value('weights', column, np.nan)

            df_averaged = wavg(df_merged, df_weights.loc['weights'].values)

            return df_averaged.rename('ismn_weighted_average')
        else:
            print('Merge wirh longst sensor series')
            # Use the longest sensor TS as reference and fill its holes from other, correlating sensors
            return merge_ts(df_sensor)



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
        data_group['ISMN-Merge'] = ts_ismnmerge * 100

        return data_group[startdate:enddate]

if __name__ == '__main__':
    from datetime import datetime

    ismn_data_path = os.path.join(root_path.u, 'datasets', 'ISMN', 'insituUSA',
                 'Data_seperate_files_19500101_20170321_2365493_xzeO_20170321')

    gpi = 700119
    timeframe = [datetime(1978, 10, 26), datetime(2016, 12, 31)]
    data = ISMNdata(None, r'C:\Temp\csv',ismn_data_path, max_depth=0.1)
    data2 = ISMNdata('CCI_41_COMBINED', r'C:\Temp\csv',ismn_data_path, max_depth=0.1)
    df = data.read_gpi(gpi, timeframe[0], timeframe[1]).rename(columns={'ISMN-Merge':'ISMN_scaleself'})
    df2 = data2.read_gpi(gpi, timeframe[0], timeframe[1]).rename(columns={'ISMN-Merge':'ISMN_scalets'})
    merged = pd.concat([df2, df], axis=1)
    merged = add_scaled(merged.dropna(), 'linreg', label_in='ISMN_scalets', label_scale='ISMN_scaleself')


