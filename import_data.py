# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 10:42:39 2016

@author: wpreimes
"""

import numpy as np
import pandas as pd
from smecv.input.common_format import CCIDs
from gldas.interface import GLDASTs
from rsdata.TRMM_TMPA.interface import Tmpa3B42_GenericIO_Ts
from merra.interface import MERRA2_Ts

from general.time_series.anomaly import calc_anomaly

import os
from datetime import datetime
from pygeogrids.netcdf import load_grid
import pytesmo.io.ismn.interface as ismn
import matplotlib.pyplot as plt
import pickle
from HomogeneityTesting.otherfunctions import merge_ts
def read_warp_ssm(ssm,ssf,gpi):
    
    import pytesmo.timedate.julian as julian
    
    ssm_ts = ssm.read_ts(gpi)
    ssf_ts = ssf.read_ts(gpi)
    
    # ceil or round?
    jd = np.round(ssm_ts['jd']-0.5)+0.5
    #jd = ssm_ts['jd']
    ts = pd.DataFrame(ssm_ts, index=julian.julian2datetimeindex(jd))
    ts['ssf'] = ssf_ts['ssf']
    
    ts = ts[(ts['proc_flag']<=2)&(ts['ssf']==1)]['sm']
    ts.index.tz = None
    
    return ts.groupby(level=0).last()

class ISMNdata_USA(object):
    def __init__(self,path_ismn,breaktime,timeframe,max_depth=0.1):
        #Create a list of gpis nearest to the stations of the dataset
        #If a gpi is nearest for multiple stations,
            #create a list of stations for these gpis that have to be merged 
            #when importing data for the gpi
        defaultfile=r"H:\workspace\HomogeneityTesting\output\ismn_usa_gpi_netsta.pkl"
        self.breaktime=datetime.strptime(breaktime, '%Y-%m-%d')
        self.timeframe=[datetime.strptime(timeframe[0], '%Y-%m-%d'),
                        datetime.strptime(timeframe[1], '%Y-%m-%d')]
        self.max_depth=max_depth
        self.path_ismn=path_ismn
        self.ISMN_reader = ismn.ISMN_Interface(path_ismn)
        networks = self.ISMN_reader.list_networks()
        
        land_grid=load_grid(r"R:\Datapool_processed\GLDAS\GLDAS_NOAH025SUBP_3H\ancillary\GLDAS_025_grid.nc")
        
        
        
        if os.path.isfile(defaultfile):
            with open(defaultfile, 'rb') as f:
                self.gpis_with_netsta=pickle.load(f)
        else:
            self.gpis_with_netsta={}
            #IDS of measurements of valid variable and depth
                                                 
            for i,network in enumerate(networks):
                print(network,'%i of %i'%(i,len(networks)-1))
                stations = self.ISMN_reader.list_stations(network = network)
                for station in stations:
                    station_obj = self.ISMN_reader.get_station(stationname=station,network=network)
                    gpi,dist=land_grid.find_nearest_gpi(station_obj.latitude,
                                                        station_obj.latitude)
                    
                    variables=station_obj.get_variables()
                    if 'soil moisture' in variables:
                        depths_from,depths_to = station_obj.get_depths('soil moisture')
                        depths_from=np.unique(depths_from)
                        depths_to=np.unique(depths_to)
                        
                        #Check if any sensor measured in the correct depth
                        if any(depths_to <= self.max_depth):
                            station_timeframe=station_obj.get_min_max_obs_timestamp()
                            #Check if station measured during the timeframe
                            
                            if (station_timeframe[0]<self.timeframe[1]) and\
                               (station_timeframe[1]>self.timeframe[0]):
                              
                                if gpi in self.gpis_with_netsta.keys():
                                   self.gpis_with_netsta[gpi].append((network,station))
                                else:
                                   self.gpis_with_netsta.update({gpi:[(network,station)]})
            
            with open(defaultfile, 'wb') as f:
                pickle.dump(self.gpis_with_netsta,f,pickle.HIGHEST_PROTOCOL) 
            
                    
    def read_all_for_stations_near_gpi(self,gpi):
        #Only bother with reading the station near the gpi if the gpi is in the objects list 
        if gpi in self.gpis_with_netsta.keys():
            DF_Stations=pd.DataFrame()
            stations_dict={}
            for i,(network,station) in enumerate(self.gpis_with_netsta[gpi]):
                stations_dict.update({station:{}})
                station_obj=self.ISMN_reader.get_station(stationname=station,
                                                         network=network)
                
                depths_from,depths_to = station_obj.get_depths('soil moisture')
                depths_from=np.unique(depths_from)
                depths_to=np.unique(depths_to)
                #Exclude depths below max depth
                depths_to=depths_to[np.where(depths_to<self.max_depth)]
                depths_from=depths_from[np.where(depths_to<self.max_depth)]
                
                #Iterate over all valid depths for all sensors in this depths
                for j,(depth_from,depth_to) in enumerate(zip(depths_from,depths_to)):
                    sensors=station_obj.get_sensors('soil moisture',depth_from,depth_to)                                     
                    stations_dict[station].update({(round(depth_from,3),round(depth_to,3)):sensors.tolist()})
                    for k,sensor in enumerate(sensors):
                        try:                                  
                            data = station_obj.read_variable('soil moisture',
                                                             depth_from=depths_from[0],
                                                             depth_to=depths_to[0],
                                                             sensor=sensor).data
                            DF_Stations['%s_%s_%s_%s'%(station,
                                                       str(round(depth_from,3)),
                                                       str(round(depth_to,3)),sensor)]=data[data['soil moisture_flag']=='G']['soil moisture']
                        except:
                            continue

            return stations_dict,DF_Stations[self.timeframe[0]:self.timeframe[1]]
        else:
            raise Exception, 'No insitu stations for gpi in valid range'
                
    def merge_station_series(self,gpi):
        stations_dict,DF_Stations=self.read_all_for_stations_near_gpi(gpi)
        #Include only meansurements +-6H around 0:00H
        #TODO: can measurements from different sensorsbut same depth be averaged?
        #TODO: can measurements from different depths be averaged?
        DF_Station=pd.DataFrame()
        for s,station in enumerate(stations_dict.keys()):
            DF_Depth=pd.DataFrame()
            for d,depth in enumerate(stations_dict[station].keys()):
                #For each depth merge TS from available sensors                
                DF_Sensor=pd.DataFrame()
                for i,sensor in enumerate(stations_dict[station][depth]):
                    ts_sensor=DF_Stations['%s_%s_%s_%s'%(station,depth[0],depth[1],sensor)]
                    DF_Sensor['Sensor%i'%i]=ts_sensor
                if DF_Sensor.empty:continue
                DF_Sensor=DF_Sensor.between_time('18:00','06:00',include_start=True,include_end=True)
                DF_Sensor=DF_Sensor.resample('12H',base=18,label='right').mean()                  
                DF_Sensor=DF_Sensor.at_time('06:00')
                
                DF_Depth['Depth%i'%d]=merge_ts(DF_Sensor)
            if DF_Depth.empty:continue
            DF_Station['Station%i'%s]=merge_ts(DF_Depth)                
                #Fill Nan Values in max_series with values from lin regression 
        
        
        
        #iterate over stations
            #iterate over depths
                #Sensor merge:
                #Find ts with most values: This is reference to enhance
        
        
            #Depth merge:
        
        
        #Station merge:
        
        
        
        #Choose longest series as reference series to enhance:
        
        #Merge the measurements of sensors for same station and same depth
        #Merge measurements of    
        DF_Station.plot()
        plt.show()
        
        DF_Stations['merged']=DF_Stations.mean(axis=1)
        DF_Stations['in_timeframe']=DF_Stations['merged'].between_time('18:00','06:00',include_start=True,include_end=True)
        DF_Stations['timeframe_mean']=DF_Stations['in_timeframe'].resample('12H',base=18,label='right').mean()                
        
        #TODO: Mean is influenced by not flagged outliers...
        DF_Time['P:%i'%gpi]=DF_Stations.iloc[DF_Stations.index.hour==6]['timeframe_mean']                          
        #Create one value from +-6h around 0H
        '''
        for point in DF_Time.columns.values
            corr,pval=stats.pearsonsr(DF_Time[point],test_ts)
            if corr >=0.8 and pval<0.01:
                
        DF_Time.index=DF_Time.index.date
        
        return DF_Time
        else:
            raise Exception, 'No ISMN station for GPI cell %i'%gpi
        '''
                    
                    
                                          

                    
path=os.path.join('U:\\','datasets','ISMN','insituUSA','Data_seperate_files_19500101_20170321_2365493_xzeO_20170321')
ismn_obj=ISMNdata_USA(path,breaktime='2007-01-01',timeframe=['2002-07-01','2011-10-01'])
DF_Time=ismn_obj.merge_station_series(722304)  
  
#%%

     
class QDEGdata_M(object):
    
    def __init__(self,products):
        self.products=products
        self.lkup=pd.read_csv(r"H:\workspace\GPI_lookup\gpi_LUT.csv",index_col=0)
        
        dayproducts=[]
        
        if 'cci_22_from_file' in products:
            if os.path.isdir(r'D:\USERS\wpreimes\datasets\CCI_22_M\060_monthlyAverage\combinedProd'):
                print('Found local files for cci22 monthly data')
                path_cci_22_monthly=r'D:\USERS\wpreimes\datasets\CCI_22_M\060_monthlyAverage\combinedProd'
            else:
                path_cci_22_monthly=r"R:\Datapool_processed\ESA_CCI_SM\ESA_CCI_SM_v02.2\060_monthlyAverage\combinedProd"
            self.cci_22 = CCIDs(path_cci_22_monthly)
            
        #There is NO MONTHLY CCI31 Product
            
        if 'merra2' in products:
            if os.path.isdir(r'D:\USERS\wpreimes\datasets\MERRA2_M'):
                print('Found local files for merra2')
                path_merra2=r'D:\USERS\wpreimes\datasets\MERRA2_M'
            else:
                path_merra2 = r"U:\datasets\MERRA\MERRA\MERRA2_MONTHLY\Timeseries_SM"
            self.merra2=MERRA2_Ts(path_merra2)
            
            
        if 'gldas-merged-from-file' in products:
            if os.path.isdir(r"D:\USERS\wpreimes\datasets\GLDAS-NOAH\GLDAS_NOAH_merged_M_at0"):
                print('Found local files for gldas-merged at 0H')
                path_gldas_monthly=r"D:\USERS\wpreimes\datasets\GLDAS-NOAH\GLDAS_NOAH_merged_M_at0"
            else:
                path_gldas_monthly=r"U:\datasets\GLDAS-NOAH\GLDAS_NOAH_merged_M_at0"
            self.gldas_merged=GLDASTs(path_gldas_monthly)
            
        if 'gldas_v2' in products:
            dayproducts.append('gldas_v2')
        if 'gldas_v1' in products:
            dayproducts.append('gldas_v1')
        if 'gldas-merged' in products:
            dayproducts.append('gldas-merged')
        if 'cci_22' in products:
            dayproducts.append('cci_22')
        if 'cci_31' in products:
            dayproducts.append('cci_31')
            
        if dayproducts:
            self.daydata=QDEGdata_D(products=dayproducts)
            
            
            
            
    def read_gpi(self,gpi,startdate,enddate):
        
        if not hasattr(self,'daydata'):
            data_group= pd.DataFrame()
        else:
            df_day=self.daydata.read_gpi(gpi,startdate,enddate)
            for cci_prod in ['cci_22','cci_31']:
                #For HomogeneityTesting make monthly merge only if there are 
                #more than 10 measurements a month
                if cci_prod in df_day.columns.values:
                    ts_resample=df_day[cci_prod].resample('M').mean()
                    ts_count=df_day[cci_prod].resample('M').count()
                    ts_resample.ix[ts_count < 10] = np.nan
                    df_day[cci_prod]=ts_resample
            df_month=df_day.resample('M').mean()
            data_group=df_month
            
        
                
        if 'gldas-merged-from-file' in self.products:
            try:
                ts_gldas_merged = self.gldas_merged.read_ts(gpi)[['SoilMoi0_10cm_inst']]
            except:
                ts_gldas_merged = pd.Series(index=pd.date_range(start=startdate,end=enddate))
        
            if (len(np.where(~np.isnan(ts_gldas_merged))[0])==0):
                #raise Exception, 'No GLDAS Data available for GPI %i' %gpi
                print 'No merged gldas data for gpi %0i' % gpi
            else:
                ts_gldas_merged=ts_gldas_merged.rename(columns={'SoilMoi0_10cm_inst':'gldas-merged-from-file'})

                '''          
                ts_gldas_merged_for_mean=pd.DataFrame(index=ts_gldas_merged.index.to_datetime(),
                                                   data=ts_gldas_merged)
                ts_gldas_merged_for_mean=ts_gldas_merged_for_mean.resample('M',how='mean')
                '''
            if data_group.empty:
                data_group=ts_gldas_merged
            else:
                data_group=pd.concat([data_group,ts_gldas_merged],axis=1)
                
            
        if 'merra2' in self.products:
            try:
                gpi_merra=self.lkup.loc[gpi].gpi_merra
                ts_merra2 = self.merra2.read_ts(gpi_merra)['GWETTOP']
                ts_merra2=ts_merra2.resample('M').mean()
                ts_merra2=ts_merra2*100
            except:
                ts_merra2 = pd.Series(index=pd.date_range(start=startdate,end=enddate))
                ts_merra2=ts_merra2.resample('M').mean()
                
            if ts_merra2.isnull().all():
                print 'No merra2 data for gpi %i' %gpi
  

            ts_merra2.index=ts_merra2.index.to_datetime().date
            ts_merra2.index=ts_merra2.index.to_datetime()
            
            if data_group.empty:
                data_group['merra2'] = ts_merra2 
            else:
                data_group=pd.concat([data_group,ts_merra2.rename('merra2')],axis=1)  
        

        return data_group[startdate:enddate]
    
 #%%           
            
class QDEGdata_D(object):
    
    #changed path in "D:\workspace\smecv\grids\ecv.py" to grid file
    #changed path definition in D:\workspace\pygrids\GLDAS_NOAH.py
    #change path in rsdata-GLDAS-interface
    #Requires ECV_CCI_gridv4.nc
    
    def __init__(self,products,resample='mean'):
        
        self.products=products
        self.resample=resample
        
        hourproducts=[]
        
        
        if 'cci_22' in products:
            if os.path.isdir(r'D:\USERS\wpreimes\datasets\CCI_22_D'):
                print('Found local files for cci22 daily data')
                path_cci_22=r'D:\USERS\wpreimes\datasets\CCI_22_D'
            else:
                path_cci_22 = r"R:\Datapool_processed\ESA_CCI_SM\ESA_CCI_SM_v02.2\050_combinedProduct"
            self.cci_22 = CCIDs(path_cci_22)
            
        if 'cci_31' in products:
            if os.path.isdir(r'D:\USERS\wpreimes\datasets\CCI_31_D'):
                print('Found local files for cci31 daily data')
                path_cci_31=r'D:\USERS\wpreimes\datasets\CCI_31_D'
            else:
                path_cci_31=r"R:\Datapool_processed\ESA_CCI_SM\ESA_CCI_SM_v03.1\050_combinedProduct"
            self.cci_31 = CCIDs(path_cci_31)
            
        if 'ascat' in products:
            if os.path.isdir(r'D:\USERS\wpreimes\datasets\ascata'):
                print('Found local files for ascat daily data')
                path_ascat=r'D:\USERS\wpreimes\datasets\ascata'
            else:
                path_ascat = r"U:\datasets\ascata"
            self.ascat=CCIDs(path_ascat)
            
        if 'gldas_v1' in products:
            hourproducts.append('gldas_v1')
        if 'gldas_v2' in products:
            hourproducts.append('gldas_v2')
        if 'gldas_v21' in products:
            hourproducts.append('gldas_v21')
        if 'gldas-merged' in products:
            hourproducts.append('gldas-merged')
        
        if hourproducts:
            self.hourdata=QDEGdata_3H(products=hourproducts)
            
    def read_gpi(self,gpi,startdate,enddate):
        """
        Read one or multiple products for selected ground point in QDEG grid
        and return SM values as dataframe from a defined start date to a
        defined end date.
        
        Parameters
        ----------
        gpi : int
            index of point in QDEG grid
            
        startdate,enddate : strings of form: 'yyyy-mm-dd'
            start and end date of sm timeseries
            
        products : *list
            one or multiple products ('ascat', 'amsre', 'amsr2', 'gldas', 'trmm')
        """
        
        
        if not hasattr(self,'hourdata'):
            data_group= pd.DataFrame(index=pd.date_range(startdate,enddate,freq='D'))
        else:
            df_hour=self.hourdata.read_gpi(gpi,startdate,enddate)
            if self.resample == 'mean':
                #Resample over all values of a day
                df_hour=df_hour.resample('D').mean()
            else:
                #Resample only over the vales of a certain time for each day
                df_hour=df_hour.at_time(self.resample)
                df_hour=df_hour.resample('D').mean()

            data_group=df_hour
            
            

        if 'ascat' in self.products:
            try:
                ts_ascat = self.ascat.read_ts(gpi)
                ts_ascat = ts_ascat[ts_ascat.flag==0]['sm'][startdate:enddate]
            except:
                ts_ascat = pd.Series(index=pd.date_range(start=startdate,end=enddate))
            
            if ts_ascat.isnull().all():
                print 'No cci data for gpi %0i' % gpi
            
            ts_ascat.index=ts_ascat.index.to_datetime().date
            ts_ascat.index=ts_ascat.index.to_datetime()
            if data_group.empty:
                data_group['ascat'] = ts_ascat
            else:
                data_group=pd.concat([data_group,ts_ascat.rename('ascat')],axis=1)
            
               
        if 'cci_22' in self.products:
            try:
                ts_cci = self.cci_22.read_ts(gpi,mask_sm_nan=True)['sm'][startdate:enddate]
            except:
                ts_cci = pd.Series(index=pd.date_range(start=startdate,end=enddate))
                
            if ts_cci.isnull().all():
                print 'No cci data for gpi %0i' % gpi
                
            ts_cci.index=ts_cci.index.to_datetime().date
            ts_cci.index=ts_cci.index.to_datetime()
            if data_group.empty:
                data_group['cci_22'] = ts_cci
            else:
                data_group=pd.concat([data_group,ts_cci.rename('cci_22')],axis=1)
            
            
        if 'cci_31' in self.products:
            try:
                ts_cci = self.cci_31.read_ts(gpi,mask_sm_nan=True)['sm'][startdate:enddate]
            except:
                ts_cci = pd.Series(index=pd.date_range(start=startdate,end=enddate))
            if ts_cci.isnull().all():
                print 'No cci data for gpi %0i' % gpi

            ts_cci.index=ts_cci.index.to_datetime().date
            ts_cci.index=ts_cci.index.to_datetime()

            if data_group.empty:
                data_group['cci_31']= ts_cci
            else:
                data_group=pd.concat([data_group,ts_cci.rename('cci_31')],axis=1)
                    
                
            
        if 'trmm' in self.products:
            try:
                ts_trmm = self.trmm.read_gp(gpi)['p'][startdate:enddate]
                ts_trmm[np.isnan(ts_trmm)]=0
            except:
                ts_trmm = pd.Series(index=pd.date_range(start=startdate,end=enddate))
            if ts_cci.isnull().all():
                print 'No trmm data for gpi %0i' % gpi
    
                
            ts_trmm.index=ts_trmm.index.to_datetime().date
            ts_trmm.index=ts_cci.index.to_datetime()
            
            if data_group.empty:
                data_group['trmm']= ts_trmm
            else:
                data_group=pd.concat([data_group,ts_trmm.rename('trmm')],axis=1)
        
        return data_group[startdate:enddate]
       
       
       
class QDEGdata_3H(object):

    
    def __init__(self,products):
        self.products=products
        
        if 'gldas_v1' in products:
            if os.path.isdir(r'D:\USERS\wpreimes\datasets\gldas_v1'):
                print('Found local files for GLDAS1 3H data')
                path_gldasv1=r'D:\USERS\wpreimes\datasets\gldas_v1'
            else:
                path_gldasv1=r"R:\Datapool_processed\GLDAS\GLDAS_NOAH025SUBP_3H\datasets\netcdf"
            self.gldas_v1 = GLDASTs(path_gldasv1)
            
        if 'gldas_v2' or 'gldas-merged' in products:
            if os.path.isdir(r'D:\USERS\wpreimes\datasets\gldas_v2'):
                print('Found local files for GLDAS2 3H data')
                path_gldasv2=r'D:\USERS\wpreimes\datasets\gldas_v2'
            else:
                path_gldasv2=r"R:\Datapool_processed\GLDAS\GLDAS_NOAH025_3H.020\datasets\netcdf"
            self.gldas_v2 = GLDASTs(path_gldasv2)
            
        if 'gldas_v21' or 'gldas-merged'in products:
            if os.path.isdir(r'D:\USERS\wpreimes\datasets\gldas_v21'):
                print('Found local files for GLDAS21 3H data')
                path_gldasv21=r'D:\USERS\wpreimes\datasets\gldas_v21'
            else:
                path_gldasv21=r"R:\Datapool_processed\GLDAS\GLDAS_NOAH025_3H.2.1\datasets"
            self.gldas_v21 = GLDASTs(path_gldasv21)
            
 
        
    def read_gpi(self,gpi,startdate,enddate):
            
            data_group= pd.DataFrame()
        
            if 'gldas_v1' in self.products:
                try:
                    ts_gldas_v1 = self.gldas_v1.read_ts(gpi)['086_L1']
                except:
                    ts_gldas_v1 = pd.Series(index=pd.date_range(start=startdate,end=enddate))
            
                if (len(np.where(~np.isnan(ts_gldas_v1))[0])==0):
                    #raise Exception, 'No GLDAS Data available for GPI %i' %gpi
                    print 'No gldas v1 data for gpi %0i' % gpi
                    
                if data_group.empty:
                    data_group['gldas_v1']=ts_gldas_v1
                else:
                    data_group=pd.concat([data_group,ts_gldas_v1.rename('gldas_v1')],axis=1)
                    
                    
            if 'gldas_v2' in self.products:
                try:
                    ts_gldas_v2 = self.gldas_v2.read_ts(gpi)['086_L1']
                except:
                    ts_gldas_v2 = pd.Series(index=pd.date_range(start=startdate,end=enddate))
            
                if (len(np.where(~np.isnan(ts_gldas_v2))[0])==0):
                    #raise Exception, 'No GLDAS Data available for GPI %i' %gpi
                    print 'No gldas v1 data for gpi %0i' % gpi

                if data_group.empty:
                    data_group['gldas_v2']=ts_gldas_v2
                else:
                    data_group=pd.concat([data_group,ts_gldas_v2.rename('gldas_v2')],axis=1)
                    
            
            if 'gldas_v21' in self.products:
                try:
                    ts_gldas_v21 = self.gldas_v21.read_ts(gpi)['SoilMoi0_10cm_inst']
                except:
                    ts_gldas_v21 = pd.Series(index=pd.date_range(start=startdate,end=enddate))
            
                if (len(np.where(~np.isnan(ts_gldas_v21))[0])==0):
                    #raise Exception, 'No GLDAS Data available for GPI %i' %gpi
                    print 'No gldas v1 data for gpi %0i' % gpi

                if data_group.empty:
                    data_group['gldas_21']=ts_gldas_v21
                else:
                    data_group=pd.concat([data_group,ts_gldas_v21.rename('gldas_21')],axis=1)
          
          
            if 'gldas-merged' in self.products:
                try:
                    ts_gldas_v2 = self.gldas_v2.read_ts(gpi)['086_L1']
                    ts_gldas_v21 = self.gldas_v21.read_ts(gpi)['SoilMoi0_10cm_inst']
                except:
                    ts_gldas_v2 = pd.Series(index=pd.date_range(start=startdate,end=enddate))
                    ts_gldas_v21 = pd.Series(index=pd.date_range(start=startdate,end=enddate))
            
                if (len(np.where(~np.isnan(ts_gldas_v2))[0])==0) or(len(np.where(~np.isnan(ts_gldas_v21))[0])==0):
                    #raise Exception, 'No GLDAS Data available for GPI %i' %gpi
                    print 'No gldas v2 or v21 data for gpi %0i' % gpi
                else:
                    df_gldas_merged=pd.concat([ts_gldas_v2.rename('gldas_v2'),
                                               ts_gldas_v21.rename('gldas_v21')],axis=1)

                    ts_gldas_merged=df_gldas_merged.mean(axis=1).rename('gldas-merged')
                    

                if data_group.empty:
                    data_group['gldas-merged']=ts_gldas_merged
                else:
                    data_group=pd.concat([data_group,ts_gldas_merged.rename('gldas-merged')],axis=1)
        
            return data_group[startdate:enddate]
#Testing
'''
ttime=['2002-07-01','2007-01-01','2015-10-01']

#test_data=QDEGdata_D(products=['gldas_v21','gldas_merged','cci_31','cci_22'],resample='0:00')
test_data=QDEGdata_M(products=['gldas-merged-from-file','merra2','cci_31','cci_22'])
ts3=test_data.read_gpi(497651,ttime[0],ttime[2])
'''
'''
'''

'''
#test_data=QDEGdata_M(products=['cci_22','cci_31','gldas_v2'])
#gpi_inh: 363662,406520
test_data=QDEGdata_M(products=['merra2','cci_22'])
ts1=test_data.read_gpi(741898,ttime[0],ttime[2])

ax=ts1.plot(title=r'Example Inhomogeneity')
ax.set_xlabel(r'YEAR')
ax.set_ylabel(r'SM $[m^3/m^3]$')

plt.savefig(r'C:\Users\wpreimes\Desktop\ex_inh1.png', dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)
'''
'''
test_data=QDEGdata_M(products=['merra2'])
ts2=test_data.read_gpi(,ttime[0],ttime[2])


ts=pd.concat([ts1,ts2],axis=1)
'''