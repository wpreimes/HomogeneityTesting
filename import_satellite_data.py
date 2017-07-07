# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 10:42:39 2016

@author: wpreimes
"""

import numpy as np
import pandas as pd
from smecv.input.common_format import CCIDs
from rsdata.ESA_CCI_SM.interface import ESA_CCI_SM

from general.timedate.julian import caldat

from gldas.interface import GLDASTs
from merra.interface import MERRA2_Ts

from general.time_series.anomaly import calc_anomaly

import os
from datetime import datetime
from pygeogrids.netcdf import load_grid
import pytesmo.io.ismn.interface as ismn
import matplotlib.pyplot as plt
import pickle
from HomogeneityTesting.otherfunctions import merge_ts,regress
import re


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


class QDEGdata_M(object):
    
    def __init__(self,products):
        self.products = products
        self.lkup = pd.read_csv(r"H:\workspace\GPI_lookup\gpi_LUT.csv",index_col=0)

        dayproducts = []
            
        if 'merra2' in self.products:
            if os.path.isdir(r'D:\USERS\wpreimes\datasets\MERRA2_M'):
                print('Found local files for merra2')
                path_merra2 = r'D:\USERS\wpreimes\datasets\MERRA2_M'
            else:
                path_merra2 = r"U:\datasets\MERRA\MERRA\MERRA2_MONTHLY\Timeseries_SM"
            self.merra2 = MERRA2_Ts(path_merra2)


        if 'gldas-merged-from-file' in self.products:
            if os.path.isdir(r"D:\USERS\wpreimes\datasets\GLDAS-NOAH\GLDAS_NOAH_merged_M_at0"):
                print('Found local files for gldas-merged at 0H')
                path_gldas_monthly = r"D:\USERS\wpreimes\datasets\GLDAS-NOAH\GLDAS_NOAH_merged_M_at0"
            else:
                path_gldas_monthly = r"U:\datasets\GLDAS-NOAH\GLDAS_NOAH_merged_M_at0"
            self.gldas_merged = GLDASTs(path_gldas_monthly)
            
        if 'gldas_v2' in self.products:
            dayproducts.append('gldas_v2')
        if 'gldas_v1' in self.products:
            dayproducts.append('gldas_v1')
        if 'gldas-merged' in self.products:
            dayproducts.append('gldas-merged')

        cci_re = re.compile("cci_.+")
        if any([cci_re.match(product) for product in self.products]):
            for cci_product in [cci_re.match(product) for product in self.products]:
                if cci_product:
                    cci_product = cci_product.group()
                    dayproducts.append(cci_product)
            
        if dayproducts:
            self.daydata = QDEGdata_D(products=dayproducts, only_sm=True)

            
    def read_gpi(self,gpi,startdate,enddate):
        
        if not hasattr(self,'daydata'):
            data_group = pd.DataFrame()
        else:
            df_day = self.daydata.read_gpi(gpi,startdate,enddate)

            cci_re = re.compile("cci_.+")

            for cci_prod in [product+'_sm' for product in self.products if cci_re.match(product) != None ]:
                #For HomogeneityTesting make monthly merge only if there are 
                #more than 10 measurements a month
                if cci_prod in df_day.columns.values:
                    df_day = df_day[df_day[cci_prod] != -999999.]
                    ts_resample = df_day[cci_prod].resample('M').mean()
                    ts_count = df_day[cci_prod].resample('M').count()
                    ts_resample.ix[ts_count < 10] = np.nan
                    df_day[cci_prod] = ts_resample


            df_month = df_day.resample('M').mean()
            data_group = df_month
                
        if 'gldas-merged-from-file' in self.products:
            try:
                ts_gldas_merged = self.gldas_merged.read_ts(gpi)[['SoilMoi0_10cm_inst']]
            except:
                ts_gldas_merged = pd.Series(index=pd.date_range(start=startdate, end=enddate))
        
            if (len(np.where(~np.isnan(ts_gldas_merged))[0])==0):
                #raise Exception, 'No GLDAS Data available for GPI %i' %gpi
                print 'No merged gldas data for gpi %0i' % gpi
            else:
                ts_gldas_merged = ts_gldas_merged.rename(columns={'SoilMoi0_10cm_inst':'gldas-merged-from-file'})

                '''          
                ts_gldas_merged_for_mean=pd.DataFrame(index=ts_gldas_merged.index.to_datetime(),
                                                   data=ts_gldas_merged)
                ts_gldas_merged_for_mean=ts_gldas_merged_for_mean.resample('M',how='mean')
                '''
            if data_group.empty:
                data_group = ts_gldas_merged
            else:
                data_group = pd.concat([data_group,ts_gldas_merged], axis=1)
                
            
        if 'merra2' in self.products:
            try:
                gpi_merra = self.lkup.loc[gpi].gpi_merra
                ts_merra2 = self.merra2.read_ts(gpi_merra)['GWETTOP']
                ts_merra2 = ts_merra2.resample('M').mean()
                ts_merra2 = ts_merra2 * 100
            except:
                ts_merra2 = pd.Series(index=pd.date_range(start=startdate, end=enddate))
                ts_merra2 = ts_merra2.resample('M').mean()
                
            if ts_merra2.isnull().all():
                print 'No merra2 data for gpi %i' %gpi

            ts_merra2.index = ts_merra2.index.to_datetime().date
            ts_merra2.index = ts_merra2.index.to_datetime()
            
            if data_group.empty:
                data_group['merra2'] = ts_merra2 
            else:
                data_group = pd.concat([data_group,ts_merra2.rename('merra2')], axis=1)

        return data_group[startdate:enddate]
    

            
class QDEGdata_D(object):
    
    #changed path in "D:\workspace\smecv\grids\ecv.py" to grid file
    #changed path definition in D:\workspace\pygrids\GLDAS_NOAH.py
    #change path in rsdata-GLDAS-interface
    #Requires ECV_CCI_gridv4.nc
    
    def __init__(self,products,resample='mean',only_sm=True):
        
        self.products=products
        self.resample=resample
        self.only_sm=only_sm
        hourproducts=[]

        #TODO: before adding versions: Add cfg file and add data pathes in cfgfile
        cci_versions = ['22','31','33']
        cci_types = ['COMBINED','ACTIVE','PASSIVE']

        cci_res = [re.compile("cci_%s.+" % version) for version in cci_versions]
        self.cci={}
        for cci_re in cci_res:
            for cci_product in [cci_re.match(product) for product in self.products]:
                if cci_product:
                    cci_product = cci_product.group()
                    prefix, version, type = cci_product.split('_')
                    type = type.upper()

                    if prefix != 'cci' or version not in cci_versions or type not in cci_types:
                        raise Exception('cci version or product not known...use format cci_XX_PRODUCT')

                    if os.path.isdir(r'D:\USERS\wpreimes\datasets\CCI_%s_D' % version):
                        print('Found local files for cci_%s data' % version)
                        file_cfg = r"H:\workspace\HomogeneityTesting\local_file_paths\ESA_CCI_SM_v0%s.%s_local.cfg"%(version[0],version[1])
                    else:
                        file_cfg = None
                        print('Try files for cci_%s data on R' % version)
                    cci = ESA_CCI_SM('ESA_CCI_SM_v0%s.%s' % (version[0], version[1]),
                                     parameter=type,
                                     cfg_path=r"H:\workspace\HomogeneityTesting\cci_cfg_local")

                    self.cci.update({cci_product:cci})

            
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
            one or multiple products ('gldas', 'trmm')
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
            


        if self.cci:
            for name, cci_product in self.cci.iteritems():
                try:
                    df_cci = pd.DataFrame(cci_product.read(gpi))
                    df_cci=df_cci.set_index('jd')
                    m,d,y=caldat(df_cci.index.values)
                    df_cci.index=pd.to_datetime(pd.DataFrame(data={'year':y,'month':m,'day':d}))
                    df_cci=df_cci[startdate:enddate]

                    df_cci = df_cci[df_cci['flag'] == 0]
                    df_cci['sm'].loc[df_cci['sm'] == -999999.] = np.nan

                except:
                    df_cci = pd.Series(index=pd.date_range(start=startdate, end=enddate))

                if df_cci['sm'].isnull().all():
                    print 'No cci sm data for gpi %0i' % gpi

                if data_group.empty:
                    if self.only_sm:
                        data_group['%s' % name] = df_cci['sm']
                    else:
                        data_group[name + '_' + df_cci.columns.values] = df_cci
                else:
                    if self.only_sm:
                        data_group = pd.concat([data_group, df_cci['sm'].rename('%s' % name)], axis=1)
                    else:
                        data_group= pd.concat([data_group,df_cci],axis=1)
            
        if 'trmm' in self.products:
            try:
                ts_trmm = self.products['trmm'].read_gp(gpi)['p'][startdate:enddate]
                ts_trmm[np.isnan(ts_trmm)]=0
            except:
                ts_trmm = pd.Series(index=pd.date_range(start=startdate,end=enddate))
            if ts_trmm.isnull().all():
                print 'No trmm data for gpi %0i' % gpi
    
                
            ts_trmm.index=ts_trmm.index.to_datetime().date
            
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
            
                if (len(np.where(~np.isnan(ts_gldas_v2))[0]) == 0) or (len(np.where(~np.isnan(ts_gldas_v21))[0]) == 0):
                    print 'No gldas v2 or v21 data for gpi %0i' % gpi
                    df_gldas_merged = pd.DataFrame(index=pd.date_range(start=startdate, end= enddate, freq='3H'),
                                                   data={'gldas-merged':np.nan})
                else:
                    df_gldas_merged = pd.concat([ts_gldas_v2.rename('gldas_v2'),
                                                 ts_gldas_v21.rename('gldas_v21')], axis=1)

                ts_gldas_merged = df_gldas_merged.mean(axis=1).rename('gldas-merged')
                    

                if data_group.empty:
                    data_group['gldas-merged'] = ts_gldas_merged
                else:
                    data_group = pd.concat([data_group,ts_gldas_merged.rename('gldas-merged')],axis=1)
        
            return data_group[startdate:enddate]

#Testing ISMN
'''
timeframe=['2002-07-01','2011-10-01']

path=os.path.join('U:\\','datasets','ISMN','insituUSA','Data_seperate_files_19500101_20170321_2365493_xzeO_20170321')
ismn_obj=ISMNdata_USA(timeframe=timeframe,breaktime='2007-01-01')

data=QDEGdata_M(products=['cci_33_combined','cci_22_passive'])
ts_cci=data.read_gpi(737687,timeframe[0],timeframe[1])


try:
    DF_Time=ismn_obj.merge_stations_around_gpi(737687,ts_cci/100)  
except:
    print('GPI failed for ISMN data')

ttime=['2007-07-01','2011-10-01','2012-07-01']
test_data=QDEGdata_M(products=['cci_31_active','cci_31_passive','cci_31'])
ts3=test_data.read_gpi(339656,ttime[0],ttime[2])



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

test_data=QDEGdata_M(products=['merra2'])
ts2=test_data.read_gpi(,ttime[0],ttime[2])


ts=pd.concat([ts1,ts2],axis=1)
'''