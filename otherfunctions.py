# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 12:58:31 2016

@author: wpreimes
"""

from scipy import stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime,timedelta

def regress(dataframe):
    '''
    Perform regression of column refdata on column testdata
    '''
    #Berechnet Korrelation zwischen Testdaten und Referenzdaten
    R,pval=stats.pearsonr(dataframe.refdata,dataframe.testdata)
    out=dataframe.refdata
    if R<0 or np.isnan(R):
        ress=[np.nan]
        return out, R, pval, ress
    
    testdata=dataframe.testdata.values
    refdata=dataframe.refdata.values
    refdata_ones = np.vstack([refdata, np.ones(len(refdata))]).T
 
    ress=np.linalg.lstsq(refdata_ones,testdata)[0][::-1]  
    dataframe['ones']=1
    xm=np.matrix(dataframe.as_matrix(columns=['ones','refdata']))
    out=np.dot(xm,np.matrix(ress).transpose())

    return out, R, pval, ress


def datetime2matlabdn(dt):
   ord = dt.toordinal()
   mdn = dt + timedelta(days = 366)
   frac = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
   return mdn.toordinal() + frac    
    
          
def fill_holes_lress(fill_series,other_series):
    '''
    Perform regression of column refdata on column testdata as used for combining
    multiple insitu measurements
    '''
    #Find holes in max_series where there is data in other_series
    df=pd.DataFrame(data={'fill_series':fill_series,
                          'other_series':other_series},
                          index=fill_series.index)
    #Drop values where both sets are nan, no interpolation possible here
    df=df.dropna(how='all')
    #Group blocks where contiguous values are Nan in the series to fill
       
    df['nanmask']=df['fill_series'].isnull().astype(bool)
    df['nanmask_shift']=df['nanmask'].shift(1) 
    df['nan_change']=df['nanmask'] != df['nanmask_shift']
    df['nangroup']=df['nan_change'].cumsum()

    holes=df[['fill_series','other_series','nangroup','nanmask']].groupby(['nangroup','nanmask'])
    #calculate regression coefficients for holes from other_series
    for count,hole in holes:
        if hole.nanmask.all():
            slope,intercept,rvalue,pvalue,stderr=stats.linregress(range(len(hole.fill_series)),
                                                                  hole.other_series)
            
            df.ix[hole.index,'fill_series']= intercept + slope * hole.other_series                                                     
        
    return df['fill_series']
    '''        
    testdata=dataframe.testdata.values
    refdata=dataframe.refdata.values
    refdata_ones = np.vstack([refdata, np.ones(len(refdata))]).T
 
    ress=np.linalg.lstsq(refdata_ones,testdata)[0][::-1]  
    dataframe['ones']=1
    xm=np.matrix(dataframe.as_matrix(columns=['ones','refdata']))
    out=np.dot(xm,np.matrix(ress).transpose())
    #df[fill_series]= intercept + slope * xi
    '''
    #fill holes in max_series from regression coefficients found

def merge_ts(dataframe):
    '''
    Merge temporally coinciding timeseries via interpolation from
    bivariate linear regression.
    Input: Dataframe of temporally coinciding timeseries.
    -Find TS with most values. Use a reference to enhance.
    -Enhance missing values in reference via linear regression from other TS.
        Only use correlating (R>0.8) TS for interpolation.
    -Return gap filled reference TS
    '''
    max_series_name=np.argmax(dataframe.notnull().sum().to_frame()[0])
    quick_corr=dataframe.corr().sort(max_series_name,ascending=False).index.values[1:]
    #Order other TSs via correlation to max_series and return list corr_ordered
    #Improve max_series first with best fitting other,than second best fitting other,...
    DF_other_series=pd.DataFrame()
    for other_series_name in dataframe[quick_corr]:
        corr,pval=stats.pearsonr(dataframe[[max_series_name,other_series_name]].dropna()[max_series_name],
                                 dataframe[[max_series_name,other_series_name]].dropna()[other_series_name])
        if corr >=0.8 and pval < 0.01:
            DF_other_series[other_series_name]=dataframe[other_series_name]
            
            #dataframe=dataframe.rename(columns={max_series_name:'refdata',other_series_name:'testdata'})
            
            #out,R,pval,ress=regress(dataframe)
            #DF_Sensor[max_series_name+'_new']=out
         
            
    max_series=dataframe[max_series_name]
    other_series_merged=DF_other_series.mean(axis=1)
    if max_series.isnull().all():
        return None
    elif other_series_merged.isnull().all():
        #If there is no second TS, return the unchanged reference series
        return max_series
    else:
        #Else perform interpolation from linear regression
        merged_ts=fill_holes_lress(max_series,other_series_merged)
        return merged_ts

            