# -*- coding: utf-8 -*-
"""
Created on Wed Dec 07 13:46:17 2016

@author: wpreimes
"""
import sys
sys.path.append("H:\workspace")

import warnings
from pygeogrids.netcdf import load_grid

import numpy as np
import pandas as pd
import csv
from datetime import datetime,timedelta
import scipy.stats as stats
import scipy.io as io


import os,glob

from HomogeneityTesting.FKtest import FKtest
from HomogeneityTesting.otherfunctions import regress

warnings.simplefilter(action = "ignore", category = RuntimeWarning)

from HomogeneityTesting.import_data import QDEGdata_D,QDEGdata_M
import matplotlib.pyplot as plt
def datetime2matlabdn(dt):
   ord = dt.toordinal()
   mdn = dt + timedelta(days = 366)
   frac = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
   return mdn.toordinal() + frac    

#%% VARIABLE WERTE
#Timeframe legt einen Zeitbereich für H-Test fest
#Choose eiter gldas_combined or merra
class homog_test(object):
    
    def __init__(self,QDEG_gpi_csv,test_prod,ref_prod,timeframe,breaktime,alpha,workpath):
        
        self.workpath=workpath
        DF_Points=pd.read_csv(QDEG_gpi_csv,index_col=0)
        if 'grid' in DF_Points.columns.values:
            del DF_Points['grid']
        DF_Points.index.name='gpi_quarter'
        DF_Points=DF_Points.sort_values('cell',axis='index')
        self.DF_Points=DF_Points
        
        self.ref_prod=ref_prod
        self.test_prod=test_prod        

        self.timeframe=[datetime.strptime(timeframe[0], '%Y-%m-%d'),
                        datetime.strptime(timeframe[1], '%Y-%m-%d')]
        self.breaktime=datetime.strptime(breaktime, '%Y-%m-%d')
        
        
        self.alpha=alpha
        
        
        #Observation and reference data
        if ref_prod == 'insitu':
            self.data=QDEGdata_M(products=[test_prod])
            self.ismndata=ISMNdata(timeframe=self.timeframe,breaktime=self.breaktime,max_depth=0.1)
        
        else:
            self.data=QDEGdata_M(products=[ref_prod,test_prod])
        if self.test_prod not in ['cci_22','cc_31']:
            self.test_prod=self.test_prod[0:6]
        self.check_valid_cci_range()
        #self.testdata_mask=self.calc_testdata_mask(filepath=r'H:\workspace\HomogeneityTesting\output\maskfiles_cci_global\%s'%self.test_prod)
    
    
    def check_valid_cci_range(self):
        valid_ranges={'cci_22':['1980-01-01','2015-12-31'],
                      'cci_31':['1980-01-01','2015-12-31']}
        
        for time in self.timeframe:
            if not valid_ranges[self.test_prod][0]<=time.strftime('%Y-%m-%d')\
                    <=valid_ranges[self.test_prod][1]:
                raise Exception, 'Selected Timeframe is not valid for product %s' %self.test_prod
            else:
                pass
        
   
    def save_as_mat(self,gpi):
        exp_data=self.data.read_gpi(gpi,
                                    self.timeframe[0].strftime('%Y-%m-%d'),
                                    self.timeframe[1].strftime('%Y-%m-%d'))
        exp_data=exp_data/100
                                    

        matdate=[]
        for dt in exp_data[self.test_prod].index:
            matdate.append(datetime2matlabdn(dt))

            
        X={'tspan':matdate,
           'sm':exp_data[self.ref_prod].values}

        Y={'tspan':matdate,
           'sm':exp_data[self.test_prod].values}     
        
        timeframe=[datetime2matlabdn(self.timeframe[0]),
                   datetime2matlabdn(self.timeframe[1])]
        #timeframe=[matdate[0],matdate[-1]]
        breaktime=datetime2matlabdn(self.breaktime)   
        
        data_dict={'X':X,'Y':Y,'timeframe':timeframe,'breaktime':breaktime}       
        io.savemat('matlab_data\SMdata_'+str(gpi),data_dict,oned_as='column')

    '''
    def calc_testdata_mask(self,filepath):

        ttime=[self.timeframe[0].strftime('%Y-%m-%d'),
               self.breaktime.strftime('%Y-%m-%d'),
               self.timeframe[1].strftime('%Y-%m-%d')]
               
        filename=str(self.test_prod)+'_'+ttime[0]+'-'+ttime[1]+'-'+ttime[2]+'_invalid_gpis.csv'
        path=os.path.join(filepath,filename)
        #If a mask file is available, load it
        #If no mask is available, calculate it
        if os.path.isfile(path):
            with open(path,'rb') as maskfile:
                reader=csv.reader(maskfile,delimiter=',')
                masklist=[list(map(int,gpi)) for gpi in list(reader)]
                return masklist[0]
        else:
            #TODO: Are there other land grids?
            
            data=QDEGdata_D(products=[self.test_prod])
            

            all_gpis=self.DF_Points.index.values
                
            bad_gpis=[]
            for i,gpi in enumerate(all_gpis):
                if i%100 == 0: 
                    print '%i of %i' %(i,all_gpis.size)
                try:
                    DF_Time=data.read_gpi(gpi,ttime[0],ttime[2])
                    DF_Time.reindex(pd.date_range(DF_Time.index[0],DF_Time.index[-1]),
                                                  fill_value=np.nan)
                    starttime=DF_Time.index[0]
                    endtime=DF_Time.index[-1]
                    DF_Should= pd.DataFrame(index=pd.date_range(start=datetime(starttime.year,starttime.month,1),
                                                                end=datetime(endtime.year,endtime.month,30),
                                                                freq='D')) 
                    DF_Should['fill']=1                     
                    DF_Time['fill']=1    
                    should_days=DF_Should.fill.groupby([DF_Should.index.year,DF_Should.index.month]).count()*0.3
                    is_days=DF_Time[self.test_prod].groupby([DF_Time.index.year,DF_Time.index.month]).count()
                    if True in (is_days<should_days).values:
                        bad_gpis.append(gpi)
                except: 
                    bad_gpis.append(gpi)

                    
            with open(path,'wb') as maskfile:
                wr=csv.writer(maskfile, quoting=csv.QUOTE_NONE)
                wr.writerow(bad_gpis)
                
            return self.calc_testdata_mask(filepath)
                   
    '''   
                
    def run_tests(self,gpi,FK_Test,WK_Test):    

        min_data_size=3 #Minimale Länge der zu testenden Datenserie
        #%%
        if 'h_WK' not in self.DF_Points.columns:
            self.DF_Points['h_WK']=str(np.nan)
        if 'h_FK' not in self.DF_Points.columns:
            self.DF_Points['h_FK']=str(np.nan) 
        
        ttime=[self.timeframe[0],self.breaktime,self.timeframe[1]]


        #Import the test data and reference datasets for the active ground point
        try:        
            if self.ref_prod == 'insitu':
                DF_Time=self.data.read_gpi(gpi,
                                           ttime[0].strftime('%Y-%m-%d'),
                                           ttime[2].strftime('%Y-%m-%d'))
                                           
                DF_Time['ismn_merge']=self.ismndata.merge_station_series(gpi,
                                           DF_Time[self.test_prod])
                                    
            else:    
                DF_Time=self.data.read_gpi(gpi,
                                           ttime[0].strftime('%Y-%m-%d'),
                                           ttime[2].strftime('%Y-%m-%d'))
                                           
            DF_Time=DF_Time/100
        except:
            raise Exception, 'Could not import data for gpi %i' %gpi
    
                                                    
                #Change NaN placeholder to numpy nans
                #if 'cci' in DF_Time.columns:
                    #DF_Time.loc[DF_Time['cci']== -999999, 'cci'] = np.nan
        
            
        DF_Time = DF_Time.rename(columns={self.ref_prod: 'refdata',
                                          self.test_prod: 'testdata'}) 
            
        #TODO: Gapfill verwenden? Vorteile? Nachteile?    
        #TODO: Anomaly verwenden?
            #for p in products:
                #sm_anomaly=calc_anomaly(DF_Time[p],window_size=31)  
                #DF_Time[p+'_anomaly'] = sm_anomaly
           
         #%%
        #Drop days where either dataset is missing     
        DF_Time=DF_Time.dropna()  
        
        #Check if any data is left for testdata and reference data
        if DF_Time.isnull().all().refdata or DF_Time.isnull().all().testdata:
            raise Exception, 'No coinciding data for the selected timeframe'
            
        #Check if lengths of remaining datasets equal
        if DF_Time['refdata'].index.size != DF_Time['testdata'].index.size:
            raise Exception, 'Test timeseries and reference timeseries do not match'
                
        #TODO: Check if Dataseries coincide in time
           
        #Add consecutive index column
        #DF_Time.loc[:,'ind'] = pd.Series(range(1, DF_Time.index.size+1,1), index=DF_Time.index)
        
        
        #Calculation of Spearman-Correlation coefficient
        corr,pval=stats.spearmanr(DF_Time['testdata'],DF_Time['refdata'])    
        
        # Check the rank correlation so that correlation is positive and significant at 0.05
        #TODO: Keep pval threshold of 0.05? Pub: 0.01 , MATLAB: 0.05
        if not (corr > 0 and pval < 0.01):
            raise Exception, 'Spearman correlation failed with correlation %f (must be >0)and pval %f (must be <0.01)' %(corr,pval)   

        
        
        
        #%% Relative Test
            #Divide Time Series into subgroubs according to timeframe and breaktime
        DF_Time['group']=np.nan
     
        iall=DF_Time.loc[self.timeframe[0]:self.timeframe[1]]
        
        i1=DF_Time.loc[self.timeframe[0]:self.breaktime]
        i2=DF_Time.loc[self.breaktime+pd.DateOffset(1):self.timeframe[1]]            

        
        DF_Time.loc[i1.index,'group']=0
        DF_Time.loc[i2.index,'group']=1
        ni1=len(i1.index)
        ni2=len(i2.index)
        
        DF_Time=DF_Time.dropna()  
        
    
        
        #Check if group data sizes are above selected minimum size
        if ni1<min_data_size or ni2<min_data_size:
           raise Exception, 'Minimum Dataseries Length not reached. Size is %i and/or %i !> %i' %(ni1,ni2,min_data_size)    

        #TODO: Eventuell scipy.stats.linregress direkt implementieren (kein großer Unterschied)
        #Calculate regress as done in myregress.m
        DF_Time['bias_corr_refdata'],rxy,pval,ress=regress(DF_Time)  
        
        if any(np.isnan(ress)):
            raise Exception, 'Negative or NaN correlation'
        
        #Calculate difference TimeSeries
        DF_Time['Q']=DF_Time.testdata-DF_Time.bias_corr_refdata
        
        #%%
        #Wilcoxon rank sum test
        #TODO: stats.wilcoxon würde gleich lange Datenserien benötigen, trotzdem anschauen
        #Alternativ use u, p_value = mannwhitneyu(group1, group2, two-sided test)
        #stats_WK,p_WK=stats.ranksums(DF_Time['Q'].loc[i1.index],DF_Time['Q'].loc[i2.index])
        #print stats_WK
        #print p_WK
        
        #SAME AS:
        
        if WK_Test:
            #Eigentlich funktionieren beide Methoden. p_val unterscheidet sich
            #in 5. Nachkommastelle, stats wird nicht weiter verwendet
            p_WK=stats.mannwhitneyu(DF_Time['Q'].loc[i1.index],DF_Time['Q'].loc[i2.index],
                                    alternative='two-sided')[1]
            stats_WK=stats.ranksums(DF_Time['Q'].loc[i1.index],DF_Time['Q'].loc[i2.index])[0]
            
            if p_WK<self.alpha:
                h_WK=1
            else:
                h_WK=0
            
            #print 'WK Test performed. Result:%i, p:%f, zval:%f' %(h_WK,p_WK,stats_WK)
 
            Wilkoxon={'h':h_WK,'zval':stats_WK,'p':p_WK}
        else:
            Wilkoxon={'h':'WK not selected'}
        

        
       
        if FK_Test:
            h_FK,stats_FK=FKtest(DF_Time[['Q','group']],
                                 mode='median',alpha=self.alpha)
            '''
            print 'FK Test performed. Result:'
            print 'stats: chi: z=%f, df=%i, pval=%f' %(stats_FK['chi']['z'],stats_FK['chi']['df'],stats_FK['chi']['pval'])
            print 'stats: f: chi: z=%f, df=%s, pval=%f' %(stats_FK['f']['z'],str(stats_FK['f']['df']),stats_FK['f']['pval'])
            '''
            #TODO: Fligner-Killeen test
            #stats.fligner('center'=median)  
            FlignerKilleen={'h':h_FK,'stats':stats_FK}
        else:
            FlignerKilleen={'h':'FK not selected'}
        
        return {'Wilkoxon':Wilkoxon, 'FlignerKilleen':FlignerKilleen}

        '''
    #def adjust_TS(self,h_WK,h_FK):
        #%%
        #Perform adjustment of Timeseries AFTER breaktime (if inhomogeneity exists data AFTER the breaktime
        #is matched to data BEFORE breaktime)  
        if (h_FK == 1 or h_WK == 1):
            if h_WK==1:
               print 'WK test found inhomogeneity' 
            if h_FK==1:
               print 'FK test found inhomogeneity'

               
            B=[]
            rval=[]
            pval=[]
            for data in [i1,i2]:
                refdata=np.stack((np.ones(data.refdata.values.size),data.refdata.values))
                testdata=data.testdata.values
                b=np.dot(np.linalg.inv(np.asmatrix(np.dot(refdata,np.transpose(refdata)))),np.matrix.transpose(np.asmatrix(np.dot(refdata,testdata))))
                B.append(b)
                r,p=stats.pearsonr(refdata[1],testdata)
                rval.append(r)
                pval.append(p)
        
            B=np.squeeze(np.asarray(B))
        
            regtest={'R':rval, 'p':pval}   
            print 'Regression coefficients are: '+str(regtest)
            
            if any(r < 0 for r in rval):
                print 'negative correlation found, adjustment NOT performed'
                raise Exception, 'negative correlation found, adjustment NOT performed'
    #TODO:WICHTIG: Hier war urspünglich > 0.05:??????
            if any(p > 0.05 for p in pval):
                print 'positive correlation not significant, adjustment NOT performed'
                raise Exception, 'positive correlation not significant, adjustment NOT performed'
            
            #Perform the linear adjustment for testdata after the breaktime
            cc=B[0][1]/B[1][1]
            dd=B[0][0] - cc*B[1][0]
            
            adjusted_part1=DF_Time[['testdata']][self.timeframe[0]:self.breaktime]
            adjusted_part2=cc*DF_Time[['testdata']][self.breaktime:self.timeframe[1]:]+dd
            DF_Time['adjust']=pd.concat([adjusted_part1, adjusted_part2])
            DF_Time['adjusted']=DF_Time.adjust[self.breaktime:]
            adj_setting={'slope':cc,'intercept':dd,'model':B}
            print 'Adjustment was performed using the following parameters:'
            print 'Slope: %f' %adj_setting['slope']
            print 'Intercept: %f' %adj_setting['intercept']
            print 'Model: {0} and {1}'.format(adj_setting['model'][0],adj_setting['model'][1])
            DF_Time[['bias_corr_refdata','testdata','adjusted']].plot()

            plt.show()
    
        else:
            print 'No inhomogeneity was found between %s and %s' %(str(self.timeframe[0]),
                                                                   str(self.timeframe[1]))
       ''' 
                       

        
    

