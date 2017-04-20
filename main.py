# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 17:05:55 2017

@author: wpreimes
"""

'''
Breaktimes sind Punkte an denen Inhomogenitäten vermutet werden
Laut Processing overview und Research letter:
    Combined Data: (ignore September 1987) , August 1991, January 1998, July 2002, Januar 2007, October 2011, July 2012, (additional May 2015)

used alpha in research letter:0.01
'''

import sys

workpath=r"H:\workspace"
if workpath not in sys.path:
    sys.path.append(workpath)
    
from HomogeneityTesting.interface import homog_test
from HomogeneityTesting.otherplots import inhomo_plot_with_stats,calc_longest_homogeneous_period,show_tested_gpis
import numpy as np
import os
from matplotlib.colors import LinearSegmentedColormap
from datetime import datetime


def start(test_prod,ref_prod,QDEG_gpi_csv,workpath):
    
    if test_prod in ['cci_22_from_file','cci_22']:
        
        breaktimes=['2002-07-01','2011-10-01','1991-08-01',
                    '2012-07-01','2007-01-01','1998-01-01']
                    
        timeframes=[['1998-01-01','2007-01-01'],
                    ['2007-01-01','2012-07-01'],
                    ['1987-07-01','1998-01-01'],
                    ['2011-10-01','2015-01-01'],
                    ['2002-07-01','2011-10-01'],
                    ['1991-08-01','2002-07-01']]
                    
    elif test_prod in ['cci_31_from_file','cci_31']:
        
        breaktimes=['2002-07-01','2007-01-01','2007-10-01',
                    '2010-07-01','2011-10-01','2012-07-01',
                    '2015-05-01','1991-08-01','1998-01-01']
                    
        timeframes=[['1998-01-01','2007-01-01'],
                    ['2002-07-01','2007-10-01'],
                    ['2007-01-01','2010-07-01'],
                    ['2007-10-01','2011-10-01'],
                    ['2010-07-01','2012-07-01'],
                    ['2011-10-01','2015-05-01'],
                    ['2012-07-01','2015-12-31'],
                    ['1987-09-01','1998-01-01'],
                    ['1991-08-01','2002-07-01']]
    else:
        raise Exception, 'Test product unknown'
        
    i=1
    while os.path.exists(os.path.join(workpath,'v'+str(i))):
        i+=1
    else:
        os.makedirs(os.path.join(workpath,'v'+str(i)))
        
    workpath=os.path.join(workpath,'v'+str(i))
        
    for breaktime,timeframe in zip(breaktimes,timeframes):
            try:
                test_obj=homog_test(QDEG_gpi_csv,
                                    test_prod=test_prod,
                                    ref_prod=ref_prod,
                                    timeframe=timeframe,                    
                                    breaktime=breaktime,
                                    alpha=0.01,
                                    workpath=workpath)
            except: continue
            
            
            test_obj.DF_Points['h_all']=np.nan
            test_obj.DF_Points['message']='Not processed'
            
            print 'Start testing'
                            
            for iteration,gpi in enumerate(test_obj.DF_Points.index.values):
                if iteration%1000 == 0: 
                    print 'Processing QDEG Point %i (iteration %i of %i)' %(gpi,iteration,test_obj.DF_Points.index.values.size)
                
                if test_obj.ref_prod == 'ISMN-merge':               
                    valid_insitu_gpis=test_obj.ismndata.gpis_with_netsta
                    
                    if gpi not in valid_insitu_gpis.keys():
                        continue
                
                try:
                    #test_obj.save_as_mat(gpi=gpi)   
                    testresult=test_obj.run_tests(gpi=gpi,
                                                 FK_Test=True,
                                                 WK_Test=True)
                                                 
                    test_obj.DF_Points.set_value(gpi,'h_FK',
                                                 testresult['FlignerKilleen']['h'])
                                                 
                    test_obj.DF_Points.set_value(gpi,'h_WK',
                                                 testresult['Wilkoxon']['h'])
                    
                    test_obj.DF_Points.set_value(gpi,'message',
                                                 'Processing OK')               
                except Exception as e:
                    #In case something went wrong
                    test_obj.DF_Points.set_value(gpi,'h_FK',np.nan)
                    test_obj.DF_Points.set_value(gpi,'h_WK',np.nan)
                    test_obj.DF_Points.set_value(gpi,'message',str(e))
            
                                
                wk=test_obj.DF_Points.h_WK.loc[gpi]
                fk=test_obj.DF_Points.h_FK.loc[gpi]
                if wk == 1 and fk==0:
                    test_obj.DF_Points.set_value(gpi,'h_all',1.0)
                elif wk == 0 and fk==1:
                    test_obj.DF_Points.set_value(gpi,'h_all',2.0)
                elif wk == 1 and fk==1:
                    test_obj.DF_Points.set_value(gpi,'h_all',3.0)
                elif wk == 0 and fk==0:
                    test_obj.DF_Points.set_value(gpi,'h_all',4.0)

            #Add Info to log file
            test_obj.add_log_line('Finished testing for timeframe %s at %s' %(timeframe,datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            
            
            test_obj.DF_Points=test_obj.DF_Points[['cell','lat','lon','h_WK','h_FK','h_all','message']]
            test_obj.DF_Points.to_csv(os.path.join(test_obj.workpath,'DF_Points_%s_%s_%s_%s.csv' %(ref_prod,timeframe[0],breaktime,timeframe[1])),
                                      index=True, na_rep='nan') 
                                      
            test_obj.add_log_line('Saved results to: DF_Points_%s_%s_%s_%s.csv' %(ref_prod,timeframe[0],breaktime,timeframe[1]))
            

    #Plotting
    show_tested_gpis(test_obj.workpath,test_obj.ref_prod)
    inhomo_plot_with_stats(test_obj.workpath)
    
    test_obj.add_log_line('Created plots for Homogeneity Testing results and Tested GPIs')
    
    
    test_obj.add_log_line('=====================================')
    calc_longest_homogeneous_period(test_obj.workpath,test_obj.test_prod,test_obj.ref_prod)   
    test_obj.add_log_line('Created Plot for Longest Homogeneous Period')  


    #TODO: Create Scatterplot for comparing RTG and RTM



#Refproduct must be one of gldas-merged,gldas-merged-from-file,merra2,ISMN-merge
'''
start('cci_31','ISMN-merge',
      r"H:\workspace\HomogeneityTesting\csv\pointlist_global_quarter.csv",
      r'H:\workspace\HomogeneityTesting\output')  
      
start('cci_22','ISMN-merge',
      r"H:\workspace\HomogeneityTesting\csv\pointlist_global_quarter.csv",
      r'H:\workspace\HomogeneityTesting\output')  

'''
start('cci_31','merra2',
      r"H:\workspace\HomogeneityTesting\csv\pointlist_global_quarter.csv",
      r'H:\workspace\HomogeneityTesting\output')                          
'''
start('cci_22','merra2',
      r"H:\workspace\HomogeneityTesting\csv\pointlist_global_quarter.csv",
      r'H:\workspace\HomogeneityTesting\output')    

start('cci_31','gldas-merged-from-file',
      r"H:\workspace\HomogeneityTesting\csv\pointlist_global_quarter.csv",
      r'H:\workspace\HomogeneityTesting\output')    
      
start('cci_22','gldas-merged-from-file',
      r"H:\workspace\HomogeneityTesting\csv\pointlist_global_quarter.csv",
      r'H:\workspace\HomogeneityTesting\output')    
'''