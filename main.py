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
        
        breaktimes=['2011-10-01','2002-07-01','1991-08-01',
                    '2012-07-01','2007-01-01','1998-01-01']
                    
        timeframes=[['2007-01-01','2012-07-01'],
                    ['1998-01-01','2007-01-01'],
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

                #if gpi in test_obj.testdata_mask: continue
                if iteration%100 == 0: print 'Processing QDEG Point %i (iteration %i of %i)' %(gpi,iteration,test_obj.DF_Points.index.values.size)
                
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
            
            test_obj.DF_Points=test_obj.DF_Points[['cell','lat','lon','h_WK','h_FK','h_all','message']]
            test_obj.DF_Points.to_csv(os.path.join(test_obj.workpath,'DF_Points_%s_%s_%s_%s.csv' %(ref_prod,timeframe[0],breaktime,timeframe[1])),
                                      index=True, na_rep='nan') 
            
            
            #Plotting
            show_tested_gpis(test_obj.workpath)
            inhomo_plot_with_stats(test_obj.workpath)
            calc_longest_homogeneous_period(test_obj.workpath,test_obj.test_prod,test_obj.ref_prod)
            
            
            
            '''
            colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1),(.5,.5,.5)] 
            cmap=LinearSegmentedColormap.from_list('HomogeneityGroups',
                                                   colors, 
                                                   N=4)
                                                   
            spatial_plot_quarter_grid(test_obj.DF_Points,
                                      tags=['h_all'],
                                      title='HomogeneityTest_%s_(breaktime-%s)'%(ref_prod,test_obj.breaktime.strftime("%Y-%m-%d")),
                                      cbrange=(1,4),
                                      cmap=cmap,
                                      cblabel='1=WK,2=FK,3=both,4=None',
                                      #continent='NA',
                                      path=test_obj.workpath)
            '''

#start('merra2',r"H:\workspace\HomogeneityTesting\csv\pointlist_global_quarter.csv")
#Refproduct: gldas-merged,gldas-merged-from-file---um0h,merra2
'''
start('cci_31','merra2',
      r"H:\workspace\HomogeneityTesting\csv\pointlist_global_quarter.csv",
      r'H:\workspace\HomogeneityTesting\output')                          
'''
start('cci_22','merra2',
      r"H:\workspace\HomogeneityTesting\csv\pointlist_global_quarter.csv",
      r'H:\workspace\HomogeneityTesting\output')    
'''
start('cci_31','gldas-merged-from-file',
      r"H:\workspace\HomogeneityTesting\csv\pointlist_United_quarter.csv",
      r'H:\workspace\HomogeneityTesting\output')    
      
start('cci_22','gldas-merged-from-file',
      r"H:\workspace\HomogeneityTesting\csv\pointlist_Australia_quarter.csv",
      r'H:\workspace\HomogeneityTesting\output')    
'''