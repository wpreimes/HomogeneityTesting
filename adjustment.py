# -*- coding: utf-8 -*-
"""
Created on Thu Jun 08 10:16:42 2017

@author: wpreimes
"""
import pandas as pd
import numpy as np
from scipy import stats


class Adjustment(object):
    def __init__(self,dataframe,testresults):
        if type(testresults['wilkoxon']['h']) == int:
            self.h_wk = testresults['wilkoxon']['h']
        if type(testresults['fligner_killeen']['h']) == int:
            self.h_fk = testresults['fligner_killeen']['h']
        self.data = dataframe


def adjust_ts(self):
    #Perform adjustment of Timeseries AFTER breaktime (if inhomogeneity exists data AFTER the breaktime
    #is matched to data BEFORE breaktime)  
    if (self.h_fk == 1 or self.h_wk == 1):

        i1 = self.data.loc[self.timeframe[0]:self.breaktime]
        i2 = self.data.loc[self.breaktime+pd.DateOffset(1):self.timeframe[1]]

        B = []
        rval = []
        pval = []
        for data in [i1,i2]:
            refdata = np.stack((np.ones(data.refdata.values.size),data.refdata.values))
            testdata = data.testdata.values
            b = np.dot(np.linalg.inv(np.asmatrix(np.dot(refdata,np.transpose(refdata)))),np.matrix.transpose(np.asmatrix(np.dot(refdata,testdata))))
            B.append(b)
            r,p = stats.pearsonr(refdata[1],testdata)
            rval.append(r)
            pval.append(p)
    
        B = np.squeeze(np.asarray(B))
    
        regtest = {'R':rval, 'p':pval}
        print('Regression coefficients are: '+str(regtest))
        
        if any(r < 0 for r in rval):
            print('negative correlation found, adjustment NOT performed')
            raise Exception('negative correlation found, adjustment NOT performed')
            # TODO:WICHTIG: Hier war urspünglich > 0.05:??????
        if any(p > 0.05 for p in pval):
            print 'positive correlation not significant, adjustment NOT performed'
            raise Exception, 'positive correlation not significant, adjustment NOT performed'
        
        #Perform the linear adjustment for testdata after the breaktime
        cc = B[0][1] / B[1][1]
        dd = B[0][0] - cc*B[1][0]
        
        adjusted_part1 = self.data[['testdata']][self.timeframe[0]:self.breaktime]
        adjusted_part2 = cc * self.data[['testdata']][self.breaktime:self.timeframe[1]:] + dd
        self.data['adjust'] = pd.concat([adjusted_part1, adjusted_part2])
        self.data['adjusted'] = self.data.adjust[self.breaktime:]
        adj_setting = {'slope':cc, 'intercept':dd, 'model':B}
        '''
        print 'Adjustment was performed using the following parameters:'
        print 'Slope: %f' %adj_setting['slope']
        print 'Intercept: %f' %adj_setting['intercept']
        print 'Model: {0} and {1}'.format(adj_setting['model'][0],adj_setting['model'][1])
        '''
        return adj_setting, self.data['adjusted']


    else:
        return False,None

