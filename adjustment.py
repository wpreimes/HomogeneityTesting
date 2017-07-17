# -*- coding: utf-8 -*-
"""
Created on Thu Jun 08 10:16:42 2017

@author: wpreimes
"""
import pandas as pd
import numpy as np
from scipy import stats
import os

from HomogeneityTesting.interface import HomogTest
from HomogeneityTesting.otherfunctions import cci_timeframes
from HomogeneityTesting.points_to_netcdf import pd_from_2Dnetcdf

class Adjust(object):
    def __init__(self, workdir, test_prod, ref_prod, anomaly):
        self.workdir = workdir
        self.test_prod = test_prod
        self.ref_prod = ref_prod
        self.anomaly = anomaly

        mergetimes = cci_timeframes(self.test_prod)
        self.breaktimes = mergetimes['breaktimes'][::-1]
        self.timeframes = mergetimes['timeframes'][::-1]

        self.test_results = self._init_import_test_results()


    def _init_import_test_results(self):
        dataframes = []
        for breaktime in self.breaktimes:
            filename = 'HomogeneityTest_%s_%s.nc' % (self.ref_prod, breaktime)
            testresults = pd_from_2Dnetcdf(os.path.join(self.workdir, filename),
                                           grid='land')[['test_results']]
            testresults= testresults.rename(columns={'test_results': breaktime})
            dataframes.append(testresults)

        return pd.concat(dataframes, axis=1)


    def load_data(self, timeframe, breaktime):
        test_obj = HomogTest(self.test_prod, self.ref_prod, timeframe,
                             breaktime, )
        pass

def create_empty_ncfile(self):
    pass


def adjust_ts(self, data, breaktime):
    # Perform adjustment of Timeseries AFTER breaktime (if inhomogeneity exists data AFTER the breaktime
    # is matched to data BEFORE breaktime)

    # data = dataframe with columns "refdata" and "testdata"

    if (self.h_fk == 1 or self.h_wk == 1):
        # Set 1= Data <=breaktime
        # Set 2= Data after breaktime

        i1 = data.loc[:self.breaktime]
        i2 = data.loc[self.breaktime + pd.DateOffset(1):]

        B = []
        rval = []
        pval = []
        for data in [i1, i2]:
            refdata = np.stack((np.ones(data['refdata'].values.size), data['refdata'].values))
            testdata = data['testdata'].values
            b = np.dot(np.linalg.inv(np.asmatrix(np.dot(refdata, np.transpose(refdata)))),
                       np.matrix.transpose(np.asmatrix(np.dot(refdata, testdata))))
            B.append(b)
            r, p = stats.pearsonr(refdata[1], testdata)
            rval.append(r)
            pval.append(p)

        B = np.squeeze(np.asarray(B))

        regtest = {'R': rval, 'p': pval}
        # print('Regression coefficients are: ' + str(regtest))

        if any(r < 0 for r in rval):
            raise Exception('negative correlation found, adjustment NOT performed')
        if any(p > 0.05 for p in pval):
            raise Exception('positive correlation not significant, adjustment NOT performed')

        # Perform the linear adjustment for testdata before the breaktime
        cc = B[0][1] / B[1][1]
        dd = B[0][0] - cc * B[1][0]

        adjusted_part1 = data[['testdata']][:self.breaktime]
        adjusted_part2 = cc * self.data[['testdata']][self.breaktime:self.timeframe[1]:] + dd
        self.data['adjust'] = pd.concat([adjusted_part1, adjusted_part2])
        self.data['adjusted'] = self.data.adjust[self.breaktime:]
        adj_setting = {'slope': cc, 'intercept': dd, 'model': B}
        '''
        print 'Adjustment was performed using the following parameters:'
        print 'Slope: %f' %adj_setting['slope']
        print 'Intercept: %f' %adj_setting['intercept']
        print 'Model: {0} and {1}'.format(adj_setting['model'][0],adj_setting['model'][1])
        '''
        return adj_setting, self.data['adjusted']


    else:
        return False, None

adjust_obj = Adjust(r"H:\workspace\HomogeneityTesting\output\CCI31EGU",'cci_31_combined','merra2',
                    anomaly=False)
