# -*- coding: utf-8 -*-
"""
Created on Wed Dec 07 13:46:17 2016

@author: wpreimes
"""
import warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
from typing import Union
from pytesmo.time_series.anomaly import calc_anomaly

from datetime import datetime

from cci_timeframes import CCITimes
from otherfunctions import regress, datetime2matlabdn
from import_satellite_data import QDEGdata_M, QDEGdata_D


warnings.simplefilter(action="ignore", category=RuntimeWarning)


class HomogTest(object):
    def __init__(self, test_prod, ref_prod, tests, alpha, anomaly, adjusted_ts_path=None):
        # type: (str,str,list,float,Union(bool,str)) -> None
        '''
        :param test_prod:
        :param ref_prod:
        :param timeframe:
        :param breaktime:
        :param tests:
        :param alpha:
        :param anomaly: To use anomaly data choose 'timeframe' or 'ccirange', else False
        '''
        self.ref_prod = ref_prod
        self.test_prod = test_prod
        self.range = CCITimes(test_prod, ignore_position=True).get_times(None, as_datetime=False)['ranges']
        self.tests = tests
        '''
        if self.ref_prod == 'ISMN-merge':
            self.data = QDEGdata_M(products=[test_prod])
            self.ismndata = ISMNdataUSA(self.timeframe, self.breaktime, max_depth=0.1)
        else:
        '''
        self.data = QDEGdata_D(products=[self.ref_prod, self.test_prod])

        self.alpha = alpha
        self.anomaly = anomaly

        if adjusted_ts_path:
            self.adjusted_ts_path = adjusted_ts_path



    def check_testresult(self, testresult):
        '''
        Checks if all tests are negative (no break)
        :param testresult: dict
        :return: bool
        '''
        if all(h == 0 for h in [testresult[test]['h'] for test in self.tests]):
            return None, True
        else:
            break_found_by=[]
            for test in self.tests:
                if testresult[test]['h'] == 1:
                    break_found_by.append(test)

            return break_found_by, False


    @staticmethod
    def wk_test(dataframe, alternative='two-sided', alpha=0.01):
        # type: (pd.DataFrame, str) -> (float,dict)

        p_wk = stats.mannwhitneyu(dataframe['Q'].loc[dataframe['group'] == 0],
                                        dataframe['Q'].loc[dataframe['group'] == 1],
                                        alternative=alternative)[1]
        stats_wk = stats.ranksums(dataframe['Q'].loc[dataframe['group'] == 0],
                                  dataframe['Q'].loc[dataframe['group'] == 1])[0] #type: dict

        if p_wk < alpha:
            h = 1
        else:
            h = 0

        return h, {'zval':stats_wk, 'pval': p_wk}

    @staticmethod
    def fk_test(dataframe, mode='median', alpha=0.01):
        # type: (pd.DataFrame, str, float) -> (int,dict)
        '''
        FKTEST Fligner-Killeen test for homogeneity of variances.

         Trujillo-Ortiz, A., R. Hernandez-Walls and N. Castro-Castro. (2009).FKtest:
           Fligner-Killeen test for homogeneity of variances. A MATLAB file. [WWW document].
           URL http://www.mathworks.com/matlabcentral/fileexchange/25040

        Input data format:
            Dataframe: 2 columns
                column1 (data): difference data Q    column2 (group): group number (1 or 2 für reference data or testdata)
        '''

        # number of measurements and datagroups
        df = dataframe.rename(columns={dataframe.ix[:, 0].name: 'data'})
        df = df.dropna()
        N = df.index.size
        K = df['group'].nunique()

        df['A'] = np.nan

        if mode == 'median':
            for i in range(K):
                subset = df.ix[df['group'] == i]
                groupmed = subset.data.median()
                df.ix[df['group'] == i, 'groupmed'] = groupmed  # groupmedians
                df.ix[df['group'] == i, 'groupme_diff'] = np.abs(
                    subset['data'] - groupmed)  # difference data-groupmedians

        if mode == 'mean':
            for i in range(K):
                subset = df.ix[df['group'] == i]
                groupmean = subset.data.mean()
                df.ix[df['group'] == i, 'groupmean'] = groupmean  # groupmeans
                df.ix[df['group'] == i, 'groupme_diff'] = np.abs(
                    subset['data'] - groupmean)  # difference data-groupmeans

        Z = stats.rankdata(df['groupme_diff'])  # score ranking ALL
        sta_norm_dist = stats.norm.ppf(0.5 + (Z / (2. * (N + 1.))))  # score standard normal distribution ALL
        df['A'] = sta_norm_dist
        M = df['A'].mean()  # overall mean

        nn = []
        mm = []
        bb = []
        for i in range(K):
            subset = df.ix[df['group'] == i]

            nn.append(subset.index.size)
            mm.append(np.mean(subset['A']))
            df.ix[df['group'] == i, 'groupAmean'] = mm[i]
            bb.append((nn[i] * (mm[i] - M) ** 2))
            df.ix[df['group'] == i, 'groupB'] = bb[i]

        B = np.array(df['groupB'].unique())
        V = df['A'].var()  # Overall Variance Score
        X2 = np.sum(B) / V  # Fligner-Killeen statistic by the Chi-squared approximation
        v = K - 1  # statistic degree of freedom
        F = (X2 / v) / ((N - 1. - X2) / (N - K))  # Fligner-Killeen statistic by the Fisher approximation

        P1 = 1 - stats.chi2.cdf(X2, v)
        P2 = 1 - stats.f.cdf(F, v, N - K)

        # TODO: Laut Chun-Hsu statt F X2??
        stats_fk = {'chi': {'z': X2, 'df': v, 'pval': P1}, 'f': {'z': F, 'df': [v, N - K], 'pval': P2}}

        if stats_fk['chi']['pval'] < alpha:
            h = 1
        else:
            h = 0

        return h, stats_fk

    @staticmethod
    def vn_test():
        #TODO: implement von Neumann test
        pass

    def read_gpi(self, gpi, start, end):
        # type: (int) -> pd.DataFrame


        # Import the test data and reference datasets for the active ground point
        if self.anomaly == 'ccirange':
            #Calculate the anomaly over the whole CCI version time frame (1978-present)
            range = [time.strftime('%Y-%m-%d') for time in self.range]
            try:
                df_time = (self.data).read_gpi(gpi, range[0], range[1])
                df_time = df_time / 100  # type: pd.DataFrame
                df_time[self.ref_prod] = calc_anomaly(df_time[self.ref_prod])
                df_time[self.test_prod] = calc_anomaly(df_time[self.test_prod])
            except:
                raise Exception('9: Could not import data for gpi %i' % gpi)

            if self.ref_prod == 'ISMN-merge':
                print('CCI range anomaly wont work with ISMN data')

        else:
            try:
                if self.ref_prod == 'ISMN-merge':
                    df_time = self.data.read_gpi(gpi, start, end)

                    df_time = df_time / 100  # type: pd.DataFrame

                    df_time['ISMN-merge'] = self.ismndata.merge_stations_around_gpi(gpi, df_time[self.test_prod])
                else:
                    df_time = self.data.read_gpi(gpi, start, end)
                    df_time = df_time / 100

                if self.anomaly == 'timeframe':
                    df_time[self.ref_prod] = calc_anomaly(df_time[self.ref_prod])
                    df_time[self.test_prod] = calc_anomaly(df_time[self.test_prod])
            except:
                raise Exception('9: Could not import data for gpi %i' % gpi)

        df_time = df_time.rename(columns={self.ref_prod: 'refdata',
                                          self.test_prod: 'testdata'})

        # Drop days where either dataset is missing
        return df_time


    def check_corr(self, df_time):
        # Check if any data is left for testdata and reference data
        if df_time.isnull().all().refdata or df_time.isnull().all().testdata:
            raise Exception('1: No coinciding data for the selected timeframe')

        # Check if lengths of remaining datasets equal TODO: it always is, can be removed
        if df_time['refdata'].index.size != df_time['testdata'].index.size:
            raise Exception('2: Test timeseries and reference timeseries do not match')

        # Calculation of Spearman-Correlation coefficient
        corr, pval = stats.spearmanr(df_time['testdata'], df_time['refdata'], nan_policy='omit')

        # Check the rank correlation so that correlation is positive and significant
        if not (corr > 0 and pval < 0.01): #TODO: stricter thresholds?
            raise Exception('3: Spearman correlation failed with correlation %f \ '
                            '(must be >0)and pval %f (must be <0.01)' % (corr, pval))

        return corr, pval

    def group_by_breaktime(self, df_time, breaktime, min_data_size, ignore_exception=False):
        '''
        Divide Time Series into 2 subgroups according to breaktime (before/after)
        :param df_time: pd.DataFrame
        :param breaktime: datetime
        :return: pd.DataFrame
        '''
        df_time['group'] = np.nan

        i1 = df_time.loc[:breaktime]
        i2 = df_time.loc[breaktime + pd.DateOffset(1):]

        df_time.loc[i1.index, 'group'] = 0
        df_time.loc[i2.index, 'group'] = 1
        ni1 = len(i1.index)
        ni2 = len(i2.index)

        # Check if group data sizes are above selected minimum size
        if not ignore_exception:
            if ni1 < min_data_size or ni2 < min_data_size:
                raise Exception('4: Minimum Dataseries Length not reached. Size is %i and/or %i !> %i'
                                % (ni1, ni2, min_data_size))

        return df_time, ni1, ni2

    def ref_data_correction(self, df_time):
        df_time['bias_corr_refdata'], rxy, pval, ress = regress(df_time)

        if any(np.isnan(ress)):
            raise Exception('5: Negative or NaN correlation after refdata correction')

        return df_time['bias_corr_refdata']

    def run_tests(self, data):
        # type: (pd.DataFrame, list) -> (pd.DataFrame, dict)
        '''
        Prepares Data for Testing. Bias Correction of Reference Data. Analyzes data sufficiency.
        :param data:
        :param breaktime:
        :param min_data_size:
        :return:
        '''
        tests = self.tests
        # Wilcoxon rank sum test
        if 'wilkoxon' in tests:
            try:
                h_wk, stats_wk = self.wk_test(data, 'two-sided')
                wilkoxon = {'h': h_wk, 'stats': stats_wk}
            except:
                wilkoxon = {'h': 'Error during WK testing'}
                pass
        else:
            wilkoxon = {'h': 'WK not selected'}

        if 'fligner_killeen' in tests:
            try:
                h_fk, stats_fk = self.fk_test(data[['Q', 'group']], mode='median', alpha=self.alpha)
                fligner_killeen = {'h': h_fk, 'stats': stats_fk}
            except:
                fligner_killeen = {'h': 'Error during FK testing'}
                pass
        else:
            fligner_killeen = {'h': 'FK not selected'}

        if type(wilkoxon['h']) is str and type(fligner_killeen['h']) is str:
            raise Exception('6: WK test and FK test failed')

        return data, {'wilkoxon': wilkoxon, 'fligner_killeen': fligner_killeen}