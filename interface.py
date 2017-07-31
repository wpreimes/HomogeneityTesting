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

from cci_timeframes import get_timeframes
from otherfunctions import regress, datetime2matlabdn
from import_satellite_data import QDEGdata_M, QDEGdata_D
from import_ismn_data import ISMNdataUSA

warnings.simplefilter(action="ignore", category=RuntimeWarning)


class HomogTest(object):
    def __init__(self, test_prod, ref_prod, alpha, anomaly):
        # type: (str,str,float,Union(bool,str)) -> None
        '''
        :param test_prod:
        :param ref_prod:
        :param timeframe:
        :param breaktime:
        :param alpha:
        :param anomaly: To use anomaly data choose 'timeframe' or 'ccirange', else False
        '''
        self.ref_prod = ref_prod
        self.test_prod = test_prod

        test_prod_times = get_timeframes(test_prod, as_datetime=True)
        self.breaktimes =test_prod_times['breaktimes']
        self.timeframes = test_prod_times['timeframes']
        self.range = test_prod_times['ranges']

        '''
        if self.ref_prod == 'ISMN-merge':
            self.data = QDEGdata_M(products=[test_prod])
            self.ismndata = ISMNdataUSA(self.timeframe, self.breaktime, max_depth=0.1)
        else:
        '''
        self.data = QDEGdata_M(products=[self.ref_prod, self.test_prod])

        self.alpha = alpha
        self.anomaly = anomaly

        self._init_validate_cci_range()

    def _init_validate_cci_range(self):
        # TODO: Not needed anymore
        # type: (None) -> None
        '''
        cci_re = {name: re.compile("%s.+" % name) for name in valid_ranges.keys()}
        if not any([cci_re[version].match(self.test_prod) for version, x in cci_re.iteritems()]):
            raise Exception('Unknown Test Product')
        else:
        
        prefix, version, type = self.test_prod.split('_')
        name = prefix + '_' + version
        '''
        for timeframe in self.timeframes:
            for startend in timeframe:
                if not self.range[0] <= startend <= self.range[1]:
                    raise Exception('Selected Timeframe is not valid for product %s' % self.test_prod)
        for breaktime in self.breaktimes:
            if not self.range[0] <= breaktime <= self.range[1]:
                raise Exception('Selected Breaktimes not valid for product %s' % self.test_prod)


    @staticmethod
    def wk_test(dataframe, alternative='two-sided'):
        # type: (pd.DataFrame, str) -> (float,dict)

        T_wk, p_wk = stats.mannwhitneyu(dataframe['Q'].loc[dataframe['group'] == 0],
                                        dataframe['Q'].loc[dataframe['group'] == 1],
                                        alternative=alternative)[1]
        stats_wk = stats.ranksums(dataframe['Q'].loc[dataframe['group'] == 0],
                                  dataframe['Q'].loc[dataframe['group'] == 1])[0]
        

        return p_wk, T_wk, stats_wk

    @staticmethod
    def fk_test(dataframe, mode='median', alpha=0):
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

    def run_tests(self, data, breaktime, min_data_size, tests):
        # type: (pd.DataFrame, datetime, int, np.array) -> (pd.DataFrame, dict)

        # Check if any data is left for testdata and reference data
        if data.isnull().all().refdata or data.isnull().all().testdata:
            raise Exception('1: No coinciding data for the selected timeframe')

        # Check if lengths of remaining datasets equal
        if data['refdata'].index.size != data['testdata'].index.size:
            raise Exception('2: Test timeseries and reference timeseries do not match')

        # Calculation of Spearman-Correlation coefficient
        corr, pval = stats.spearmanr(data['testdata'], data['refdata'], nan_policy='omit')  #TODO: omit or propagate?

        # Check the rank correlation so that correlation is positive and significant
        if not (corr > 0 and pval < 0.01): #TODO: stricter thresholds?
            raise Exception('3: Spearman correlation failed with correlation %f \ '
                            '(must be >0)and pval %f (must be <0.01)' % (corr, pval))

        # Relative Test
        # Divide Time Series into subgroubs according to timeframe and breaktime
        data['group'] = np.nan

        i1 = data.loc[:breaktime]
        i2 = data.loc[breaktime + pd.DateOffset(1):]

        data.loc[i1.index, 'group'] = 0
        data.loc[i2.index, 'group'] = 1
        ni1 = len(i1.index)
        ni2 = len(i2.index)

        df_time = data.dropna() #TODO: again?

        # Check if group data sizes are above selected minimum size
        if ni1 < min_data_size or ni2 < min_data_size:
            raise Exception('4: Minimum Dataseries Length not reached. Size is %i and/or %i !> %i'
                            % (ni1, ni2, min_data_size))

        df_time['bias_corr_refdata'], rxy, pval, ress = regress(df_time)

        if any(np.isnan(ress)):
            raise Exception('5: Negative or NaN correlation after refdata correction')

        # Calculate difference TimeSeries
        df_time['Q'] = df_time.testdata - df_time.bias_corr_refdata

        # Wilcoxon rank sum test
        if 'wk' in tests:
            try:
                p_wk, T_wk, stats_wk = self.wk_test(df_time, 'two-sided')

                if p_wk < self.alpha:
                    h_wk = 1
                else:
                    h_wk = 0

                wilkoxon = {'h': h_wk, 'zval': stats_wk, 'T': T_wk, 'p': p_wk}
            except:
                wilkoxon = {'h': 'Error during WK testing'}
                pass
        else:
            wilkoxon = {'h': 'WK not selected'}

        if 'fk' in tests:
            try:
                h_fk, stats_fk = self.fk_test(df_time[['Q', 'group']], mode='median', alpha=self.alpha)
                fligner_killeen = {'h': h_fk, 'stats': stats_fk}
            except:
                fligner_killeen = {'h': 'Error during FK testing'}
                pass
        else:
            fligner_killeen = {'h': 'FK not selected'}

        if type(wilkoxon['h']) is str and type(fligner_killeen['h'] is str):
            raise Exception('6: WK test and FK test failed')

        return df_time, {'wilkoxon': wilkoxon, 'fligner_killeen': fligner_killeen}
