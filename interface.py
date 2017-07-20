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

from datetime import datetime

import scipy.io as io

from general.time_series.anomaly import calc_anomaly

from HomogeneityTesting.otherfunctions import regress, datetime2matlabdn, cci_timeframes
from HomogeneityTesting.import_satellite_data import QDEGdata_M, QDEGdata_D
from HomogeneityTesting.import_ismn_data import ISMNdataUSA

warnings.simplefilter(action="ignore", category=RuntimeWarning)


class HomogTest(object):
    def __init__(self, test_prod, ref_prod, timeframe, breaktime, alpha, anomaly):
        # type: (str,str,list,str,float,Union(bool,str)) -> None
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
        if breaktime and timeframe.size == 2:
            self.timeframe = [datetime.strptime(timeframe[0], '%Y-%m-%d'),
                              datetime.strptime(timeframe[1], '%Y-%m-%d')]

            self.breaktime = datetime.strptime(breaktime, '%Y-%m-%d')

            if self.ref_prod == 'ISMN-merge':
                self.data = QDEGdata_M(products=[test_prod])
                self.ismndata = ISMNdataUSA(self.timeframe, self.breaktime, max_depth=0.1)
            else:
                self.data = QDEGdata_M(products=[self.ref_prod, self.test_prod])

        else: # if timeframe and breaktime are None
            self.timeframe = timeframe
            self.breaktime = breaktime
            self.data = QDEGdata_M(products=[self.ref_prod, self.test_prod])

        self.alpha = alpha
        self.anomaly = anomaly
        self.valid_range = self._init_validate_cci_range()

    def _init_validate_cci_range(self):
        # type: (None) -> list
        valid_ranges = cci_timeframes(self.test_prod)['ranges']
        '''
        cci_re = {name: re.compile("%s.+" % name) for name in valid_ranges.keys()}
        if not any([cci_re[version].match(self.test_prod) for version, x in cci_re.iteritems()]):
            raise Exception('Unknown Test Product')
        else:
        
        prefix, version, type = self.test_prod.split('_')
        name = prefix + '_' + version
        '''
        if not self.timeframe: # for adjustment no timeframes are given
            return valid_ranges
        else: # if timeframes are given, do some testing
            for time in self.timeframe:
                if not valid_ranges[0] <= time.strftime('%Y-%m-%d') <= valid_ranges[1]:
                    raise Exception('Selected Timeframe is not valid for product %s' % self.test_prod)
                return valid_ranges

    @staticmethod
    def wk_test(dataframe, alternative='two-sided'):
        # type: (pd.DataFrame, str) -> (float,dict)

        p_wk = stats.mannwhitneyu(dataframe['Q'].loc[dataframe['group'] == 0],
                                  dataframe['Q'].loc[dataframe['group'] == 1],
                                  alternative=alternative)[1]
        stats_wk = stats.ranksums(dataframe['Q'].loc[dataframe['group'] == 0],
                                  dataframe['Q'].loc[dataframe['group'] == 1])[0]

        return p_wk, stats_wk

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

    def read_gpi(self, gpi):
        # type: (int) -> pd.DataFrame
        # Observation and reference data

        if not self.timeframe and not self.breaktime:
            # For adjustment read dara from first to last valid date
            # TODO: Change this to fit with lower ttime
            ttime = [datetime.strptime(self.valid_range[0],'%Y-%m-%d'),
                     None,
                     datetime.strptime(self.valid_range[1], '%Y-%m-%d')]
        else:
            # For break detection data reading
            ttime = [self.timeframe[0], self.breaktime, self.timeframe[1]]

        # Import the test data and reference datasets for the active ground point
        if self.anomaly == 'ccirange':
            if self.ref_prod == 'ISMN-merge':
                print('CCI range anomaly wont work with ISMN data')
                raise Exception('CCI range anomaly wont work with ISMN data')
            else:
                try:
                    df_time = (self.data).read_gpi(gpi, self.valid_range[0], self.valid_range[1])
                    df_time = df_time / 100  # type: pd.DataFrame
                    df_time[self.ref_prod] = calc_anomaly(df_time[self.ref_prod])
                    df_time[self.test_prod] = calc_anomaly(df_time[self.test_prod])
                    df_time = df_time.loc[ttime[0].strftime('%Y-%m-%d'):ttime[2].strftime('%Y-%m-%d')]
                except:
                    raise Exception('9: Could not import data for gpi %i' % gpi)
        else:
            try:
                if self.ref_prod == 'ISMN-merge':
                    df_time = self.data.read_gpi(gpi, ttime[0].strftime('%Y-%m-%d'), ttime[2].strftime('%Y-%m-%d'))

                    df_time = df_time / 100  # type: pd.DataFrame

                    df_time['ISMN-merge'] = self.ismndata.merge_stations_around_gpi(gpi, df_time[self.test_prod])
                else:
                    df_time = self.data.read_gpi(gpi,
                                                 ttime[0].strftime('%Y-%m-%d'),
                                                 ttime[2].strftime('%Y-%m-%d'))
                    df_time = df_time / 100

                if self.anomaly == 'timeframe':
                    # TODO: Gapfill verwenden? Vorteile? Nachteile?
                    df_time[self.ref_prod] = calc_anomaly(df_time[self.ref_prod])
                    df_time[self.test_prod] = calc_anomaly(df_time[self.test_prod])
            except:
                raise Exception('9: Could not import data for gpi %i' % gpi)

        df_time = df_time.rename(columns={self.ref_prod: 'refdata',
                                          self.test_prod: 'testdata'})

        # Drop days where either dataset is missing
        return df_time.dropna()

    def run_tests(self, gpi, tests):
        # type: (int, np.array) -> (pd.DataFrame, dict)

        # Minimale Länge der zu testenden Datenserie
        # TODO: choose larger threshold?
        min_data_size = 3

        df_time = self.read_gpi(gpi)

        # Check if any data is left for testdata and reference data
        if df_time.isnull().all().refdata or df_time.isnull().all().testdata:
            raise Exception('1: No coinciding data for the selected timeframe')

        # Check if lengths of remaining datasets equal
        if df_time['refdata'].index.size != df_time['testdata'].index.size:
            raise Exception('2: Test timeseries and reference timeseries do not match')

        # TODO: Check if Dataseries coincide in time

        # Calculation of Spearman-Correlation coefficient
        corr, pval = stats.spearmanr(df_time['testdata'], df_time['refdata'])

        # Check the rank correlation so that correlation is positive and significant at 0.05
        if not (corr > 0 and pval < 0.01):
            raise Exception('3: Spearman correlation failed with correlation %f \ '
                            '(must be >0)and pval %f (must be <0.01)' % (corr, pval))

        # Relative Test
        # Divide Time Series into subgroubs according to timeframe and breaktime
        df_time['group'] = np.nan

        i1 = df_time.loc[self.timeframe[0]:self.breaktime]
        i2 = df_time.loc[self.breaktime + pd.DateOffset(1):self.timeframe[1]]

        df_time.loc[i1.index, 'group'] = 0
        df_time.loc[i2.index, 'group'] = 1
        ni1 = len(i1.index)
        ni2 = len(i2.index)

        df_time = df_time.dropna()

        # Check if group data sizes are above selected minimum size
        if ni1 < min_data_size or ni2 < min_data_size:
            raise Exception('4: Minimum Dataseries Length not reached. Size is %i and/or %i !> %i'
                            % (ni1, ni2, min_data_size))

        df_time['bias_corr_refdata'], rxy, pval, ress = regress(df_time)

        if any(np.isnan(ress)):
            raise Exception('5: Negative or NaN correlation')

        # Calculate difference TimeSeries
        df_time['Q'] = df_time.testdata - df_time.bias_corr_refdata

        # Wilcoxon rank sum test
        if 'wk' in tests:
            try:
                p_wk, stats_wk = self.wk_test(df_time, 'two-sided')

                if p_wk < self.alpha:
                    h_wk = 1
                else:
                    h_wk = 0

                wilkoxon = {'h': h_wk, 'zval': stats_wk, 'p': p_wk}
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

    def save_as_mat(self, gpi):
        # type: (int) -> None
        # Saves the SM data for the active time frame for the selected gpi as .mat

        print('Exporting Testdata and Referencedata to .mat')
        exp_data = self.data.read_gpi(gpi,
                                      self.timeframe[0].strftime('%Y-%m-%d'),
                                      self.timeframe[1].strftime('%Y-%m-%d'))
        exp_data = exp_data / 100  # type: pd.DataFrame

        matdate = []
        for dt in exp_data[self.test_prod].index:
            matdate.append(datetime2matlabdn(dt))

        x = {'tspan': matdate,
             'sm': exp_data[self.ref_prod].values}

        y = {'tspan': matdate,
             'sm': exp_data[self.test_prod].values}

        timeframe = [datetime2matlabdn(self.timeframe[0]),
                     datetime2matlabdn(self.timeframe[1])]

        breaktime = datetime2matlabdn(self.breaktime)

        data_dict = {'X': x, 'Y': y, 'timeframe': timeframe, 'breaktime': breaktime}
        io.savemat('matlab_data\SMdata_' + str(gpi), data_dict, oned_as='column')
