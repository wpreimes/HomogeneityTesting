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
from collections import Counter
from datetime import datetime

from scipy.stats import fligner, levene
from cci_timeframes import CCITimes
from otherfunctions import regress, datetime2matlabdn
from import_satellite_data import QDEGdata_D
from import_ismn_data import ISMNdataUSA

warnings.simplefilter(action="ignore", category=RuntimeWarning)


class BreakTestData(object):
    '''
    Class containing properties of data used for break detection
    '''
    def __init__(self, test_prod, ref_prod, anomaly):
        # type: (str,str,Union(bool,str)) -> None
        '''
        :param test_prod: str
            Name of the product as in QDEGdata to test for breaks against the reference
            e.g for CCI: CCI_*version*_*PRODUCT*
        :param ref_prod: str
            Name of the reference product as in QDEGdata
            e.g merra2
        :param anomaly: bool
            Set True to use anomaly data for testing
        '''
        self.ref_prod = ref_prod
        self.test_prod = test_prod
        self.range = CCITimes(self.test_prod, ignore_position=True).get_times(None, as_datetime=False)['ranges']
        self.anomaly = anomaly
        if self.ref_prod == 'ISMN-Merge':
            self.ismndata = ISMNdataUSA('merra2', max_depth=0.1)
            self.data = QDEGdata_D(products=[self.test_prod])
        else:
            self.data = QDEGdata_D(products=[self.ref_prod, self.test_prod])


    def read_gpi(self, gpi, start=None, end=None):
        '''
        Read time series for the class products from start date to end date

        :param gpi: int
            index of ground point for which data is read
        :param start: string (%Y-%m-%d)
            Start date of the time series to read
        :param end: string (%Y-%m-%d)
            End date of the time series to read
        :return: pd.DataFrame
            DataFrame of SM values of the test product and ref product
        '''

        # Import the test data and reference datasets for the active ground point
        if not start:
            start = self.range[0]
        if not end:
            end = self.range[1]

        if self.anomaly == 'ccirange':
            # Calculate the anomaly over the whole CCI version time frame (1978-present)
            range = [time.strftime('%Y-%m-%d') for time in self.range]
            try:
                df_time = (self.data).read_gpi(gpi, range[0], range[1])
                df_time = df_time / 100  # type: pd.DataFrame

                df_time[self.ref_prod] = calc_anomaly(df_time[self.ref_prod])
                df_time[self.test_prod] = calc_anomaly(df_time[self.test_prod])
            except:
                raise Exception('9: Could not import data for gpi %i' % gpi)

            if self.ref_prod == 'ISMN-Merge':
                print('CCI range anomaly wont work with ISMN data')

        else:
            try:
                if self.ref_prod == 'ISMN-Merge':
                    df_time = self.data.read_gpi(gpi, start, end)

                    df_time = df_time / 100  # type: pd.DataFrame

                    df_time['ISMN-Merge'] = self.ismndata.read_gpi(gpi, start, end)
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

    def group_by_breaktime(self, df_time_in, breaktime, min_data_size, ignore_exception=False):
        '''
        Divide Time Series into 2 subgroups according to breaktime (before/after)

        :param df_time: pd.DataFrame
        :param breaktime: datetime
        :return: pd.DataFrame
        '''
        df_time = df_time_in.copy()
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

    @staticmethod
    def temp_resample(df_in, how='M', threshold=None):
        '''
        Resample a dataframe to monthly values, if the number of valid values (not nans) in a month
        is smaller than the defined threshold, the monthly resample will be NaN

        :param how: str
            Time frame for temporal resampling, M = monthly, 10D = 10daily,...
        :param threshold: float
            % of valid days (not nan) in timeframe defined in 'how'
        :return: pd.DataFrame
            The monthly resampled Data
        '''
        df = df_in.copy()

        #TODO: Move this outside the function
        # Check if any data is left for testdata and reference data
        if df.isnull().all().refdata or df.isnull().all().testdata:
            raise Exception('1: No data for the selected timeframe')

        if not threshold:
            return df.resample(how).mean()
        else:
            if how!='M':
                raise NotImplementedError

            years, months = df.index.year, df.index.month
            startday = datetime(years[0], months[0], 1)
            last_year, last_month = years[-1], months[-1]

            if last_month == 12:
                next_month, next_year = 1 , last_year + 1
            else:
                next_month, next_year = last_month + 1, last_year

            days_last_month =  (datetime(next_year, next_month, 1) - datetime(last_year, last_month, 1)).days
            endday = datetime(last_year, last_month, days_last_month)

            index_full = pd.DatetimeIndex(start=startday, end = endday, freq = 'D')
            df_alldays = pd.DataFrame(index=index_full,
                                      data = {'dCount' : range(index_full.size)})

            df = pd.concat([df, df_alldays],axis=1)
            concat = []
            for column in df.columns.values:
                if column == 'dCount':continue
                resampled = df[[column, 'dCount']].groupby(pd.TimeGrouper(how)) \
                    .filter(lambda g: g.count()[column] > round(g.count()['dCount'] * threshold)) \
                    .resample(how).mean()
                concat.append(resampled[[column]])

            return pd.concat(concat, axis=1)

    def ref_data_correction(self, df_time):
        '''
        Scale the column "refdata" in the given dataframe to values of "testdata" in the dataframe

        :param df_time: pd.DataFrame
            SM values, columns must be named 'refdata' and 'testdata'
        :return: pd.DataFrame
            DataFrame with bias corrected reference data
        '''
        data = df_time.copy()
        adjusted_data, rxy, pval, ress = regress(data)

        data['bias_corr_refdata'] = adjusted_data

        if any(np.isnan(ress)):
            raise Exception('5: Negative or NaN correlation after refdata correction')

        return data['bias_corr_refdata']

    @staticmethod
    def filter_by_quantiles(df_in, lower=0.1, upper=.9):
        '''
        Mask data outside of the definded quantile range
        '''
        df = df_in.copy()
        upper_threshold = df.quantile(upper)
        lower_threshold = df.quantile(lower)
        df.loc[:, 'diff_flag'] = np.nan
        index_masked=df.query('Q < %f & Q > % f' % (upper_threshold, lower_threshold)).index
        #index_masked = df[(df[colname] < upper_threshold) & (df[colname] > lower_threshold)].index
        df.loc[index_masked, 'diff_flag'] = 0
        return df['diff_flag']

    def quantile_filtering(self, data_in, breaktime, parts, lower_quantile, upper_quantile):
        '''
        quantile filtering based of differences between refdata and testdata
        '''
        data = data_in.copy()
        if parts == 'both':
            filter_mask = pd.concat([self.filter_by_quantiles(data.loc[:breaktime,'Q'].to_frame(),
                                                             lower_quantile,
                                                             upper_quantile),
                                    self.filter_by_quantiles(data.loc[breaktime + pd.DateOffset(1):,'Q'].to_frame(),
                                                             lower_quantile,
                                                             upper_quantile)],
                                    axis=0)
        elif parts == 'first':
            filter_mask = pd.concat([self.filter_by_quantiles(data.loc[:breaktime,'Q'].to_frame(),
                                                             lower_quantile,
                                                             upper_quantile),
                                    pd.Series(index = data.loc[breaktime + pd.DateOffset(1):].index, data=0).to_frame()],
                                    axis=0)
        elif parts == 'last':
            filter_mask = pd.concat([pd.Series(index = data.loc[:breaktime].index, data=0).to_frame(),
                                     self.filter_by_quantiles(data.loc[breaktime + pd.DateOffset(1):,'Q'].to_frame(),
                                                             lower_quantile,
                                                             upper_quantile)],
                                    axis=0)

        data['diff_flag'] = filter_mask
        return data.loc[data['diff_flag'] == 0]


class BreakTestBase(BreakTestData):
    '''
    Class containing functions and properties for relative homogeneity testing
    '''
    def __init__(self, test_prod, ref_prod, tests, alpha, anomaly):
        # type: (str,str,dict,float,Union(bool,str)) -> None
        '''
        :param test_prod: str
            as in BreakTestData
        :param ref_prod: str
            as in BreakTestData
        :param tests: dict
            dictionary of test types and test names as implemented in self.run_tests()
            eg. {'mean':wilkoxon, 'var':'scipy_fligner_killeen'}
        :param alpha: float
            significance level for all tests
        :param anomaly: str
            as in BreakTestData
        '''
        BreakTestData.__init__(self, test_prod, ref_prod, anomaly)
        self.tests = tests
        self.fligner_approx = 'chi' #TODO: defined by user, implement fisher
        self.alpha = alpha

    def get_testresults(self):
        return self.testresults

    def check_testresult(self, testresult):
        '''
        Checks if all tests are negative (no break)

        :param testresult: dict
            test results as returned by run_tests
        :return: bool
        '''
        test_h = {test: testresult['h_%s' % test] for test in self.tests.values()}
        if all([h == 0 for h in test_h.values()]):
            return None, None, False
        else:
            break_found_by = []
            failed_tests = []
            for test in self.tests.values():
                if test_h[test] == 1:
                    break_found_by.append(test)
                elif test_h[test] == 99:
                    failed_tests.append(test)

            has_break = True if break_found_by else False

            return failed_tests if failed_tests else None, \
                   break_found_by if break_found_by else None, \
                   has_break

    def compare_testresults(self, reference_result, last_result, priority='mean'):
        ref_tests_failed, ref_found_by, ref_has_break = self.check_testresult(reference_result)
        last_tests_failed, last_found_by, last_has_break = self.check_testresult(last_result)

        if ref_has_break and not last_has_break:
            # Break was removed --> better
            return True
        elif ref_has_break and last_has_break:
            # Both found break, use priority to decide
            # compare tests
            if Counter(ref_found_by) == Counter(last_found_by):  # Same tests found break
                return False  # nothing improved --> worse
            elif priority:
                for p in priority:
                    if (self.tests[p] in ref_found_by) and (self.tests[p] not in last_found_by):
                        return True  # --> Priority was removed, other might have been added
                    else:
                        continue  # Move to next Priority
                return False  # --> Priority was not removed
            else:
                return False # --> no priority defined


    @staticmethod
    def wk_test(dataframe, alternative='two-sided', alpha=0.01):
        # type: (pd.DataFrame, str) -> (float,dict)

        p_wk = stats.mannwhitneyu(dataframe['Q'].loc[dataframe['group'] == 0],
                                  dataframe['Q'].loc[dataframe['group'] == 1],
                                  alternative=alternative)[1]
        stats_wk = stats.ranksums(dataframe['Q'].loc[dataframe['group'] == 0],
                                  dataframe['Q'].loc[dataframe['group'] == 1])[0]  # type: dict

        if p_wk < alpha:
            h = 1
        else:
            h = 0

        return h, {'zval': stats_wk, 'pval': p_wk}

    @staticmethod
    def fk_test(dataframe_in, mode='median', alpha=0.01):
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
        dataframe = dataframe_in.copy()
        df = dataframe.rename(columns={'Q': 'data'})
        df = df.dropna()
        N = df.index.size
        K = df['group'].nunique()

        df['A'] = np.nan

        if mode == 'median':
            for i in range(K):
                subset = df.ix[df['group'] == i]
                groupmed = subset.data.median()
                df.ix[df['group'] == i, 'groupmed'] = groupmed  # group medians
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
    def scipy_fk_test(dataframe, mode='median', alpha=0.1):
        df = dataframe.rename(columns={'Q': 'data'})
        df = df.dropna()
        sample1 = df[df['group'] == 0.]['data'].values
        sample2 = df[df['group'] == 1.]['data'].values

        stats, pval = fligner(sample1, sample2, center=mode)

        stats_fk = {'chi': {'z': stats, 'pval': pval}}

        if stats_fk['chi']['pval'] < alpha:
            h = 1
        else:
            h = 0
        return h, stats_fk

    @staticmethod
    def lv_test(dataframe, mode='median', alpha=0.1):
        #TODO: Not tested
        df = dataframe.rename(columns={'Q': 'data'})
        df = df.dropna()
        sample1 = df[df['group' == 0.].index]['Q'].values
        sample2 = df[df['group' == 1.].index]['Q'].values

        stats, pval = levene(sample1, sample2, center=mode)

        stats_lv = {'chi': {'z': stats, 'pval': pval}}

        if stats_lv['chi']['pval'] < alpha:
            h = 1
        else:
            h = 0
        return h, stats_lv

    def check_corr(self, df_time):
        # Calculation of Spearman-Correlation coefficient
        corr, pval = stats.spearmanr(df_time['testdata'], df_time['refdata'], nan_policy='omit')

        # Check the rank correlation so that correlation is positive and significant
        if not (corr > 0 and pval < 0.01):  # TODO: stricter thresholds?
            raise Exception('3: Spearman correlation failed with correlation %f ' 
                            '(must be >0) and pval %f (must be <0.01)' % (corr, pval))

        return corr, pval

    def restructure_test_results(self, test_results):
        # type: (dict) -> dict

        restructured_results = {}

        for test in self.tests.values():
            restructured_results['h_%s' % test] = test_results[test]['h']
            if self.fligner_approx in test_results[test]['stats'].keys():
                stats = test_results[test]['stats'][self.fligner_approx]
                restructured_results['z_%s' % test] = stats['z']
            else:
                stats = test_results[test]['stats']
            restructured_results['p_%s' % test] = stats['pval']

        test_status = 0 # '0: Testing successful'

        if 'mean' in self.tests.keys():
            mean_test_result = restructured_results['h_%s' % self.tests['mean']]
            if np.isnan(mean_test_result):
                test_status = 6 # '6: WK was selected but failed.'
                restructured_results['h_%s' % self.tests['mean']] = 99
        else:
            mean_test_result = np.nan
        if 'var' in self.tests.keys():
            var_test_result = restructured_results['h_%s' % self.tests['var']]
            if np.isnan(var_test_result):
                test_status = 7 # '7: FK was selected but failed.'
                restructured_results['h_%s' % self.tests['var']] = 99
        else:
            var_test_result = np.nan

        # Combine results to single value
        #TODO: Do this when plotting data
        if mean_test_result == 1:
            if var_test_result == 1:
                all = 3.0
            else:
                all = 1.0
        elif var_test_result == 1:
            all = 2.0
        else:
            all = 4.0

        restructured_results.update({'test_results': all, 'test_status': test_status})

        return restructured_results

    def run_tests(self, data):
        # type: (pd.DataFrame, list) -> (pd.DataFrame, dict)
        '''
        Prepares Data for Testing. Bias Correction of Reference Data. Analyzes data sufficiency.
        :param data:
        :param breaktime:
        :param min_data_size:
        :return:
        '''
        tests = self.tests.values()
        testresults = {}

        # Wilcoxon rank sum test
        if 'wilkoxon' in tests:
            try:
                h_wk, stats_wk = self.wk_test(data, 'two-sided')
                wilkoxon = {'h': h_wk, 'stats': stats_wk}
            except:
                wilkoxon = {'h': np.nan, 'stats': np.nan}
                pass
            testresults['wilkoxon'] = wilkoxon
        if 'fligner_killeen' in tests:
            try:
                h_fk, stats_fk = self.fk_test(data[['Q', 'group']], mode='median', alpha=self.alpha)
                fligner_killeen = {'h': h_fk, 'stats': stats_fk}
            except:
                fligner_killeen = {'h': np.nan, 'stats': np.nan}
                pass
            testresults['fligner_killeen'] = fligner_killeen
        if 'scipy_fligner_killeen' in tests:
            try:
                h_fk, stats_fk = self.scipy_fk_test(data[['Q', 'group']], mode='median', alpha=self.alpha)
                fligner_killeen = {'h': h_fk, 'stats': stats_fk}
            except:
                fligner_killeen = {'h': np.nan, 'stats': np.nan}
                pass
            testresults['scipy_fligner_killeen'] = fligner_killeen

        self.testresults = self.restructure_test_results(testresults)

        return self.testresults

    def get_meta(self):
        test_meta = {1 : 'WK only',
                     2 : 'FK only',
                     3 : 'WK and FK',
                     4 : 'None'}
        status_meta = {0 : 'Testing successful',
                       1 : 'No data for the selected timeframe',
                       2 : ' ',
                       3 : 'Spearman correlation too low',
                       4 : 'Min. Dataseries len. not reached',
                       5 : 'neg/nan correl. aft. bias corr.',
                       6 : 'WK was selected but failed',
                       7 : 'FK was selected but failed',
                       8 : 'WK test and FK test failed',
                       9 : 'Could not import data for gpi'}
        return test_meta, status_meta

