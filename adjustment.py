# -*- coding: utf-8 -*-
"""
Created on Aug 18 12:27 2017

@author: wpreimes
"""
from datetime import datetime
import numpy as np
import scipy.stats as stats
import pandas as pd
import os
import matplotlib.pyplot as plt


class LinearAdjustment(object):
    def __init__(self, data, breaktime,
                 adjust_part='first', adjust_param='both', model_plots=False, maxplotiter=4):

        self.data = data.copy()
        self.breaktime = breaktime
        self.adjust_part = adjust_part
        self.adjust_param = adjust_param
        self.maxplotliter = maxplotiter

        self.model_plots = self._load_model_plot(model_plots)

        self.residuals, self.B = self.regress_model() # create B
        self.calc_correction_values()

    def _load_model_plot(self, model_plots):
        if model_plots:
            if model_plots==True:
                fig, axs = plt.subplots(self.maxplotliter+1, 2,
                                        figsize=(10, self.maxplotliter*4),
                                        facecolor='w', edgecolor='k')
                self.plotcounter = 0
                return (fig, axs)
            elif isinstance(model_plots, tuple):
                self.plotcounter = model_plots[0]
                return model_plots[1]
            else: return False
        else:
            return False

    def get_lin_model_params(self):
        lin_model_params = {'SlopeRatio': self.slope_ratio,
                            'InterceptDiff': self.intercept_diff,
                            'Slope_before_break': self.B[0, 0],
                            'Intercept_before_break': self.B[1, 0],
                            'Slope_after_break': self.B[0, 1],
                            'Intercept_after_break': self.B[1, 1],
                            'NormalizedResidualsSquareSum_before_break': self.norm_residuals_square_sums[0],
                            'NormalizedResidualsSquareSum_after_break': self.norm_residuals_square_sums[1],
                            'MaxResidual_before_break': self.max_abs_residuals[0],
                            'MaxResidual_after_break': self.max_abs_residuals[1]}

        return lin_model_params


    def check_residuals(self):
        return
        # TODO: THis works only for dropping elements of the first group
        previous_residuals = adjust_obj.order_sq_e
        # drop_el = int(len(previous_residuals) * 0.1)
        # drop_el = drop_el if drop_el > 0 else 1
        drop_el = retries  # drop one more element for each retrie
        drop_elements = adjust_obj.order_sq_e[-drop_el:]
        index_drop_elements = self.data.index[drop_elements]
        data_resampled = self.data.drop(index_drop_elements)

    def update_lin_model_params(self):
        # cannot be implemented in this class, main?
        '''
        Dont overwrite results but combine results of multiple iterations
        :return:
        '''
        pass

    def save_plots(self, path, name):
        if self.model_plots:
            plt.tight_layout()
            plt.savefig(os.path.join(path, name + '.png'), dpi=50)
            plt.close()
        else:
            return

    def calc_correction_values(self):
        '''
        Calculate Differences in the 2 models (before/after breaktime)
        Save the model parameters
        :return: dit
            return Model parameters and difference parameters
        '''
        if self.adjust_part == 'last':
            # Perform the linear adjustment for testdata after the breaktime
            self.slope_ratio = self.B[0, 0] / self.B[0, 1]  # k1/k2
            self.intercept_diff = self.B[1, 0] - self.slope_ratio * self.B[1, 1] # d1 - cc*d2

        elif self.adjust_part == 'first':
            # Perform the linear adjustment for testdata before the breaktime
            # B =   | k1  k2|
            #       | d1  d2 |
            self.slope_ratio = self.B[0, 1] / self.B[0, 0]  # k2/k1
            self.intercept_diff = self.B[1, 1] - self.slope_ratio * self.B[1, 0] # d2 - cc*d1
        else:
            raise Exception("select 'first' or 'last' for part to adjust")

    def lms_model(self):
        pass

    def regress_model(self, data='M'):
        '''
        Model Data before/after breaktime via linear model, return Model parameters
        :param data: pd.DataFrame
        :param breaktime: datetime.datetime or str
        :return: np.matrix, dict
        '''
        dataframe = self.data.copy()

        dataframe = dataframe.dropna()
        if isinstance(self.breaktime, str):
            breaktime = datetime.strptime(self.breaktime, '%Y-%m-%d')


        i1 = dataframe[:self.breaktime]
        i2 = dataframe[self.breaktime+pd.DateOffset(1):]

        B_dict = {'b1':None, 'b2':None}
        rval = []
        pval = []

        norm_residuals_square_sums=[]
        max_abs_residuals=[]
        s02 = []


        for i,data in enumerate([i1, i2]):
            testdata = data['testdata']
            refdata = data['refdata']

            r, p = stats.pearsonr(refdata.values, testdata.values)
            rval.append(r)
            pval.append(p)

            if any(r < 0 for r in rval):
                raise Exception('2: negative corr %s (%s)' % ('before breaktime' if i==0  else 'after breaktime',
                                                                     str(rval)))
            if any(p > 0.1 for p in pval):  # todo: war 0.05
                raise Exception('3: positive corr %s not significant (%s)' % ('before breaktime' if i==0  else 'after breaktime',
                                                                              str(pval)))

            n = len(testdata)
            X = np.transpose(np.stack(( data['refdata'].values, np.ones(data['refdata'].values.size))))
            N = np.dot(np.transpose(X), X)
            try:
                b = np.dot(np.dot(np.linalg.inv(N), np.transpose(X)), testdata.values)
            except np.linalg.LinAlgError:
                raise Exception('6: N Matrix singular')

            y_model = np.dot(b, np.transpose(X))
            e = np.dot(-1 * X, b) + testdata.values
            #if i==0:
            #    order_sq_e = np.argsort(e**2) # last index is element with largest squared residual of group1

           #TODO: Check the residuals and rerun the modelling if some residuals are too large
           # self.check_residuals(testdata, e, breaktime)

            norm_residuals_square_sums.append(np.sum(np.square(e))/len(e))
            max_abs_residuals.append(max(abs(e)))
            s02.append(np.dot(np.transpose(e),e) / (n-2))

            B_dict['b%i' % i] = b

            dataframe.loc[testdata.index, 'e'] = e
            dataframe.loc[testdata.index, 'y_model'] = y_model

            if self.model_plots:
                if self.plotcounter > self.maxplotliter:
                    self.plotcounter = self.maxplotliter
                    self.model_plots[1][self.plotcounter, 0].clear()
                    self.model_plots[1][self.plotcounter, 1].clear()

                self.model_plots[1][self.plotcounter, i].scatter(refdata, testdata)
                self.model_plots[1][self.plotcounter, i].plot(refdata, y_model, 'b')
                #self.model_axs[self.model_plot_iteration, i].plot((np.min(testdata), np.max(refdata)), "r--")
                self.model_plots[1][self.plotcounter, i].axis('equal')
                if self.plotcounter == 0:
                   self.model_plots[1][self.plotcounter, i].set_title('Model %s' % 'before break' if i == 0 else 'after break')
                self.model_plots[1][self.plotcounter, i].set_xlabel('Reference Data')
                self.model_plots[1][self.plotcounter, i].set_ylabel('Test Data CCI SM')

                textbox ='n: %i \n' \
                         'b1: %f \n' \
                         'b2: %f \n' % (n, b[0], b[1])

                self.model_plots[1][self.plotcounter, i].annotate(textbox,
                                                                  fontsize=10,
                                                                  xy=(0.7, 0.1),
                                                                  xycoords='axes fraction',
                                                                  bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 3})



        self.B = np.matrix([[B_dict['b0'][0], B_dict['b1'][0]],
                            [B_dict['b0'][1], B_dict['b1'][1]]])
        self.corr = {'R': rval, 'p': pval}
        self.norm_residuals_square_sums = norm_residuals_square_sums
        self.max_abs_residuals = max_abs_residuals
        self.s02 = s02
        #self.order_sq_e = order_sq_e

        return dataframe[['e', 'y_model']], self.B


    def adjust_data(self, data_to_adjust, return_part=None):
        '''
        Adjust the object data based on the parameters in B

        :param data_to_adjust:
        :param return_part: string
            'all' to return adjusted dataframe over whole time frame
            'first' to return only data before breaktime
            'last' to return only data after breaktime
        :return: pd.DataFrame, dict
        '''
        if not return_part:
            return_part = self.adjust_part

        dataframe = data_to_adjust.copy()


        part1 = dataframe.loc[:self.breaktime, 'testdata']
        part2 = dataframe.loc[self.breaktime + pd.DateOffset(1):, 'testdata']

        cc = self.slope_ratio
        dd = self.intercept_diff

        if self.adjust_part == 'last':
            if self.adjust_param == 'd': #TODO: This is deprecated
                #meandiff = part1.mean() - part2.mean()
                meandiff = dd
                adjusted_part2 = part2 + meandiff
            else:
                adjusted_part2 = cc * part2 + dd
            dataframe['adjusted'] = pd.concat([part1, adjusted_part2])

        elif self.adjust_part == 'first':
            if self.adjust_param == 'd': #TODO: This is deprecated
                #meandiff = part2.mean() - part1.mean()
                meandiff = dd
                adjusted_part1 = part1 + meandiff
            else:
                adjusted_part1 = cc * part1 + dd
            dataframe['adjusted'] = pd.concat([adjusted_part1, part2])
        else:
            raise Exception("select 'first' or 'last' for part to adjust")


        if return_part == 'both':
            return dataframe['adjusted']
        elif return_part == 'last':
            return dataframe['adjusted'][self.breaktime+pd.DateOffset(1):]
        elif return_part == 'first':
            return dataframe['adjusted'][:self.breaktime]
        else:
            raise Exception("Select 'first' , 'last' or 'both' for part to return")


    def extract_adjust_infos(self, data_dict):
        #TODO: Delete this
        # type: (dict) -> dict
        status = int(data_dict['adj_status'][0])
        if status not in [0]:
            return {'slope': np.nan, 'intercept': np.nan,
                    'k1_before': np.nan, 'k2_before': np.nan,
                    'd1_before': np.nan, 'd2_before': np.nan,
                    'k1_after': np.nan, 'k2_after': np.nan,
                    'd1_after': np.nan, 'd2_after': np.nan,
                    'adj_status': status}
        else:
            data_dict['adj_status'] = status

            return data_dict  # DataDict from Adjustment incl. status

    def get_meta(self):
        adjustment_meta = {0 : 'no adjustment performed (initial)',
                           1 : 'adjusted data for time frame saved',
                           2 : 'negative correlation',
                           3 : 'positive correlation not significant',
                           4 : 'max number of iterations reached',
                           5 : 'max. iter. reached w.o improvements',
                           6 : 'N Matrix singular',
                           7 : ' ',
                           8 : 'adjusted Data is not being stored',
                           9 : ' '}
        return adjustment_meta


if __name__=='__main__':
    model_plots=True

    xx = np.array([-0.51, 51.2])
    yy = np.array([0.33, 51.6])
    means = [xx.mean(), yy.mean()]
    stds = [xx.std() / 3, yy.std() / 3]
    corr = 0.8  # correlation
    covs = [[stds[0] ** 2, stds[0] * stds[1] * corr],
            [stds[0] * stds[1] * corr, stds[1] ** 2]]

    m = np.random.multivariate_normal(means, covs, 367).T
    for retries in [0,1]:
        data = pd.DataFrame(index=pd.DatetimeIndex(start=datetime(2000,1,1), end=datetime(2001,1,1), freq='D'),
                            data = {'refdata': m[0], 'testdata': m[1]} )
        adjust_obj = LinearAdjustment(data,
                                      data.resample('M').mean(),
                                      datetime(2000,6,1),
                                      'first',
                                      'both',
                                      model_plots if retries == 0 else (retries,adjust_obj.model_plot))