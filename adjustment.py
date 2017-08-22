# -*- coding: utf-8 -*-
"""
Created on Aug 18 12:27 2017

@author: wpreimes
"""
from datetime import datetime
import numpy as np
import scipy.stats as stats
import pandas as pd

def adjustment(data_daily, B, breaktime, adjust_param = 'both', adjust_part='first', return_part='all',
               check_adjustment=False):
    '''
    :param data: pd.DataFrame
        from Homogeneity Testing, must not contain nan
    :param B: np.matrix
    :param breaktime: string or datetime
    :param adjust_param: 'add' for adjusting only additive bias, 'mult' for adjusting only multiplicative bias
                         'both' for adjusting both model parameters
    :param adjust_part: string
        'first' to adjust values BEFORE break time
        'last' to adjust values AFTER break time
    :param return_part: string
        'all' to return adjusted dataframe over whole time frame
        'first' to return only data before breaktime
        'last' to return only data after breaktime
    :return: pd.DataFrame, dict
    '''

    dataframe = data_daily.copy()


    if adjust_part == 'last':
        # Perform the linear adjustment for testdata after the breaktime
        cc = B[0][1] / B[1][1]
        dd = B[0][0] - cc * B[1][0]
        adjusted_part1 = dataframe[['testdata']][:breaktime]

        if adjust_param == 'both':
            adjusted_part2 = cc * dataframe[['testdata']][breaktime:] + dd
        elif adjust_param == 'add':
            adjusted_part2 = dataframe[['testdata']][breaktime:] + dd
        elif adjust_param == 'multi':
            adjusted_part2 = cc * dataframe[['testdata']][breaktime:]
        adjusted_part2 = adjusted_part2[breaktime + pd.DateOffset(1):]


    elif adjust_part == 'first':
        # Perform the linear adjustment for testdata before the breaktime
        cc = B[1][1] / B[0][1]
        dd = B[1][0] - cc * B[0][0]
        adjusted_part2 = dataframe[['testdata']][breaktime + pd.DateOffset(1):]
        if adjust_param == 'both':
            adjusted_part1 = cc * dataframe[['testdata']][:breaktime] + dd
        elif adjust_param == 'add':
            adjusted_part1 = dataframe[['testdata']][:breaktime] + dd
        elif adjust_param == 'multi':
            adjusted_part1 = cc * dataframe[['testdata']][:breaktime]
        adjusted_part1 = adjusted_part1

    else:
        raise Exception("Select 'first' or 'last' for part to adjust")

    dataframe['adjusted'] = pd.concat([adjusted_part1, adjusted_part2])

    if check_adjustment:
        tolerance = 0.05 # TODO: für monthly weniger streng als für daily
        B, corr = adjustment_params(dataframe[['adjusted','refdata']].rename(columns={'adjusted':'testdata'}),
                                    breaktime)

        print 'B after adjustment: %s' %str(B.flatten())

        p1_B1_after = B[0][0]
        p1_B2_after = B[0][1]

        p2_B1_after = B[1][0]
        p2_B2_after = B[1][1]

        if not np.isclose(p1_B1_after, p2_B1_after, atol=tolerance) or \
                not np.isclose(p1_B2_after, p2_B2_after, atol=tolerance):
            # !!!!!!!!!!RECURSIVE PART!!!!!!!!!!!!!!!
            print 'B tolerance after adjustment not reached, retrying'
            # Take adjusted data from first iteration
            dataframe = dataframe[['adjusted', 'refdata']].rename(columns={'adjusted': 'testdata'})
            # Retry adjustment with adjusted data and new B
            adjusted, adj_setting = adjustment(dataframe, B, breaktime, adjust_param, adjust_part, return_part,
                                               check_adjustment=True)

            dataframe['adjusted'] = adjusted

            if return_part == 'all':
                return dataframe['adjusted'], adj_setting
            elif return_part == 'last':
                return dataframe['adjusted'][breaktime+pd.DateOffset(1):], adj_setting
            elif return_part == 'first':
                return dataframe['adjusted'][:breaktime], adj_setting
            else:
                raise Exception("Select 'first' , 'last' or 'all' for part to return")

    # !!!!!!!!!!ESCAPE PART!!!!!!!!!!!!!!!
    adj_setting = {'slope': cc, 'intercept': dd, 'part1_B1': B[0][0], 'part1_B2': B[0][1],
                   'part2_B1': B[1][0], 'part2_B2': B[1][1]}

    if return_part == 'all':
        return dataframe['adjusted'], adj_setting
    elif return_part == 'last':
        return dataframe['adjusted'][breaktime+pd.DateOffset(1):], adj_setting
    elif return_part == 'first':
        return dataframe['adjusted'][:breaktime], adj_setting
    else:
        raise Exception("Select 'first' , 'last' or 'all' for part to return")

def adjustment_params(data, breaktime):
    '''
    Model Data before/after breaktime via linear model, return Model paramters
    :param data: pd.DataFrame
    :param breaktime: datetime.datetime or str
    :return: np.matrix, dict
    '''
    dataframe = data.copy()
    dataframe = dataframe.dropna()
    if isinstance(breaktime, str):
        breaktime = datetime.strptime(breaktime, '%Y-%m-%d')


    i1 = dataframe[:breaktime]
    i2 = dataframe[breaktime:]

    # TODO: What happens if timeframe before / after breaktime is not equally long

    # if i1.index.size != i2.index.size:
    #     if i1.index.size < i2.index.size:
    #         i2 = i2[:i1.index.size]
    #     elif i1.index.size > i2.index.size:
    #         i1 = i1[-1 * i2.index.size:]

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

    #TODO: activate tests

    if any(r < 0 for r in rval):
        raise Exception('2: negative correlation found (%s)' % str(rval))
    if any(p > 0.1 for p in pval): # todo: war 0.05
        raise Exception('3: positive correlation not significant (%s)' % str(pval))

    return B, regtest



if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')

    from interface import HomogTest
    from otherfunctions import regress
    import matplotlib.dates as dates
    import matplotlib.pyplot as plt
    import pandas as pd
    import scipy.stats as stats
    import numpy as np
    from cci_timeframes import get_timeframes

    test_obj = HomogTest('cci_31_combined',
                         'merra2',
                         ['wilkoxon', 'fligner_killeen'],
                         0.01,
                         False,
                         None)

    df_full = test_obj.read_gpi(433291, start=test_obj.range[0], end=test_obj.range[1])
    for i, (timeframe, breaktime) in enumerate(zip(test_obj.timeframes, test_obj.breaktimes)):
        print ' '
        print ' '
        print '================%s===============' % str(breaktime)
        sett = False
        j=0
        while (not sett) and j < 5:
            df = df_full[['testdata','refdata']][timeframe[0]:timeframe[1]]
            df = df.dropna()
            try:
                corr, pval = test_obj.check_corr(df)
                df, _, _ = test_obj.group_by_breaktime(df, breaktime, min_data_size=3,
                                                         ignore_exception=False)  # TODO: different value?
            except:
                print 'TESTING FAILED'
                break
            # Correct bias in reference data
            df['refdata'] = test_obj.ref_data_correction(df[['testdata', 'refdata']])
            # Calculate difference TimeSeries
            df['Q'] = df['testdata'] - df['refdata']
            corr, p = test_obj.check_corr(df)
            print '------------BEFORE ADJUSTMENT-----------'
            print ('Correlation total: corr: %f , pval: %f' %(corr, p) )
            # Run Tests
            _, testresult = test_obj.run_tests(data=df)
            print testresult
            if testresult['fligner_killeen']['h'] == 0 and testresult['wilkoxon']['h'] == 0:
                print 'Adjustment not necessary'
                break
            x1 = dates.date2num([pd.Timestamp(date).to_datetime() for date in df[:breaktime].index.values])
            x2 = dates.date2num([pd.Timestamp(date).to_datetime() for date in df[breaktime:].index.values])

            slope1, inter1, r1, p1, stderr1 = stats.linregress(x1, df['testdata'][:breaktime])
            part1 = slope1 * x1 + inter1
            print 'Variance part 1: %f' % np.var(part1)
            print 'Mean part 1: %f' % np.mean(part1)
            slope2, inter2, r2, p2, stderr2 = stats.linregress(x2, df['testdata'][breaktime:])
            part2 = slope2 * x2 + inter2
            print 'Variance part 2: %f' % np.var(part2)
            print 'Mean part 2: %f' % np.mean(part2)
            if part1.size+part2.size != df.index.size:
                part2 = part2[1:]
            df['testdata_ress_before'] = np.concatenate((part1, part2)) # breaktime overlaps, drop it for first set
            df[['testdata', 'testdata_ress_before']].plot()

            print('Testdata_preadjust: PART 1: slope: %f, intercept: %f' % (slope1, inter1))
            print('Testdata_preadjust: PART 2: slope: %f, intercept: %f' % (slope2, inter2))
            print '------------AFTER ADJUSTMENT-----------'
            try:
                B, regtest = adjustment_params(df[['testdata', 'refdata']],
                                               breaktime)


                ts_adjusted, sett = adjustment(df[['testdata', 'refdata']], B, breaktime,
                                               adjust_param='both', adjust_part='first',
                                               return_part='all' if i == 0 else 'first' , check_adjustment=True)
                df_full.loc[ts_adjusted.index, 'testdata'] = ts_adjusted
                j+=1

            except Exception as e:
                print 'Adjustment failed: %s' %e
                continue

        if not sett:
            continue

        print ('regtest: %s' % str(regtest))

        testdata_before = df['testdata'].copy()



        df = df_full[['testdata','refdata']][timeframe[0]:timeframe[1]]
        df=df.dropna()
        try:
            corr, pval = test_obj.check_corr(df)
            df, _, _ = test_obj.group_by_breaktime(df, breaktime, min_data_size=3,
                                                     ignore_exception=False)  # TODO: different value?
        except:
            print 'TESTING FAILED'
            pass
        df['refdata'] = test_obj.ref_data_correction(df[['testdata', 'refdata']])
        df['Q'] = df['testdata'] - df['refdata']
        corr, p = test_obj.check_corr(df)
        print ('Correlation total: corr: %f , pval: %f' %(corr, p) )
        _, testresult = test_obj.run_tests(data=df)  # wilkoxon removed, fk not!!
        print testresult

        x1 = dates.date2num([pd.Timestamp(date).to_datetime() for date in df[:breaktime].index.values])
        x2 = dates.date2num([pd.Timestamp(date).to_datetime() for date in df[breaktime:].index.values])
        slope1, inter1, r1, p1, stderr1 = stats.linregress(x1, df['testdata'][:breaktime])
        part1 = slope1 * x1 + inter1
        print 'Variance part 1: %f' % np.var(part1)
        print 'Mean part 1: %f' % np.mean(part1)
        slope2, inter2, r2, p2, stderr2 = stats.linregress(x2, df['testdata'][breaktime:])
        part2 = slope2 * x2 + inter2

        if part1.size+part2.size != df.index.size:
            part2=part2[1:]
        print 'Variance part 2: %f' % np.var(part2)
        print 'Mean part 2: %f' % np.mean(part2)
        df['testdata_ress_after'] = np.concatenate((part1, part2))
        df[['testdata', 'testdata_ress_after']].plot()
        df['testdata_before'] =testdata_before
        print('Testdata_aftadjust: PART 1: slope: %f, intercept: %f' % (slope1, inter1))
        print('Testdata_aftadjust: PART 2: slope: %f, intercept: %f' % (slope2, inter2))
        #df_roll = pd.rolling_mean(df[['testdata','testdata_before']], 5, 3)
