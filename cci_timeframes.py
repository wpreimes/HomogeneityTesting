# -*- coding: utf-8 -*-
"""
Created on Jul 28 16:03 2017

@author: wpreimes
"""
import numpy as np
from datetime import datetime

def get_timeframes(product, skip_times=None, as_datetime=False):
    # TODO: Add timeframes for active and passive products for 22 and 33
    # TODO: Check timeframe and breaktimes for CCI33

    #Order matters!!!!
    timeframes = {'cci_22_combined': np.array(
                                      [['2011-10-01', '2014-12-31'], ['2007-01-01', '2012-07-01'],
                                       ['2002-07-01', '2011-10-01'], ['1998-01-01', '2007-01-01'],
                                       ['1991-08-01', '2002-07-01'], ['1987-07-01', '1998-01-01']]),
                  'cci_31_combined': np.array(
                                      [['2011-10-01', '2015-05-01'], ['2010-07-01', '2012-07-01'],
                                       ['2007-10-01', '2011-10-01'], ['2007-01-01', '2010-07-01'],
                                       ['2002-07-01', '2007-10-01'], ['1998-01-01', '2007-01-01'],
                                       ['1991-08-01', '2002-07-01'], ['1987-09-01', '1998-01-01']]),
                  'cci_31_passive': np.array(
                                     [['2012-07-01', '2015-12-31'], ['2011-10-01', '2015-05-01'],
                                      ['2010-07-01', '2012-07-01'], ['2007-10-01', '2011-10-01'],
                                      ['2002-07-01', '2010-07-01'], ['1998-01-01', '2007-10-01'],
                                      ['1987-09-01', '2002-07-01']]),
                  'cci_31_active': np.array(
                                     [['2007-01-01', '2015-12-01'], ['2003-02-01', '2012-11-01'],
                                      ['1997-05-01', '2007-01-01'], ['1991-08-01', '2003-02-01']]),
                  'cci_33_combined': np.array(
                                     [['2011-10-05', '2016-12-31'], ['2010-01-15', '2012-07-01'],
                                      ['2007-10-01', '2011-10-05'], ['2002-06-18', '2010-01-15'],
                                      ['1998-01-01', '2007-10-01'], ['1987-07-09', '2002-06-18'],
                                      ['1978-10-26', '1998-01-01']]),
                  'cci_40_combined': np.array(
                                    [['2011-10-05', '2016-12-31'], ['2010-01-15', '2012-07-01'],
                                     ['2007-10-01', '2011-10-05'], ['2007-01-01', '2010-01-15'],
                                     ['2002-06-19', '2007-10-01'], ['1998-01-01', '2007-01-01'],
                                     ['1991-08-05', '2002-06-19'], ['1987-07-09', '1998-01-01'],
                                     ['1978-10-26', '1991-08-05']])
               }

    breaktimes = {'cci_22_combined': np.array(
                                      ['2012-07-01', '2011-10-01', '2007-01-01', '2002-07-01', '1998-01-01',
                                       '1991-08-01']),
                  'cci_31_combined': np.array(
                                      ['2012-07-01', '2011-10-01', '2010-07-01', '2007-10-01', '2007-01-01',
                                       '2002-07-01', '1998-01-01', '1991-08-01']),
                  'cci_31_passive': np.array(
                                      ['2015-05-01', '2012-07-01', '2011-10-01', '2010-07-01', '2007-10-01',
                                       '2002-07-01', '1998-01-01']),
                  'cci_31_active': np.array(['2012-11-01', '2007-01-01', '2003-02-01', '1997-05-01']),
                  'cci_33_combined': np.array(
                                      ['2012-07-01', '2011-10-05', '2010-01-15', '2007-10-01',
                                       '2002-06-18', '1998-01-01', '1987-07-09']),
                  'cci_40_combined': np.array(
                                      ['2012-07-01', '2011-10-05', '2010-01-15', '2007-10-01', '2007-01-01',
                                       '2002-06-19', '1998-01-01', '1991-08-05', '1987-07-09']
                  )
                }


    ranges = {'cci_22_combined': ['1980-01-01', '2014-12-31'],
              'cci_31_combined': ['1980-01-01', '2015-12-31'],
              'cci_31_passive': ['1980-01-01', '2015-12-31'],
              'cci_31_active': ['1991-08-01', '2015-12-31'],
              'cci_32_combined': ['1980-01-01', '2015-12-31'],
              'cci_33_combined': ['1978-10-26', '2016-12-31'],
              'cci_40_combined': ['1978-10-26', '2016-12-31']}

    # TODO: Delete [+4]!!
    timeframes['adjusted_cci']=timeframes['cci_31_combined']
    ranges['adjusted_cci'] = ranges['cci_31_combined']
    breaktimes['adjusted_cci'] = breaktimes['cci_31_combined']
    #---------------------


    if product not in timeframes.keys() or product not in breaktimes.keys() or product not in ranges.keys():
        raise Exception('No test times for selected cci product')
    return_timeframes = timeframes[product]
    return_breaktimes = breaktimes[product]
    return_ranges = ranges[product]

    if return_timeframes.shape[0] != return_breaktimes.size:
        raise Exception('Number of timeframes %i and breaktimes %i not equal' % (return_timeframes.size,
                                                                                 return_breaktimes.size))

    if skip_times:
        return_timeframes = np.array([return_timeframes[i] for i in range(return_timeframes.shape[0]) if i not in skip_times])
        return_breaktimes = np.array([return_breaktimes[i] for i in range(return_breaktimes.size) if i not in skip_times])


    if as_datetime:
        return_breaktimes = [datetime.strptime(breaktime, '%Y-%m-%d') for breaktime in return_breaktimes]
        return_timeframes = [[datetime.strptime(time, '%Y-%m-%d') for time in timeframe] for timeframe in return_timeframes]
        return_ranges = [datetime.strptime(range, '%Y-%m-%d') for range in return_ranges]

    return {'breaktimes': return_breaktimes, 'timeframes': return_timeframes, 'ranges': return_ranges}
