# -*- coding: utf-8 -*-
"""
Created on Jul 28 16:03 2017

@author: wpreimes
"""
import numpy as np
from datetime import datetime
from smecv_grid.grid import SMECV_Grid_v042

# TODO: Add timeframes for active and passive products for 22 and 33
# TODO: Check timeframe and breaktimes for CCI33

#Order matters!!!!
class CCITimes(object):
    def __init__(self, product):
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
                      'cci_33_active': np.array(
                                         [['1991-07-01','2016-12-31']]),
                      'cci_33_combined': np.array(
                                         [['2011-10-05', '2016-12-31'], ['2010-01-15', '2012-07-01'],
                                          ['2007-10-01', '2011-10-05'], ['2007-01-01', '2010-01-15'],
                                          ['2002-06-19', '2007-10-01'],
                                          {'lat < -37 or lat > 37':
                                               [['1991-08-05', '2007-01-01'],
                                                ['1987-07-09', '2002-06-19']],
                                           'lat >= -37 and lat <= 37':
                                               [['1998-01-01','2007-01-01'],
                                                ['1991-08-05', '2002-06-19'],
                                                ['1987-07-09', '1998-01-01']]
                                           },
                                          ['1978-10-26', '1991-08-05']]),
                      'cci_33_passive': np.array(
                                        [['2011-10-05', '2016-12-31'], ['2010-01-15', '2012-07-01'],
                                         ['2007-10-01', '2011-10-05'], ['2002-06-19', '2010-01-15'],
                                         {'lat < -37 or lat > 37':
                                              [['1987-07-09', '2007-10-01'],
                                               ['1978-10-26', '2002-06-19']],
                                          'lat >= -37 and lat <= 37':
                                              [['1998-01-01', '2007-10-01'],
                                               ['1987-07-09', '2002-06-19'],
                                               ['1978-10-26', '1998-01-01']]
                                          }]),
                      }
        timeframes['cci_42_combined'] = timeframes['cci_33_combined']
        timeframes['cci_42_passive'] = timeframes['cci_33_passive']
        timeframes['cci_42_active'] = timeframes['cci_33_active']

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
                      'cci_33_active': np.array(['2007-01-01']),
                      'cci_33_combined': np.array(
                                          ['2012-07-01', '2011-10-05', '2010-01-15', '2007-10-01',
                                           '2007-01-01', '2002-06-19',
                                           {'lat >= -37 and lat <= 37' : ['1998-01-01']},
                                           '1991-08-05','1987-07-09']),
                      'cci_33_passive': np.array(
                                          ['2012-07-01', '2011-10-05', '2010-01-15', '2007-10-01','2002-06-19',
                                           {'lat >= -37 and lat <= 37' : ['1998-01-01']},
                                           '1987-07-09'])
                    }
        breaktimes['cci_42_combined'] = breaktimes['cci_33_combined']
        breaktimes['cci_42_passive'] = breaktimes['cci_33_passive']
        breaktimes['cci_42_active'] = breaktimes['cci_33_active']


        ranges = {'cci_22_combined': np.array(['1980-01-01', '2014-12-31']),
                  'cci_31_combined': np.array(['1980-01-01', '2015-12-31']),
                  'cci_31_passive': np.array(['1980-01-01', '2015-12-31']),
                  'cci_31_active': np.array(['1991-08-01', '2015-12-31']),
                  'cci_32_combined': np.array(['1980-01-01', '2015-12-31']),
                  'cci_33_combined': np.array(['1978-10-26', '2016-12-31']),
                  'cci_42_combined': np.array(['1978-10-26', '2016-12-31'])}

        if product not in timeframes.keys() or product not in breaktimes.keys() or product not in ranges.keys():
            raise Exception('No test times for selected cci product')

        self.timeframes = timeframes[product]
        self.breaktimes = breaktimes[product]
        self.ranges = ranges[product]

        self.grid_points = SMECV_Grid_v042().get_grid_points()


    @staticmethod
    def as_datetimes(datestrings):
        return [datetime.strptime(time, '%Y-%m-%d') for time in datestrings]

    @staticmethod
    def del_times(datestrings, skip_times):
        return np.array([datestrings[i] for i in range(datestrings.shape[0]) if i not in skip_times])

    def gpi_dep_times(self):
        if any([isinstance(x, dict) for x in self.breaktimes]) \
            or any([isinstance(x, dict) for x in self.timeframes]):
            return True
        else:
            return False

    def get_times(self,gpi=None, ignore_conditions=False, as_datetime=False, skip_times=None):
        '''
        Extract Times from data array like breaktime arrays or timeframe arrays
        :param data_array: self.timeframes or self.breaktimes
        :param ignore_conditions: set True to get the maximum of all possible times
        :return: dict
        '''
        return_times = {'ranges': self.ranges}
        if self.gpi_dep_times():
            if ignore_conditions:
                times = self.max_times(as_datetime=as_datetime)
                return_times['breaktimes'] = times['breaktimes']
                return_times['timeframes'] = times['timeframes']
            elif not gpi:
                raise Exception('Choose gpi or ignore_conditions')
            else:
                lat = self.grid_points[2][np.where(self.grid_points[0] == gpi)[0][0]]
                lon = self.grid_points[1][np.where(self.grid_points[0] == gpi)[0][0]]

                for name,timeset in {'timeframes':self.timeframes, 'breaktimes':self.breaktimes}.iteritems():
                    return_times[name] = []
                    for i, time in enumerate(timeset):
                        if isinstance(time, dict):
                            for condition, value in time.iteritems():
                                if ignore_conditions or eval(condition):
                                    for v in value:
                                        if v not in return_times[name]:
                                            return_times[name].append(v)
                        else:
                            return_times[name].append(time)

        else:
            return_times['timeframes'] = self.timeframes
            return_times['breaktimes'] = self.breaktimes

        if skip_times:
            return_times['timeframes'] = self.del_times(return_times['timeframes'], skip_times)
            return_times['breaktimes'] = self.del_times(return_times['breaktimes'], skip_times)

        if as_datetime:
            return_times['breaktimes'] = self.as_datetimes(return_times['breaktimes'])
            return_times['timeframes'] = [self.as_datetimes(timeframe) for timeframe in return_times['timeframes']]
            return_times['ranges'] = self.as_datetimes(return_times['ranges'])

        for name, data in return_times.iteritems():
            return_times[name] = np.array(data)

        return return_times


    def max_times(self, as_datetime=False):

        # Return the maximum of times (select the condition that contains most values)

        return_data = {'breaktimes':{}, 'timeframes':{}, 'ranges':self.ranges}
        for name, data in {'breaktimes': self.breaktimes, 'timeframes': self.timeframes}.iteritems():
            union =[]
            for i,d in enumerate(data):
                if isinstance(d,dict):
                    longest_condition_name = None
                    longest_condition_size = 0
                    for condition, values in d.iteritems():
                        if len(values) >= longest_condition_size:
                            longest_condition_name = condition
                            longest_condition_size = len(values)
                    for value in d[longest_condition_name]:
                        union.append(value)
                else:
                    union.append(d)
            return_data[name] = union

        return return_data



if __name__ == '__main__':
    ds = CCITimes('cci_42_combined')
    for gpi in [731807,730367]: # 37 deg border between gpis!! --> different times for cci33 etc
        times = ds.get_times(ignore_conditions=True, as_datetime=True)
        print 'breaktimes:'
        print times['breaktimes']
        print 'timeframes:'
        print times['timeframes']
        print 'ranges:'
        print times['ranges']


