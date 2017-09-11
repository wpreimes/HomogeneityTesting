# -*- coding: utf-8 -*-
"""
Created on Jul 28 16:03 2017

@author: wpreimes
"""
import numpy as np
from datetime import datetime
from smecv_grid.grid import SMECV_Grid_v042
from otherfunctions import cci_extract, cci_string_combine


# TODO: Add timeframes for active and passive products for 22 and 33
# TODO: Check timeframe and breaktimes for CCI33

# Order matters!!!!
class CCITimes(object):
    def __init__(self, product, ignore_position=True, skip_times=None):

        timeframes = {'CCI_22_COMBINED': np.array(
            [['2011-10-01', '2014-12-31'], ['2007-01-01', '2012-07-01'],
             ['2002-07-01', '2011-10-01'], ['1998-01-01', '2007-01-01'],
             ['1991-08-01', '2002-07-01'], ['1987-07-01', '1998-01-01']]),
            'CCI_31_COMBINED': np.array(
                [['2011-10-01', '2015-05-01'], ['2010-07-01', '2012-07-01'],
                 ['2007-10-01', '2011-10-01'], ['2007-01-01', '2010-07-01'],
                 ['2002-07-01', '2007-10-01'], ['1998-01-01', '2007-01-01'],
                 ['1991-08-01', '2002-07-01'], ['1987-09-01', '1998-01-01']]),
            'CCI_31_PASSIVE': np.array(
                [['2012-07-01', '2015-12-31'], ['2011-10-01', '2015-05-01'],
                 ['2010-07-01', '2012-07-01'], ['2007-10-01', '2011-10-01'],
                 ['2002-07-01', '2010-07-01'], ['1998-01-01', '2007-10-01'],
                 ['1987-09-01', '2002-07-01']]),
            'CCI_31_ACTIVE': np.array(
                [['2007-01-01', '2015-12-01'], ['2003-02-01', '2012-11-01'],
                 ['1997-05-01', '2007-01-01'], ['1991-08-01', '2003-02-01']]),
            'CCI_33_ACTIVE': np.array(
                [['1991-07-01', '2016-12-31']]),
            'CCI_33_COMBINED': np.array(
                [['2011-10-05', '2016-12-31'], ['2010-01-15', '2012-07-01'],
                 ['2007-10-01', '2011-10-05'], ['2007-01-01', '2010-01-15'],
                 ['2002-06-19', '2007-10-01'],
                 {'lat < -37 or lat > 37':
                      [['1991-08-05', '2007-01-01'],
                       ['1987-07-09', '2002-06-19']],
                  'lat >= -37 and lat <= 37':
                      [['1998-01-01', '2007-01-01'],
                       ['1991-08-05', '2002-06-19'],
                       ['1987-07-09', '1998-01-01']]
                  },
                 ['1978-10-26', '1991-08-05']]),
            'CCI_33_PASSIVE': np.array(
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
        timeframes['CCI_41_COMBINED'] = timeframes['CCI_33_COMBINED']
        timeframes['CCI_41_PASSIVE'] = timeframes['CCI_33_PASSIVE']
        timeframes['CCI_41_ACTIVE'] = timeframes['CCI_33_ACTIVE']

        breaktimes = {'CCI_22_COMBINED': np.array(
            ['2012-07-01', '2011-10-01', '2007-01-01', '2002-07-01', '1998-01-01',
             '1991-08-01']),
            'CCI_31_COMBINED': np.array(
                ['2012-07-01', '2011-10-01', '2010-07-01', '2007-10-01', '2007-01-01',
                 '2002-07-01', '1998-01-01', '1991-08-01']),
            'CCI_31_PASSIVE': np.array(
                ['2015-05-01', '2012-07-01', '2011-10-01', '2010-07-01', '2007-10-01',
                 '2002-07-01', '1998-01-01']),
            'CCI_31_ACTIVE': np.array(['2012-11-01', '2007-01-01', '2003-02-01', '1997-05-01']),
            'CCI_33_ACTIVE': np.array(['2007-01-01']),
            'CCI_33_COMBINED': np.array(
                ['2012-07-01', '2011-10-05', '2010-01-15', '2007-10-01',
                 '2007-01-01', '2002-06-19',
                 {'lat >= -37 and lat <= 37': ['1998-01-01']},
                 '1991-08-05', '1987-07-09']),
            'CCI_33_PASSIVE': np.array(
                ['2012-07-01', '2011-10-05', '2010-01-15', '2007-10-01', '2002-06-19',
                 {'lat >= -37 and lat <= 37': ['1998-01-01']},
                 '1987-07-09'])
        }
        breaktimes['CCI_41_COMBINED'] = breaktimes['CCI_33_COMBINED']
        breaktimes['CCI_41_PASSIVE'] = breaktimes['CCI_33_PASSIVE']
        breaktimes['CCI_41_ACTIVE'] = breaktimes['CCI_33_ACTIVE']

        ranges = {'CCI_22_COMBINED': np.array(['1980-01-01', '2014-12-31']),
                  'CCI_31_COMBINED': np.array(['1980-01-01', '2015-12-31']),
                  'CCI_31_PASSIVE': np.array(['1980-01-01', '2015-12-31']),
                  'CCI_31_ACTIVE': np.array(['1991-08-01', '2015-12-31']),
                  'CCI_32_COMBINED': np.array(['1980-01-01', '2015-12-31']),
                  'CCI_33_COMBINED': np.array(['1978-10-26', '2016-12-31']),
                  'CCI_33_ACTIVE': np.array(['1991-08-05', '2016-12-31']),
                  'CCI_33_PASSIVE': np.array(['1978-10-26', '2016-12-31']),
                  'CCI_41_COMBINED': np.array(['1978-10-26', '2016-12-31']),
                  'CCI_41_ACTIVE': np.array(['1991-08-05', '2016-12-31']),
                  'CCI_41_PASSIVE': np.array(['1978-10-26', '2016-12-31'])}

        if product not in timeframes.keys() or product not in breaktimes.keys() or product not in ranges.keys():
            self.product = self._init_timeframes_for_adjusted(product)
        else:
            self.product = product

        self.skip_times = skip_times

        self.timeframes = timeframes[self.product]
        self.breaktimes = breaktimes[self.product]
        self.ranges = ranges[self.product]

        if not self.gpi_dep_times():
            self.ignore_position = True
        else:
            self.ignore_position = ignore_position

        if not self.ignore_position:
            self.grid_points = SMECV_Grid_v042().get_grid_points()

        if self.skip_times and not self.ignore_position and self.gpi_dep_times():
            raise Exception('Skipping breaktimes for gpi dependent timeframes ambiguous')

    @staticmethod
    def as_string(datetimes):
        '''
        Datetime objects to strings
        :param datetimes: np.array
        :return: np.array
        '''
        return [str(time.date()) for time in datetimes]

    @staticmethod
    def as_datetimes(datestrings):
        '''
        Turn strings to datetime objects
        :param datestrings: np.array
        :return: np.array
        '''
        return [datetime.strptime(time, '%Y-%m-%d') for time in datestrings]

    @staticmethod
    def del_times(datestrings, skip_times):
        '''
        returns input without elements selected by index in skipt_times
        :param datestrings: np.array
        :param skip_times: list
        :return: np.array
        '''
        return [datestrings[i] for i in range(len(datestrings)) if i not in skip_times]

    def _init_timeframes_for_adjusted(self, product):
        info = cci_extract(product)
        if info['adjust']:
            return cci_string_combine(info)
        else:
            raise Exception('No test times for selected cci product')

    def gpi_dep_times(self):
        '''
        Checks if current breaktimes and timeframes are positional dependent
        :return: bool
        '''
        if any([isinstance(x, dict) for x in self.breaktimes]) \
                or any([isinstance(x, dict) for x in self.timeframes]):
            return True
        else:
            return False

    def max_times(self):
        # Return the maximum of times (select the condition that contains most values)
        return_data = {'breaktimes': {}, 'timeframes': {}, 'ranges': self.ranges}
        for name, data in {'breaktimes': self.breaktimes, 'timeframes': self.timeframes}.iteritems():
            union = []
            for i, d in enumerate(data):
                if isinstance(d, dict):
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

    def get_times(self, gpi=None, as_datetime=False):
        '''
        Extract Times from data array like breaktime arrays or timeframe arrays
        :param data_array: self.timeframes or self.breaktimes
        :param ignore_conditions: set True to get the maximum of all possible times
        :return: dict
        '''
        return_times = {'ranges': self.ranges}

        if not self.gpi_dep_times:
            return_times['timeframes'] = self.timeframes
            return_times['breaktimes'] = self.breaktimes
        else:
            if not gpi or self.ignore_position:
                times = self.max_times()
                return_times['breaktimes'] = times['breaktimes']
                return_times['timeframes'] = times['timeframes']
            else:
                lat = self.grid_points[2][np.where(self.grid_points[0] == gpi)[0][0]]
                lon = self.grid_points[1][np.where(self.grid_points[0] == gpi)[0][0]]

                for name, timeset in {'timeframes': self.timeframes, 'breaktimes': self.breaktimes}.iteritems():
                    return_times[name] = []
                    for i, time in enumerate(timeset):
                        if isinstance(time, dict):
                            for condition, value in time.iteritems():
                                if eval(condition):
                                    for v in value:
                                        if v not in return_times[name]:
                                            return_times[name].append(v)
                        else:
                            return_times[name].append(time)

        if self.skip_times:
            return_times['timeframes'] = self.del_times(return_times['timeframes'], self.skip_times)
            return_times['breaktimes'] = self.del_times(return_times['breaktimes'], self.skip_times)

        if as_datetime:
            return_times['breaktimes'] = self.as_datetimes(return_times['breaktimes'])
            return_times['timeframes'] = [self.as_datetimes(timeframe) for timeframe in return_times['timeframes']]
            return_times['ranges'] = self.as_datetimes(return_times['ranges'])

        for name, data in return_times.iteritems():
            return_times[name] = np.array(data)

        return return_times

    def get_index(self, gpi, time, gpi_times=None):
        if type(time) is np.ndarray:
            if all([isinstance(t, datetime) for t in time]):
                time = self.as_string(time)
            if gpi_times:
                times = gpi_times['timeframes']
            else:
                times = self.get_times(gpi, as_datetime=False)['timeframes']
            return np.where((times == time)[:, 1])[0][0]
        else:
            if isinstance(time, datetime):
                time = self.as_string([time])[0]
            if gpi_times:
                times = gpi_times['breaktimes']
            else:
                times = self.get_times(gpi, as_datetime=False)['breaktimes']
            return np.where(times == time)[0][0]

    def breaktime_for_timeframe(self, gpi, timeframe):
        if all([isinstance(time, datetime) for time in timeframe]):
            return self.as_datetimes([self.breaktimes[self.get_index(gpi, np.array(self.as_string(timeframe)))]])[0]
        else:
            return self.breaktimes[self.get_index(gpi, timeframe)]

    def timeframe_for_breaktime(self, gpi, breaktime):
        if isinstance(breaktime, datetime):
            return self.as_datetimes(self.timeframes[self.get_index(gpi, breaktime)])
        else:
            return self.timeframes[self.get_index(gpi, breaktime)]

    def get_adjacent(self, gpi, reference, shift):
        '''
        :param timeframe: the initial position/timeframe
        :param shift: index change for timeframe to return, eg 1 -> return next timeframe, -1 return previous timeframe
        :return: the selected timeframe in same format as input timeframe
        '''
        times = self.get_times(gpi, as_datetime=False)
        index = self.get_index(gpi, reference, times)
        if type(reference) is np.ndarray:
            if all([isinstance(time, datetime) for time in reference]):
                return self.as_datetimes(times['timeframes'][index + shift])
            else:
                return times['timeframes'][index + shift]
        else:
            if isinstance(reference, datetime):
                return self.as_datetimes([times['breaktimes'][index + shift]])[0]
            else:
                return times['breaktimes'][index + shift]


if __name__ == '__main__':
    ds = CCITimes('CCI_33_COMBINED', ignore_position=True, skip_times=[0,1,2,3,4,5,6,8])
    for gpi in [620964, 805227]:  # First Point < 37 Lat, second > 37 Lat
        times = ds.get_times(gpi, as_datetime=True)
        print 'breaktimes:'
        print times['breaktimes']
        print 'timeframes:'
        print times['timeframes']
        print 'ranges:'
        print times['ranges']
        print (ds.timeframe_for_breaktime(gpi, times['breaktimes'][0]))
        print (ds.breaktime_for_timeframe(gpi, times['timeframes'][0]))
        #print (ds.get_adjacent(gpi, times['timeframes'][1], 2))
        #print (ds.get_adjacent(gpi, times['breaktimes'][5], -2))
    ds = CCITimes('CCI_33_COMBINED', ignore_position=True)
    for gpi in [620964, 805227]:  # First Point < 37 Lat, second > 37 Lat
        times = ds.get_times(gpi, as_datetime=False)
        print 'breaktimes:'
        print times['breaktimes']
        print 'timeframes:'
        print times['timeframes']
        print 'ranges:'
        print times['ranges']
        print (ds.timeframe_for_breaktime(gpi, times['breaktimes'][0]))
        print (ds.breaktime_for_timeframe(gpi, times['timeframes'][0]))
        #print (ds.get_adjacent(gpi, times['timeframes'][1], 2))
        #print (ds.get_adjacent(gpi, times['breaktimes'][5], -2))
