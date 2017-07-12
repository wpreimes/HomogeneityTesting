# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 13:28:09 2017

@author: wpreimes
"""
import pandas as pd
import numpy as np
from points_to_netcdf import points_to_netcdf
import os

class save_Log(object):
    def __init__(self, workfolder, test_prod, ref_prod, anomaly, cells):
        # type: (str, str, str, bool, list) -> None
        self.workfolder = workfolder

        with open(os.path.join(workfolder, 'log.txt'), 'w') as f:
            f.write('Log file for HomogeneityTesting \n')
            f.write('=====================================\n')
            f.write('Test Product: %s \n' % test_prod)
            f.write('Reference Product: %s \n' % ref_prod)
            f.write('Anomaly data: %s \n' % str(anomaly))
            f.write('Processed Cells: %s \n' % cells)
            f.write('=====================================\n')
            f.write('\n')

    def add_line(self, string):
        # type: (str) -> None
        # Initialize and add line to log file for current process
        with open(os.path.join(self.workfolder, 'log.txt'), 'a') as f:
            f.write(string + '\n')

class load_Log(object):
    def __init__(self, workfolder):
        self.workfolder = workfolder

    def get_products(self):
        with open(os.path.join(self.workfolder,'log.txt'), mode='r') as f:
            lines = [line.rstrip('\n') for line in f]

        test_prod = lines[2].split(' ')[2]
        ref_prod = lines[3].split(' ')[2]

        return {'test_prod': test_prod, 'ref_prod': ref_prod}



class SaveResults(object):
    def __init__(self, path, grid, filename, buffer_size=200):
        # type: (str,object,str,int) -> None
        self.total_number = [0, ]
        self.filename = filename
        self.path = path
        self.grid = grid
        self.buffer_size = buffer_size
        self.data = {'gpi': [], 'h_fk': [], 'h_wk': [], 'test_results': [], 'status': []}

    def reset_data(self):
        self.data = {'gpi': [], 'h_fk': [], 'h_wk': [], 'test_results': [], 'status': []}

    def extract_infos(self, data_dict):
        # type: (dict) -> dict
        status = int(data_dict['status'][0])
        if status in range(3, 10):
            return {'h_wk': np.nan, 'h_fk': np.nan, 'test_results': np.nan, 'status': status}
        else:
            wk = data_dict['wilkoxon']['h']
            fk = data_dict['fligner_killeen']['h']

            if type(wk) is str:
                # wk test failed and contains error message
                wk = np.nan
                status = 1
            if type(fk) is str:
                # fk failed and contains error message
                fk = np.nan
                status = 2

            all = 4.0
            if wk == 1:
                if fk == 1:
                    all = 3.0
                else:
                    all = 1.0
            elif fk == 1:
                all = 2.0

            return {'h_wk': wk, 'h_fk': fk, 'test_results': all, 'status': status}

    def fill_buffer(self, gpi, data_dict):
        # type: (int,dict) -> None
        data_dict = self.extract_infos(data_dict)
        self.data['gpi'].append(gpi)
        for name, val in data_dict.iteritems():
            self.data[name].append(val)

        if len(self.data['gpi']) == self.buffer_size:
            self.save_to_netcdf()

    def save_to_netcdf(self):
        var_meta_dicts = {
            'test_results': {
                'Description': 'Homogeneity Test Results Classified',
                'Values': '1 = WK only, 2 = FK only, 3 = WK and FK, 4 = None'},
            'status': {
                'Description': 'Homogeneity Testing Status',
                'Values': '0 = Processing OK, \
                         1 = Error during WK testing,\
                         2 = Error during FK testing,\
                         3 = No coinciding data for timeframe,\
                         4 = Test TS and Ref TS do not match,\
                         5 = Spearman correlation failed,\
                         6 = Minimum Dataseries Length not reached,\
                         7 = Negative or NaN correlation,\
                         8 = WK test and FK test failed,\
                         9 = Error reading gpi'
            }
        }
        df = pd.DataFrame.from_dict(self.data)
        df = df.set_index('gpi')

        points_to_netcdf(dataframe=df, path=self.path,
                         filename=self.filename,
                         var_meta_dicts=var_meta_dicts)

        self.reset_data()


'''     
test_obj.DF_Points['test_results']=np.nan
test_obj.DF_Points['message']='Not processed'


test_obj.save_test_results(gpi,'h_fk',testresult['FlignerKilleen']['h'])
test_obj.save_test_results.set_value(gpi,'h_wk',testresult['Wilkoxon']['h'])

test_obj.save_test_results.set_value(gpi,'message','Processing OK')  
        
#In case something went wrong
test_obj.DF_Points.set_value(gpi,'h_fk',np.nan)
test_obj.DF_Points.set_value(gpi,'h_wk',np.nan)
test_obj.DF_Points.set_value(gpi,'message',str(e))
'''
