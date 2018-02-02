# -*- coding: utf-8 -*-
"""
Created on Jän 03 13:56 2018

@author: wpreimes
"""
import os
import csv
import numpy as np

class dict_csv_wrapper(object):

    def __init__(self, data):
        '''
        Initialize the csv handler class

        :param data: dict or filepath
            Either takes a dictionary of the correct form or uses the input as a path
            to read an existing csv file of a correct form.
                dict: {'var1':[val1, val2,...], 'var2':[val1,val2,...], ...}

                file:
                    header: var1, var2, ...
                    data:   val1, val2, ...
                            ...   ...

        '''
        if isinstance(data, dict):
            self.content = data
        else:
            self.content = self._read(data) #type: dict
        self.header = sorted(self.content.keys()) #type: list

    @staticmethod
    def merge_dicts(*dict_args):
        result = None
        for dictionary in dict_args:
            if not result:
                result = dictionary
            else:
                for key, cont in dictionary.iteritems():
                    result[key] += cont
        return result

    def asint(self):
        for n, d in self.content.iteritems():
            self.content[n] = map(int, self.content[n])


    def asfloat(self):
        for n, d in self.content.iteritems():
            self.content[n] = map(float, self.content[n])

    def check_headers(self, other):
        for key in other.header:
            if key not in self.header:
                raise Exception("%s only exists in one file header" % key)
        return True

    def join_files(self, output_file=None, remove_old=True,  *others):
        for other in others:
            self.append(other)
            if remove_old:
                os.remove(other)
        if output_file:
            self.write(output_file)

        return self.content


    def _read(self, file_to_read):
        '''
        :param file_to_read: str
            full path to file that should be read
        :return: None
        '''
        with open(file_to_read, mode='r') as infile:
            reader = csv.reader(infile)
            header = None
            file_data = {}
            for i, row in enumerate(reader):
                if i == 0:
                    header = row
                    for var in header:
                        file_data[var] = []
                else:
                    for h, r in zip(header, row):
                        if r:
                            file_data[h].append(r)
        return file_data

    def append(self, other):
        other_file = dict_csv_wrapper(other)
        if self.check_headers(other_file):
            self.content = self.merge_dicts(self.content, other_file.content)

    def write(self, output_file):
        '''
        Write the content of a dictionary as csv file, keys are used as file header
        :param dicts_to_write: dict
            data to write, will be transformed to string
        :param output_file: string
            full path and filename for csv file to create
        :return: None
        '''
        with open(output_file, 'wb') as outfile:
            w = csv.writer(outfile, delimiter = ",")
            w.writerow(self.header)
            w.writerows(zip(*[self.content[key] for key in self.header]))


if __name__ == '__main__':
    file = dict_csv_wrapper({'var1':[1,2,3], 'var2':[1,2,3], 'var3': [1,2,3]})
    joined = file.join_files(r"C:\Temp\csv\merged.csv", True, r"C:\Temp\csv\test2.csv", r"C:\Temp\csv\test3.csv")



