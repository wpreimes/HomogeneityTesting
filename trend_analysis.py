# -*- coding: utf-8 -*-
"""
Created on Sep 16 11:33 2017

@author: wpreimes
"""

from scipy.stats import theilslopes
from interface import BreakTestData
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates

class Trends(BreakTestData):
    def __init__(self, test_prod, ref_prod, anomaly, use_adjusted=False):
        BreakTestData.__init__(self, test_prod, ref_prod, anomaly)
    def TheilSenTrend(self, y, alpha=0.95):
        return theilslopes(y=y,alpha=alpha)


if __name__ == "__main__":
    Trends_obj = Trends('CCI_41_COMBINED_ADJUSTED', 'merra2', False, False)
    df = Trends_obj.read_gpi(702473)['1990-01-01':'2010-12-31']
    fig, axs = plt.subplots(2, sharex=False)
    for i, data in enumerate(['testdata','refdata']):
        df_sub = df[[data]]
        df_sub = Trends_obj.temp_resample(df_sub, 'M', 0.33).dropna()
        y = df_sub[data].dropna()
        res = Trends_obj.TheilSenTrend(y.values)
        medslope = res[0]
        inter = res[1]
        lowlope = res[2]
        upslope = res[3]


        x = np.array(range(y.size))
        axs[i].plot(y.index, res[1] + res[0] * x, 'r-')
        axs[i].plot(y.index, y, 'b-')
        axs[i].set_title(data)

    plt.tight_layout()
    plt.show()