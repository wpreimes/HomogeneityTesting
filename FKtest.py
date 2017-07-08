# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:51:58 2016

@author: wpreimes
"""


def FKtest(DF_X, mode='median', alpha=0, createlog=False):
    '''
    FKTEST Fligner-Killeen test for homogeneity of variances.
    
     Trujillo-Ortiz, A., R. Hernandez-Walls and N. Castro-Castro. (2009).FKtest:
       Fligner-Killeen test for homogeneity of variances. A MATLAB file. [WWW document].
       URL http://www.mathworks.com/matlabcentral/fileexchange/25040
    '''

    '''
    Format der Eingangsdaten:
        Dataframe: 2 Spalten
            Spalte1 (data): Data    Spalte2 (group): Messgruppenzahl (1 oder 2 für Referenzdaten bzw. Testdaten)
    '''
    # number of measurements and datagroups

    import numpy as np
    import scipy.stats as stats

    DF_X = DF_X.rename(columns={DF_X.ix[:, 0].name: 'data'})
    DF_X = DF_X.dropna()
    N = DF_X.index.size
    K = DF_X['group'].nunique()

    DF_X['A'] = np.nan

    if mode == 'median':
        for i in range(K):
            subset = DF_X.ix[DF_X['group'] == i]
            groupmed = subset.data.median()
            DF_X.ix[DF_X['group'] == i, 'groupmed'] = groupmed  # groupmedians
            DF_X.ix[DF_X['group'] == i, 'groupme_diff'] = np.abs(
                subset['data'] - groupmed)  # difference data-groupmedians

    if mode == 'mean':
        for i in range(K):
            subset = DF_X.ix[DF_X['group'] == i]
            groupmean = subset.data.mean()
            DF_X.ix[DF_X['group'] == i, 'groupmean'] = groupmean  # groupmeans
            DF_X.ix[DF_X['group'] == i, 'groupme_diff'] = np.abs(
                subset['data'] - groupmean)  # difference data-groupmeans

    Z = stats.rankdata(DF_X['groupme_diff'])  # score ranking ALL
    sta_norm_dist = stats.norm.ppf(0.5 + (Z / (2. * (N + 1.))))  # score standard normal distribution ALL
    DF_X['A'] = sta_norm_dist
    M = DF_X['A'].mean()  # overall mean

    nn = []
    mm = []
    bb = []
    for i in range(K):
        subset = DF_X.ix[DF_X['group'] == i]

        nn.append(subset.index.size)
        mm.append(np.mean(subset['A']))
        DF_X.ix[DF_X['group'] == i, 'groupAmean'] = mm[i]
        bb.append((nn[i] * (mm[i] - M) ** 2))
        DF_X.ix[DF_X['group'] == i, 'groupB'] = bb[i]

    B = np.array(DF_X['groupB'].unique())
    V = DF_X['A'].var()  # Overall Variance Score
    X2 = np.sum(B) / V  # Fligner-Killeen statistic by the Chi-squared approximation
    v = K - 1  # statistic degree of freedom
    F = (X2 / v) / ((N - 1. - X2) / (N - K))  # Fligner-Killeen statistic by the Fisher approximation

    P1 = 1 - stats.chi2.cdf(X2, v)
    P2 = 1 - stats.f.cdf(F, v, N - K)

    # TODO: Laut Chun-Hsu statt F X2??
    stats = {'chi': {'z': X2, 'df': v, 'pval': P1}, 'f': {'z': F, 'df': [v, N - K], 'pval': P2}}

    if stats['chi']['pval'] < alpha:
        h = 1
    else:
        h = 0

    return h, stats
