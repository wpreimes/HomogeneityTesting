# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 10:42:39 2016

@author: wpreimes
"""

import numpy as np
import pandas as pd
from rsdata.ESA_CCI_SM.interface import ESA_CCI_SM
from pytesmo.timedate.julian import caldat

from gldas.interface import GLDASTs
from merra.interface import MERRA2_Ts
from pynetcf.time_series import GriddedNcIndexedRaggedTs
import os
from datetime import datetime
from otherfunctions import merge_ts, regress, cci_extract, cci_string_combine
from read_adjusted_ts import cciAdjustedTs
import re
from rsroot import root_path


def read_warp_ssm(ssm, ssf, gpi):
    import pytesmo.timedate.julian as julian

    ssm_ts = ssm.read_ts(gpi)
    ssf_ts = ssf.read_ts(gpi)

    # ceil or round?
    jd = np.round(ssm_ts['jd'] - 0.5) + 0.5
    # jd = ssm_ts['jd']
    ts = pd.DataFrame(ssm_ts, index=julian.julian2datetimeindex(jd))
    ts['ssf'] = ssf_ts['ssf']

    ts = ts[(ts['proc_flag'] <= 2) & (ts['ssf'] == 1)]['sm']
    ts.index.tz = None

    return ts.groupby(level=0).last()


class QDEGdata_M(object):
    def __init__(self, products):
        self.monthproducts = {}
        m_otherproduts = []
        self.m_otherdata = None

        # Data that must be resampled
        if 'gldas_v2' in products:
            m_otherproduts.append('gldas_v2')
        if 'gldas_v1' in products:
            m_otherproduts.append('gldas_v1')
        if 'gldas-merged' in products:
            m_otherproduts.append('gldas-merged')
        if 'trmm' in products:
            m_otherproduts.append('trmm')
        if 'merra2' in products:
            m_otherproduts.append('merra2')

        cci_re = re.compile("cci_.+", re.IGNORECASE)
        if any([cci_re.match(product) for product in products]):
            for cci_product in [cci_re.match(product) for product in products]:
                if cci_product:
                    cci_product = cci_product.group()
                    m_otherproduts.append(cci_product)
        if m_otherproduts:
            self.m_otherdata = QDEGdata_D(products=m_otherproduts, only_sm=True)


        if 'gldas-merged-from-file' in products:  # Monthly merged GLDAS Data for faster loading
            if os.path.isdir(r"D:\USERS\wpreimes\datasets\GLDAS-NOAH\GLDAS_NOAH_merged_M_at0"):
                print('Found local files for gldas-merged at 0H')
                path_gldas_monthly = r"D:\USERS\wpreimes\datasets\GLDAS-NOAH\GLDAS_NOAH_merged_M_at0"
            else:
                path_gldas_monthly = r"U:\datasets\GLDAS-NOAH\GLDAS_NOAH_merged_M_at0"
            self.monthproducts['gldas_merged'] = GLDASTs(path_gldas_monthly)

    def read_gpi(self, gpi, startdate, enddate):

        # Read sub-monthly data
        if isinstance(startdate, datetime):
            startdate = str(startdate)
        if isinstance(enddate, datetime):
            enddate = str(enddate)

        if not self.m_otherdata:
            data_group = pd.DataFrame()
        else:
            df_day = self.m_otherdata.read_gpi(gpi, startdate, enddate)
            df_month = df_day.resample('M').mean()
            data_group = df_month

        # Read monthly data
        if 'gldas-merged-from-file' in self.monthproducts.keys():
            try:
                ts_gldas_merged = self.monthproducts['gldas-merged-from-file'].read_ts(gpi)[['SoilMoi0_10cm_inst']]
            except:
                ts_gldas_merged = pd.Series(index=pd.date_range(start=startdate, end=enddate))

            if (len(np.where(~np.isnan(ts_gldas_merged))[0]) == 0):
                # raise Exception, 'No GLDAS Data available for GPI %i' %gpi
                #print 'No merged gldas data for gpi %0i' % gpi
                pass
            else:
                ts_gldas_merged = ts_gldas_merged.rename(columns={'SoilMoi0_10cm_inst': 'gldas-merged-from-file'})

                '''          
                ts_gldas_merged_for_mean=pd.DataFrame(index=ts_gldas_merged.index.to_datetime(),
                                                   data=ts_gldas_merged)
                ts_gldas_merged_for_mean=ts_gldas_merged_for_mean.resample('M',how='mean')
                '''
            if data_group.empty:
                data_group = ts_gldas_merged
            else:
                data_group = pd.concat([data_group, ts_gldas_merged], axis=1)

        return data_group[startdate:enddate]


class QDEGdata_D(object):
    # changed path in "D:\workspace\smecv\grids\ecv.py" to grid file
    # changed path definition in D:\workspace\pygrids\GLDAS_NOAH.py
    # change path in rsdata-GLDAS-interface
    # Requires ECV_CCI_gridv4.nc

    def __init__(self, products, resample='mean_over_day', only_sm=True):

        self.resample = resample
        self.only_sm = only_sm

        self.dayproducts = {}

        d_otherproducts = []
        self.d_otherdata = None

        # Add 3H products
        if 'gldas_v1' in products:
            d_otherproducts.append('gldas_v1')
        if 'gldas_v2' in products:
            d_otherproducts.append('gldas_v2')
        if 'gldas_v21' in products:
            d_otherproducts.append('gldas_v21')
        if 'gldas-merged' in products:
            d_otherproducts.append('gldas-merged')

        if d_otherproducts:
            self.d_otherdata = QDEGdata_3H(products=d_otherproducts)

        # Add Daily Products
        # TODO: before adding versions: Add cfg file and add data pathes in cfgfile
        cci_versions = ['22', '31', '33', '41']
        cci_types = ['COMBINED', 'ACTIVE', 'PASSIVE']

        cci_res = [re.compile("cci_%s.+" % version, re.IGNORECASE) for version in cci_versions]
        for cci_re in cci_res:
            for cci_product in [cci_re.match(product) for product in products]:
                if cci_product:
                    cci_product = cci_product.group()
                    info = cci_extract(cci_product)
                    cci_product = cci_string_combine(info)
                    if info['prefix'] != 'CCI' or info['version'] not in cci_versions or info['type'] not in cci_types:
                        raise Exception('cci version or product not known...use format cci_XX_PRODUCT')

                    if os.path.isdir(os.path.join(root_path.d, 'USERS', 'wpreimes', 'datasets', 'CCI_%s_D' % info['version'])):
                        print('Found local files for %s daily cci_%s %s data' % (info['type'], info['version'],
                                                                                 info['adjust'] if info['adjust'] else 'UNADJUSTED'))
                    else:
                        print ('Use CCI data from R')

                    cfg_path = os.path.join(root_path.d, 'USERS', 'wpreimes', 'datasets', 'HomogeneityTesting_data',
                                            'cci_cfg_local')
                    if not os.path.isdir(cfg_path):
                        cfg_path = os.path.join(root_path.h, 'HomogeneityTesting_data', 'cci_cfg_local')


                    cci = ESA_CCI_SM('ESA_CCI_SM_v0%s.%s' % (info['version'][0], info['version'][1]),
                                     parameter=info['type'],
                                     cfg_path=cfg_path)
                    self.dayproducts[cci_product] = cci

        if 'merra2' in products:
            lu_table_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'merra_gpi_LUT.csv')
            self.lkup = pd.read_csv(lu_table_file, index_col=0)

            path_merra2 = os.path.join(root_path.d, 'USERS', 'wpreimes', 'datasets', 'MERRA2_D', 'ts_daily_0030')
            if os.path.isdir(path_merra2):
                print('Found local files for daily merra2')
            else:
                path_merra2 = os.path.join(root_path.r, 'Datapool_processed', 'Earth2Observe', 'MERRA2',
                                           'M2T1NXLND.5.12.4', 'datasets', 'ts_daily_0030')
                print('Found files for daily merra2 on R')

            self.dayproducts['merra2'] = MERRA2_Ts(path_merra2)


    def read_gpi(self, gpi, startdate, enddate):
        """
        Read one or multiple products for selected ground point in QDEG grid
        and return SM values as dataframe from a defined start date to a
        defined end date.
        
        Parameters
        ----------
        gpi : int
            index of point in QDEG grid
            
        startdate,enddate : strings of form: 'yyyy-mm-dd'
            start and end date of sm timeseries
            
        products : *list
            one or multiple products ('gldas', 'trmm')
        """
        if isinstance(startdate, datetime):
            startdate = str(startdate)
        if isinstance(enddate, datetime):
            enddate = str(enddate)
        # Read and resample hour data
        if not self.d_otherdata:
            data_group = pd.DataFrame(index=pd.date_range(startdate, enddate, freq='D'))
        else:
            df_hour = self.d_otherdata.read_gpi(gpi, startdate, enddate)
            if self.resample == 'mean_over_day':
                # Resample over all values of a day
                df_hour = df_hour.resample('D').mean()
            else:
                # Resample only over the vales of a certain time for each day
                df_hour = df_hour.at_time(self.resample)
                df_hour = df_hour.resample('D').mean()
            data_group = df_hour

        # read cci data (all versions in product)
        for name in [cci_name for cci_name in self.dayproducts.keys() if 'CCI' in cci_name]:
            cci_product = self.dayproducts[name]
            try:
                if isinstance(cci_product, GriddedNcIndexedRaggedTs):
                    cci_product = cciAdjustedTs(cci_product.path,os.path.join(cci_product.path, 'ECV_CCI_gridv4.nc'))
                    df_cci = pd.DataFrame(cci_product.read(gpi))
                    if df_cci.columns.values[0]!=name:
                        df_cci=df_cci.rename(columns={df_cci.columns.values[0]:name})
                    df_cci = df_cci*100
                else:
                    df_cci = pd.DataFrame(cci_product.read(gpi))
                    df_cci = df_cci.set_index('jd')
                    m, d, y = caldat(df_cci.index.values)
                    df_cci.index = pd.to_datetime(pd.DataFrame(data={'year': y, 'month': m, 'day': d}))
                    df_cci = df_cci[startdate:enddate]
                    df_cci = df_cci[df_cci['flag'] == 0]
                    df_cci['sm'].loc[df_cci['sm'] == -999999.] = np.nan
                    df_cci=df_cci.rename(columns={'sm': name})
            except:
                df_cci[name] = pd.Series(index=pd.date_range(start=startdate, end=enddate))

            if df_cci[name].isnull().all():
                # print 'No cci sm data for gpi %0i' % gpi
                pass
            if data_group.empty:
                if self.only_sm:
                    data_group['%s' % name] = df_cci[name]
                else:
                    data_group[name + '_' + df_cci.columns.values] = df_cci
            else:
                if self.only_sm:
                    data_group = pd.concat([data_group, df_cci[name]], axis=1)
                else:
                    data_group = pd.concat([data_group, df_cci], axis=1)

        if 'merra2' in self.dayproducts.keys():
            try:
                gpi_merra = self.lkup.loc[gpi].gpi_merra
                ts_merra2 = self.dayproducts['merra2'].read(gpi_merra)['GWETTOP']
                ts_merra2 = ts_merra2.resample('D').mean()
                ts_merra2 = ts_merra2 * 100
            except:
                ts_merra2 = pd.Series(index=pd.date_range(start=startdate, end=enddate))
                ts_merra2 = ts.resample('D').mean()

            if ts_merra2.isnull().all():
                # print 'No merra2 data for gpi %i' % gpi
                pass
            ts_merra2.index = ts_merra2.index.to_datetime().date
            ts_merra2.index = ts_merra2.index.to_datetime()

            if data_group.empty:
                data_group['merra2'] = ts_merra2
            else:
                data_group = pd.concat([data_group, ts_merra2.rename('merra2')], axis=1)

        if 'trmm' in self.dayproducts.keys():
            try:
                ts_trmm = self.dayproducts['trmm'].read_gp(gpi)['p'][startdate:enddate]
                ts_trmm[np.isnan(ts_trmm)] = 0
            except:
                ts_trmm = pd.Series(index=pd.date_range(start=startdate, end=enddate))
            if ts_trmm.isnull().all():
                # print 'No trmm data for gpi %0i' % gpi
                pass
            ts_trmm.index = ts_trmm.index.to_datetime().date

            if data_group.empty:
                data_group['trmm'] = ts_trmm
            else:
                data_group = pd.concat([data_group, ts_trmm.rename('trmm')], axis=1)

        return data_group[startdate:enddate]


class QDEGdata_3H(object):
    def __init__(self, products):
        self.H3products = {}

        if 'gldas_v1' in products:
            if os.path.isdir(r'D:\USERS\wpreimes\datasets\gldas_v1'):
                print('Found local files for GLDAS1 3H data')
                path_gldasv1 = r'D:\USERS\wpreimes\datasets\gldas_v1'
            else:
                path_gldasv1 = r"R:\Datapool_processed\GLDAS\GLDAS_NOAH025SUBP_3H\datasets\netcdf"
            self.H3products['gldas_v1'] = GLDASTs(path_gldasv1)

        if 'gldas_v2' or 'gldas-merged' in products:
            if os.path.isdir(r'D:\USERS\wpreimes\datasets\gldas_v2'):
                print('Found local files for GLDAS2 3H data')
                path_gldasv2 = r'D:\USERS\wpreimes\datasets\gldas_v2'
            else:
                path_gldasv2 = r"R:\Datapool_processed\GLDAS\GLDAS_NOAH025_3H.020\datasets\netcdf"
                self.H3products['gldas_v2'] = GLDASTs(path_gldasv2)

        if 'gldas_v21' or 'gldas-merged' in products:
            if os.path.isdir(r'D:\USERS\wpreimes\datasets\gldas_v21'):
                print('Found local files for GLDAS21 3H data')
                path_gldasv21 = r'D:\USERS\wpreimes\datasets\gldas_v21'
            else:
                path_gldasv21 = r"R:\Datapool_processed\GLDAS\GLDAS_NOAH025_3H.2.1\datasets"
                self.H3products['gldas_v21'] = GLDASTs(path_gldasv21)

    def read_gpi(self, gpi, startdate, enddate):
        if isinstance(startdate, datetime):
            startdate = str(startdate)
        if isinstance(enddate, datetime):
            enddate = str(enddate)
        data_group = pd.DataFrame()

        if 'gldas_v1' in self.H3products.keys():
            try:
                ts_gldas_v1 = self.H3products['gldas_v1'].read_ts(gpi)['086_L1']
            except:
                ts_gldas_v1 = pd.Series(index=pd.date_range(start=startdate, end=enddate))

            if (len(np.where(~np.isnan(ts_gldas_v1))[0]) == 0):
                # raise Exception, 'No GLDAS Data available for GPI %i' %gpi
                # print 'No gldas v1 data for gpi %0i' % gpi
                pass
            if data_group.empty:
                data_group['gldas_v1'] = ts_gldas_v1
            else:
                data_group = pd.concat([data_group, ts_gldas_v1.rename('gldas_v1')], axis=1)

        if 'gldas_v2' in self.H3products.keys():
            try:
                ts_gldas_v2 = self.H3products['gldas_v2'].read_ts(gpi)['086_L1']
            except:
                ts_gldas_v2 = pd.Series(index=pd.date_range(start=startdate, end=enddate))

            if (len(np.where(~np.isnan(ts_gldas_v2))[0]) == 0):
                # raise Exception, 'No GLDAS Data available for GPI %i' %gpi
                # print 'No gldas v1 data for gpi %0i' % gpi
                pass
            if data_group.empty:
                data_group['gldas_v2'] = ts_gldas_v2
            else:
                data_group = pd.concat([data_group, ts_gldas_v2.rename('gldas_v2')], axis=1)

        if 'gldas_v21' in self.H3products.keys():
            try:
                ts_gldas_v21 = self.H3products['gldas_v21'].read_ts(gpi)['SoilMoi0_10cm_inst']
            except:
                ts_gldas_v21 = pd.Series(index=pd.date_range(start=startdate, end=enddate))

            if (len(np.where(~np.isnan(ts_gldas_v21))[0]) == 0):
                # raise Exception, 'No GLDAS Data available for GPI %i' %gpi
                # print 'No gldas v1 data for gpi %0i' % gpi
                pass
            if data_group.empty:
                data_group['gldas_21'] = ts_gldas_v21
            else:
                data_group = pd.concat([data_group, ts_gldas_v21.rename('gldas_21')], axis=1)

        if 'gldas-merged' in self.H3products.keys():
            try:
                ts_gldas_v2 = self.H3products['gldas_v2'].read_ts(gpi)['086_L1']
                ts_gldas_v21 = self.H3products['gldas_v21'].read_ts(gpi)['SoilMoi0_10cm_inst']
            except:
                ts_gldas_v2 = pd.Series(index=pd.date_range(start=startdate, end=enddate))
                ts_gldas_v21 = pd.Series(index=pd.date_range(start=startdate, end=enddate))

            if (len(np.where(~np.isnan(ts_gldas_v2))[0]) == 0) or (len(np.where(~np.isnan(ts_gldas_v21))[0]) == 0):
                # print 'No gldas v2 or v21 data for gpi %0i' % gpi
                df_gldas_merged = pd.DataFrame(index=pd.date_range(start=startdate, end=enddate, freq='3H'),
                                               data={'gldas-merged': np.nan})
            else:
                df_gldas_merged = pd.concat([ts_gldas_v2.rename('gldas_v2'),
                                             ts_gldas_v21.rename('gldas_v21')], axis=1)

            ts_gldas_merged = df_gldas_merged.mean(axis=1).rename('gldas-merged')

            if data_group.empty:
                data_group['gldas-merged'] = ts_gldas_merged
            else:
                data_group = pd.concat([data_group, ts_gldas_merged.rename('gldas-merged')], axis=1)

        return data_group[startdate:enddate]

if __name__ == '__main__':
    from grid_functions import cells_for_continent
    gpi = 346840
    timeframe = [datetime(1978,10,26), datetime(2016,12,31)]
    data = QDEGdata_D(products=['merra2','CCI_31_COMBINED'])
    ts = data.read_gpi(gpi, timeframe[0], timeframe[1])
    bias_corr_refdata, rxy, pval, ress = regress(
        ts[['CCI_31_COMBINED', 'merra2']].rename(columns={'CCI_31_COMBINED': 'testdata', 'merra2': 'refdata'}))
    ts['merra2_corr'] = bias_corr_refdata
    print ts

