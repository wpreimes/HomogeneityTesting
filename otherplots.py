# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 17:06:28 2017

@author: wpreimes
"""
import numpy as np
import math
import os,glob
from datetime import datetime
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

import re
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def spatial_plot_quarter_grid(data, tags, \
                                llcrnrlat=-90.,\
                                urcrnrlat=90.,\
                                llcrnrlon=-180.,
                                urcrnrlon=180.,\
                                continent=None,\
                                figsize=(20,10),\
                                title='',\
                                cblabel='',\
                                textbox=None,\
                                cbrange=(0,1),\
                                cmap='jet',\
                                labelpad=-70,\
                                fontsize=20,\
                                tight=False,\
                                path=r"H:\Maps\Ascat"):
    """
    die routine plottet aus einem data frame, der als index die gpis hat, 
    jene spalten, die als "tags" übergeben werden
    """
    if continent == 'NA':
        llcrnrlat=25.
        urcrnrlat=50.
        llcrnrlon=-130.
        urcrnrlon=-60.
    if continent == 'AUS':
        llcrnrlat=-40.
        urcrnrlat=-10.
        llcrnrlon=110.
        urcrnrlon=155.
    
    # quarter degree meshgrid generieren
    lons = (np.arange(360*4)*0.25)-179.875
    lats = (np.arange(180*4)*0.25)-89.875
    lons,lats = np.meshgrid(lons,lats)
    
    
    
    for tag in tags:
        
        img = np.empty(lons.size,dtype='float32')
        img.fill(None)
        # img[gpi] = data[spalte]
        
        #img[data['gpi_quarter'].values] = data[tag] 
        img[data.index.values] = data[tag].values
        
        # mask array where invalid values (nans) occur
        img_masked = np.ma.masked_invalid(img.reshape((180*4,360*4)))

        f = plt.figure(num=None, figsize=figsize, dpi=90, facecolor='w', edgecolor='k')
    
        m = Basemap(projection='mill',llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,\
                llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon,resolution='c',)

        m.drawcoastlines()
        m.drawcountries()
        
        im = m.pcolormesh(lons,lats, img_masked, cmap=cmap, latlon=True)
        
        
        im.set_clim(vmin=cbrange[0],vmax=cbrange[1])
        
        cb = m.colorbar(im,"bottom", size="5%", pad="8%")

        for t in cb.ax.get_xticklabels():
             t.set_fontsize(fontsize)
             
        if cblabel is None:
            cb.set_label(tag,fontsize=fontsize, labelpad=labelpad)
        else:
            cb.set_label(cblabel,fontsize=fontsize, labelpad=labelpad)
        
        if title == 'tag':
            tmp_title = tag
        else:
            tmp_title = title
            
        plt.title(tmp_title,fontsize=fontsize)
        
        if tight is True:
            plt.tight_layout()
        
        if textbox:
            plt.annotate(textbox, xy=(0.01, 0.9), xycoords='axes fraction')
        
        
        if path is None:
            plt.show()
        else:
            if not os.path.exists(path):
                os.makedirs(path)
            plt.savefig(path + '\\'+ title + '.png', dpi = f.dpi)
            plt.close()
            
 
def compare_RTM_RTG(workdir,model_prod):
    
    times={'ISMN-merge':{'starttimes':[],'endtimes':[],'breaktimes':[]},
                model_prod:{'starttimes':[],'endtimes':[],'breaktimes':[]}}
    files={'ISMN-merge':[],model_prod:[]}            
    for ref_data in ['ISMN-merge',model_prod]:
        fileslist=glob.glob(os.path.join(workdir,"DF_Points_%s*.csv" %ref_data))
        files[ref_data]=fileslist
        for filename in fileslist:
            filename=filename.replace(workdir,'')
            filename=filename.replace('.','_')
            splitname=filename.split('_')
            times[ref_data]['starttimes'].append(splitname[3])
            times[ref_data]['endtimes'].append(splitname[5])
            times[ref_data]['breaktimes'].append(splitname[4])
    ''' 
    if (times['ISMN_merge']['starttimes']==times[model_prod]['starttimes']) and \
        (times['ISMN_merge']['endtimes']==times[model_prod]['endtimes']) and \
        (times['ISMN_merge']['breaktimes']==times[model_prod]['breaktimes']):
    '''

    rows=int(math.ceil(float(len(files['ISMN-merge']))/2.))
    fig,axs=plt.subplots(rows,2,figsize=(15, 6), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace=1,wspace=1)
    axs=axs.ravel()
        
    for i,(ismn_file, model_file) in enumerate(zip(files['ISMN-merge'],files[model_prod])):
        if times['ISMN-merge']['breaktimes'][i]==times[model_prod]['breaktimes'][i]:
            DF_Points_ismn=pd.read_csv(os.path.join(ismn_file),index_col=0)
            DF_Points_model=pd.read_csv(os.path.join(model_file),index_col=0)
            
            if (DF_Points_ismn.index==DF_Points_model.index).all():
                DF_Points_merged=DF_Points_model[['lat','lon','h_all']].rename(columns={'h_all':model_prod})
                DF_Points_merged['ISMN-merge']=DF_Points_ismn['h_all']               
                DF_Points_merged=DF_Points_merged.dropna(how='any')
                DF_Points_merged['diff']=DF_Points_merged[model_prod]-DF_Points_merged['ISMN-merge']
                
                DF_Points_merged.hist(column='diff',bins=[-3,-2,-1,0,1,2,3],ax=axs[i])
                #axs[i].set_xlim((-3,3))
                axs[i].set_title(str(times[model_prod]['breaktimes'][i]))
                
    fig.savefig(os.path.join(workdir,'RTM_vs_RTG'))


def show_tested_gpis(workdir,ref_prod):
    
    '''
    Calculate spatial plots for the areas where Homogeneity Tests were (not)
    performed
    '''
    
    fileslist=glob.glob(os.path.join(workdir,"DF_Points_%s*.csv"%ref_prod))
    starttimes=[]
    endtimes=[]
    breaktimes=[]

    for filename in fileslist:
        filename=filename.replace(workdir,'')
        filename=filename.replace('.','_')
        splitname=filename.split('_')
        starttimes.append(datetime.strptime(splitname[3], '%Y-%m-%d'))
        endtimes.append(datetime.strptime(splitname[5], '%Y-%m-%d'))
        breaktimes.append(datetime.strptime(splitname[4], '%Y-%m-%d'))
       
    for f,breaktime,starttime,endtime in zip(fileslist,breaktimes,starttimes,endtimes):
        DF_Points=pd.read_csv(f,index_col=0,usecols=[0,6,7],header=0,
                                       names=['gpi_quarter','h_all','message'])
        
        DF_Points['TestGroups']=np.nan
        
        for gpi in DF_Points.index.values:
            h=pd.isnull(DF_Points.h_all.loc[gpi])
            message=DF_Points.message.loc[gpi]
            if h != True and message == 'Processing OK': 
                DF_Points.set_value(gpi,'TestGroups',0)
            elif h == True and 'Dataseries Length' in message:
                DF_Points.set_value(gpi,'TestGroups',1)
            elif h == True and 'Spearman' in message:
                DF_Points.set_value(gpi,'TestGroups',2)
            else: 
                DF_Points.set_value(gpi,'TestGroups',3)
    
        colors = ['green','pink','red','white'] 
        cmap=LinearSegmentedColormap.from_list('TestGroups',
                                               colors, 
                                               N=4)
                                               
        spatial_plot_quarter_grid(DF_Points,
                                  tags=['TestGroups'],
                                  title='RTM_coverage_%s_(breaktime-%s,timeframe-%s-%s)'%(ref_prod,breaktime.strftime("%Y-%m-%d"),
                                                                                          starttime.strftime("%Y-%m-%d"),
                                                                                          endtime.strftime("%Y-%m-%d")),
                                  cbrange=(0,3),
                                  cmap=cmap,
                                  cblabel='0=Tested,2=Data length, 3=correlation, 4=other',
                                  #continent='NA',
                                  path=workdir)


'''
[ops[0],4.,ops[1],4.,ops[2],4.,ops[3],4.,ops[4]]

   for values,ops in zip(periods,calcs):
        
        ops=[ops[0],4.,ops[1],4.,ops[2],4.,ops[3],4.,ops[4]]
        r=[]

        while ops:
            if ops[0]!=4.0:
                r.append(values[0])
                del values[0]
                del ops[0]
            elif ops[0] == 4.0:
                values[1]+=values[0]
                del ops[0]
                del values[0]
        r.append(values[-1])
        results.append(r)      
'''


def inhomo_plot_with_stats(workdir):
        fileslist=glob.glob(os.path.join(workdir,"DF_Points*.csv"))
        starttimes=[]
        breaktimes=[]
        endtimes=[]
        
        for filename in fileslist:
            DF_Points=pd.read_csv(filename,index_col=0)
            
            filename=filename.replace(workdir,'')
            filename=filename.replace('.','_')
            splitname=filename.split('_')
            ref_prod=splitname[2]
            starttime=datetime.strptime(splitname[3], '%Y-%m-%d')
            endtime=datetime.strptime(splitname[5], '%Y-%m-%d')
            breaktime=datetime.strptime(splitname[4], '%Y-%m-%d')
        
            
            all_gpis=DF_Points.h_all.size
            tested_gpis=DF_Points['h_all'].loc[DF_Points['h_all'].isin([1.,2.,3.,4.])].size
            hwk_gpis=DF_Points['h_all'].loc[DF_Points['h_all'].isin([1.,3.])].size
            hfk_gpis=DF_Points['h_all'].loc[DF_Points['h_all'].isin([2.,3.])].size
            hboth_gpis=DF_Points['h_all'].loc[DF_Points['h_all'].isin([3.])].size
            
            
            
            all_tested=(float(tested_gpis)/float(all_gpis))*100.
            try:
                wk_tested=(float(hwk_gpis)/float(tested_gpis))*100.         
                fk_tested=(float(hfk_gpis)/float(tested_gpis))*100.
                both_tested=(float(hboth_gpis)/float(tested_gpis))*100.
            except:
                wk_tested,fk_tested,both_tested=0,0,0

                
            colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1),(.5,.5,.5)] 
            cmap=LinearSegmentedColormap.from_list('HomogeneityGroups',
                                                   colors, 
                                                   N=4)
                                                   
            spatial_plot_quarter_grid(DF_Points,
                                      tags=['h_all'],
                                      title='HomogeneityTest_%s_(breaktime-%s)'%(ref_prod,breaktime.strftime("%Y-%m-%d")),
                                      cbrange=(1,4),
                                      cmap=cmap,
                                      textbox='%% Pixel tested:%.2f\n'%all_tested+\
                                              '%% WK of tested:%.2f\n'%wk_tested +\
                                              '%% FK of tested:%.2f\n'%fk_tested +\
                                              '%% Both of tested:%.2f'%both_tested,
                                      cblabel='1=WK,2=FK,3=both,4=None',
                                      #continent='NA',
                                      path=workdir)


  
def calc_longest_homogeneous_period(workdir,test_prod,ref_prod):
    #TODO: This function is bad...
    fileslist=glob.glob(os.path.join(workdir,"DF_Points_%s*.csv"%ref_prod))
    starttimes=[]
    endtimes=[]
    breaktimes=[]
    periodsizes=[]
    for filename in fileslist:
        filename=filename.replace(workdir,'')
        filename=filename.replace('.','_')
        splitname=filename.split('_')
        starttimes.append(datetime.strptime(splitname[3], '%Y-%m-%d'))
        endtimes.append(datetime.strptime(splitname[5], '%Y-%m-%d'))
        breaktimes.append(datetime.strptime(splitname[4], '%Y-%m-%d'))
    
    DF_Period = pd.concat([pd.read_csv(f,index_col=0,usecols=[0,6],header=0,
                                       names=['gpi_quarter',breaktimes[i]]) for i,f in enumerate(fileslist)],axis=1)   
    
    i=1
    for starttime,breaktime,endtime in zip(starttimes,breaktimes,endtimes):
        #DF_Period['%i.1'%i] = np.where(DF_Period[breaktime]==4.,(endtime-starttime).days, ((breaktime-starttime).days),(endtime-breaktime).days))
        DF_Period['%i'%i] = (breaktime-starttime).days
        periodsizes.append((breaktime-starttime).days)
        #DF_Period['%i.2'%i] = (endtime-breaktime).days
        i+=1
    DF_Period['%i'%i]=(endtimes[-1]-breaktimes[-1]).days
    periods=map(list,DF_Period[DF_Period.columns.values[-i:]].values)
    calcs=map(list,DF_Period[breaktimes].values)
    results=[]    
    for values,ops in zip(periods,calcs):
        ops.append([4.0,np.nan])
        r=[]

        while ops:
            if ops[0]!=4.0:
                r.append(values[0])
                del values[0]
                del ops[0]
            elif ops[0] == 4.0:
                values[1]+=values[0]
                del ops[0]
                del values[0]
        #r.append(values[-1])
        results.append(r)        
    
    max_period=[]
    for periods in results:
        if max(periods) in periodsizes:
            max_period.append(np.nan)
        else:
            max_period.append(max(periods)/365.)

    DF_Period['max_Period']=max_period



    #colors = [(1, 0, 0),(1, 0, 0),(0.5, 1, 1),(0.5, 1, 1),(0,0,1),(0,0,1)] 

                   
    #cmap=LinearSegmentedColormap.from_list('HomogeneityGroups',colors,N=10)

    cmap=cmap = plt.cm.jet

    cmaplist=[(0.5,0,0,1),(0.7,0,0,1),(1,0,0,1),(1,0.2,0.2,1),(1,.5,.5,1),
              #(0.3,0.3,0.3,1),(0.5,0.5,0.5,1),(0.5,0.5,0.5,1),
              (0,0.7,0.7,1),(0,0.7,0.7,1),(0,0.7,0.7,1),(0,0.7,0.7,1),
              (0,0.6,0.7,1),(0,0.3,1,1),(0,0,0.8,1),(0,0,0.5,1)]

    if test_prod=='cci_22':
        cbrange=(8,24)
    elif test_prod=='cci_31':
        cbrange=(6,20)
              
    cmap = cmap.from_list('homogmap', cmaplist, cmap.N)
    #cmaplist = [cmap(i) for i in range(cmap.N)]
    
    spatial_plot_quarter_grid(DF_Period,
                              tags=['max_Period'],
                              title='RTM_maxHomogPeriod_%s starttime-%s endtime-%s'%(ref_prod,
                                                                                     starttimes[0].strftime("%Y %m, %d"),
                                                                                     endtimes[-1].strftime("%Y %m, %d")),
                              cbrange=cbrange,
                              cmap=cmap,
                              cblabel='Length of longest homogeneous period [years]',
                              #continent='NA',
                              path=workdir)                                 
    

#inhomo_plot_with_stats(r'H:\workspace\HomogeneityTesting\output\plottest')
#calc_longest_homogeneous_period(r'H:\workspace\HomogeneityTesting\output\v11','cci_22','merra2')
#show_tested_gpis(r'H:\workspace\HomogeneityTesting\output\v11','merra2')
#inhomo_plot_with_stats(r'H:\workspace\HomogeneityTesting\output\v12')
#compare_RTM_RTG(r'H:\workspace\HomogeneityTesting\output\v25','merra2')
'''
show_tested_gpis(r'H:\workspace\HomogeneityTesting\output\plottest','merra2')
'''