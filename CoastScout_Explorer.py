"""
CoastScout Explorer

Python tools for playing with CoastScout buoy data
Focus on buoy-buoy comparsion with CDIP 201 "Scripps Nearshore" off La Jolla, CA

Contact: Bryan Murray, brmurray@mailbox.org

Change Log:
    + 8/28/2020: First commit

"""


import pandas as pd
import numpy as np
import datetime
import netCDF4
import time
import calendar
import matplotlib.pyplot as plt

def report_timestamp2datetime(report_timestamp):
    import datetime
    fmt = '%Y-%m-%dT%H:%M:%S+00:00'
    return datetime.datetime.strptime(report_timestamp,fmt)

def iso2posix(iso_timestamp):
    '''Convert CoastScout's ISO-formated date strings to datetime, then posix timestamp'''
    import datetime
    dt = datetime.datetime.fromisoformat(iso_timestamp)
    return dt.timestamp()

#%% Import CoastScout Data
# CS "report_time" field is in UTC
    
#p2csdata = r'G:\My Drive\00_CalWave\01_CW_Demo\03_Demo_Test_Site\00_Scripps\00_MetOcean_Data\CoastScout\CoastScout Data\CS2_thruMay18.csv'
p2csdata = r'G:\My Drive\00_CalWave\01_CW_Demo\03_Demo_Test_Site\00_Scripps\00_MetOcean_Data\CoastScout\CoastScout Data\Scripps02_FebAug2020.csv'

cs = pd.read_csv(p2csdata,
                  parse_dates=[2],
                  date_parser=report_timestamp2datetime)

# create "tepoch" column to hold posix time stamps
cs['tepoch'] = cs.report_timestamp.apply(lambda dt: int(dt.timestamp()))

#%% Import CDIP data
import sys
sys.path.insert(1, r'C:\Users\bryan\Documents\GitHub\CDIPtools')
import CDIPtools

# Only download once
#CDIPtools.download_archive('201',r'G:\My Drive\00_CalWave\01_CW_Demo\03_Demo_Test_Site\00_Scripps\00_MetOcean_Data\CDIP201_ScrippsNearshore\Archive_dump')
#CDIPtools.download_realtime('201',r'G:\My Drive\00_CalWave\01_CW_Demo\03_Demo_Test_Site\00_Scripps\00_MetOcean_Data\CDIP201_ScrippsNearshore\Archive_dump')

stn = '201'
dstfolder = r'G:\My Drive\00_CalWave\01_CW_Demo\03_Demo_Test_Site\00_Scripps\00_MetOcean_Data\CDIP201_ScrippsNearshore\Archive_dump\\'
dstfile = r'tempnc' + r'.nc'
dst = dstfolder + r'\\' + dstfile

# Retreive the file
import urllib.request
#url = r'http://thredds.cdip.ucsd.edu/thredds/fileServer/cdip/archive/' + stn + r'/' + stn + 'p1_historic.nc'
url = r'http://thredds.cdip.ucsd.edu/thredds/fileServer/cdip/realtime/' + stn + 'p1_rt.nc'
#urllib.request.urlretrieve(url, dst)


#%%
# Open local copy of netCDF archive file
p2archivefile = r'G:\My Drive\00_CalWave\01_CW_Demo\03_Demo_Test_Site\00_Scripps\00_MetOcean_Data\CDIP201_ScrippsNearshore\Archive_dump\CDIP201_feb-may2020.nc'
ds = netCDF4.Dataset(p2archivefile)
#p2CDIP201_rt = r'G:\My Drive\00_CalWave\01_CW_Demo\03_Demo_Test_Site\00_Scripps\00_MetOcean_Data\CDIP201_ScrippsNearshore\Archive_dump\CDIP201_feb-may2020.nc'
#ds = netCDF4.Dataset(p2CDIP201_rt)
ds.set_always_mask(False)

# Create pandas DataFrame for these guys
dsTime = ds.variables['waveTime'][:]
timeall = [datetime.datetime.fromtimestamp(t) for t in dsTime] # Convert ncTime variable to datetime stamps
Hs = ds.variables['waveHs'][:]
Tp = ds.variables['waveTp'][:]
Dp = ds.variables['waveDp'][:] 

data = {'t':dsTime,'dt':timeall,'Hs':Hs,'Tp':Tp,'Dp':Dp}
ns = pd.DataFrame(data)
ns = ns.set_index('t')


#%% figure out how to align CoastScout and CDIP201 times
# round coastscout down to half-hour
#cs.['t'] = cs['t'].apply(lambda dt: datetime.datetime(dt.year,dt.month,dt.day,dt.hour,30*round(((float(dt.minute)) + float(dt.second)/60) / 30)))
#cs.index = cs.index.apply(lambda ti: ti - (ti%1800)) # can't APPLY to an index
cs['thh'] = cs.tepoch
cs.thh= cs.thh.apply(lambda ti: ti - (ti%1800)) 
#cs['dt'] = cs.thh.apply(unix2datetime)
cs = cs.set_index('thh')


# round buoy 201 'nearshore' down to half-hour
#ns['t'] = ns['t'].apply(lambda dt: datetime.datetime(dt.year, dt.month, dt.day, dt.hour,30*round(((float(dt.minute)) + float(dt.second)/60) / 30)))
ns['thh'] = ns.index#-8*60*60
ns.thh= ns.thh.apply(lambda ti: ti - (ti%1800))
#ns['dt'] = ns.thh.apply(unix2datetime)
ns = ns.set_index('thh')



#%%
df = pd.merge(ns,cs,left_index=False,right_index=True,how='inner',on=['thh','thh']) #<-- returns empty
df.dropna(inplace=True)
df['dt'] = df.tepoch.apply(lambda thh: datetime.datetime.fromtimestamp(thh))
#common = ns.index.intersection(cs.index)
#cs.thh[cs.thh.isin(ns.thh)]

# Identify SWS2 times
sws2 = df[(df.Hs>1.3) & (df.Hs<1.4) & (df.Tp>7.5) & (df.Tp<8.5)]

#%% Plot Tp time series
fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(df.index,df.Tp)
ax1.plot(df.index,df.wave_peak_period_seconds)
ax1.plot(df.index-8*60*60,df.wave_peak_period_seconds)
ax1.plot(df.index+8*60*60,df.wave_peak_period_seconds)
ax1.legend([r'CDIP201',r'CoastScout2',r'CS2 - 8hr',r'CS2 + 8hr'])

#%% Plot Hs time series
fig=plt.figure()
pHs=fig.add_subplot(111)
pHs.plot(df.dt,df.Hs,'b-',zorder=2)
pHs.plot(df.dt,df.wave_significant_height_m,'g-',zorder=3)
pHs.legend([r'CDIP201',r'CoastScout2'])
pHs.set_ylabel('Hs, m',fontsize=18)
pHs.set_ylim(0,3)
pHs.set_yticks([0,1.,2.,3.])
plt.ylabel('Significant Wave Height Hs (m)')
pHs.grid(b=True, which='major', color='dimgrey', linestyle='--')

# Set title
plt.suptitle('CDIP201 vs CoastScout, with times of interest',fontsize=22)

# March 18 Storm
plt.axvline(x=datetime.datetime(2020,3,18,18,0,0),c='r',zorder=1)
plt.axvline(x=datetime.datetime(2020,3,18,22,0,0),c='r',zorder=1)

# May 16 SWS2
plt.axvline(x=datetime.datetime(2020,5,16,14,0,0),c='m',zorder=1)
plt.axvline(x=datetime.datetime(2020,5,16,17,0,0),c='m',zorder=1)


#%% Time Series Compendium
f, (pHs,pTp,pDp) = plt.subplots(3,1,sharex=True,figsize=(15,10),num='CDIP201 vs CoastScout, Feb-Aug 2020')

# Set title
plt.suptitle('CDIP201 vs CoastScout, Feb-Aug 2020',fontsize=22)

# Create 3 stacked subplots for three parameters: (Hs,Tp,Dp)
pHs.plot(df.dt,df.Hs,'b-')
pHs.plot(df.dt,df.wave_significant_height_m,'g-')
pHs.legend(['CDIP201','Scripps02'],fontsize=14,loc='upper right')
pHs.set_ylabel('Hs, m',fontsize=18)
pHs.set_ylim(0,3)
pHs.set_yticks([0,1.,2.,3.])
pHs.grid(b=True, which='major', color='dimgrey', linestyle='--')

pTp.plot(df.dt,df.Tp,'b-')
pTp.plot(df.dt,df.wave_peak_period_seconds,'g-')
pTp.legend(['CDIP201','Scripps02'],fontsize=14,loc='upper right')
pTp.set_ylabel('Tp, sec',fontsize=18)
pTp.set_ylim(2,22)
pTp.set_yticks([4,8,12,16,20])
pTp.grid(b=True, which='major', color='dimgrey', linestyle='--')

pDp.plot(df.dt,df.Dp,'b-')
pDp.plot(df.dt,df.wave_mean_direction_degrees,'g-')
pDp.legend(['CDIP201','Scripps02'],fontsize=14,loc='upper right')
pDp.set_ylabel('Dp, degrees', fontsize=18) 
pDp.set_ylim(180,360)
pDp.set_yticks([180,225,270,315,360])
pDp.grid(b=True, which='major', color='dimgrey', linestyle='--')


#%% Comparison Compendium
fig = plt.figure()

import numpy.polynomial.polynomial as poly
# poly.polyfit returns [A,B,C] for A + Bx + Cx^2 + Dx^3...

# Hs
coefs_Hs = poly.polyfit(df.Hs,df.wave_significant_height_m, 1)
ffit_Hs = poly.polyval(df.Hs, coefs_Hs)

ax1 = fig.add_subplot(221)
ax1.plot(df.Hs,df.wave_significant_height_m,'b*')
ax1.plot([0,2.25],[0,2.25],'m-')
ax1.plot(df.Hs, ffit_Hs,'r-')
ax1.set_xlim([0,2.25])
ax1.set_ylim([0,2.25])
ax1.set_xlabel('CDIP 201 Nearshore')
ax1.set_ylabel('CoastScout 2')
ax1.legend(['Data','Linear','y='+str(round(coefs_Hs[1],2)) + 'x + ' + str(round(coefs_Hs[0],3))])
ax1.text(.85,.25,'Negative slope implies\n that CDIP over-estimates Hs,\n or CoastScout under-estimates')

# Tp
coefs_Tp = poly.polyfit(df.Tp,df.wave_peak_period_seconds, 1)
ffit_Tp = poly.polyval(df.Tp, coefs_Tp)

ax2 = fig.add_subplot(222)
ax2.plot(df.Tp,df.wave_peak_period_seconds,'b*')
ax2.plot([0,22],[0,22],'m-')
ax2.plot(df.Tp, ffit_Tp,'r-')
ax2.set_xlim([0,22])
ax2.set_ylim([0,22])
ax2.set_xlabel('CDIP 201 Nearshore')
ax2.set_ylabel('CoastScout 2')
ax2.legend(['Data','Linear','y='+str(round(coefs_Tp[1],2)) + 'x + ' + str(round(coefs_Tp[0],3))])
ax2.text(10,18,'Poor agreement due to Tp binning. \n Must dig further into swell resolution codes.')

# Dp
coefs_Dp = poly.polyfit(df.Dp,df.wave_mean_direction_degrees, 1)
ffit_Dp = poly.polyval(df.Dp, coefs_Dp)

ax3 = fig.add_subplot(223)
ax3.plot(df.Dp,df.wave_mean_direction_degrees,'b*')
ax3.plot([200,360],[200,360],'m-')
ax3.plot(df.Dp, ffit_Dp,'r-')
ax3.set_xlim([200,360])
ax3.set_ylim([200,360])
ax3.set_xlabel('CDIP 201 Nearshore')
ax3.set_ylabel('CoastScout 2')
ax3.legend(['Data','Linear','y='+str(round(coefs_Dp[1],2)) + 'x + ' + str(round(coefs_Dp[0],3))])
ax3.text(220,300,'Poor agreement due to Tp binning. \n Must dig further into swell resolution codes.')

#%% Just Compare Hs
fig = plt.figure()

import numpy.polynomial.polynomial as poly
# poly.polyfit returns [A,B,C] for A + Bx + Cx^2 + Dx^3...

# Hs
coefs_Hs = poly.polyfit(df.Hs,df.wave_significant_height_m, 1)
ffit_Hs = poly.polyval(df.Hs, coefs_Hs)

ax1 = fig.add_subplot(111)
ax1.plot(df.Hs,df.wave_significant_height_m,'b*')
ax1.plot([0,2.25],[0,2.25],'m-')
ax1.plot(df.Hs, ffit_Hs,'r-')
ax1.set_xlim([0,2.25])
ax1.set_ylim([0,2.25])
ax1.set_xlabel('Hs from CDIP 201 Nearshore',fontsize=14,weight='bold')
ax1.set_ylabel('Hs from CoastScout',fontsize=14,weight='bold')
ax1.set_title('CDIP201 vs CoastScout 2, Feb-May 2020',fontsize=18,weight='bold')
ax1.legend(['Data','Ideal','y='+str(round(coefs_Hs[1],2)) + 'x + ' + str(round(coefs_Hs[0],3))])
ax1.text(.85,.25,'Negative slope implies\n that CDIP over-estimates Hs,\n or CoastScout under-estimates',fontsize=14,weight='bold',color='r')


#%% JPD for CoastScout 2

import sys
sys.path.insert(1, r'C:\Users\bryan\Documents\GitHub\WIStools')
import WIStools
import mhkit
coastscout2 = WIStools.WISstation('00002')

tp_edges = np.arange(2.,16.5,0.5)
tp_centers = tp_edges[0:-1]+0.25
hs_edges = np.arange(0.0,3.25,0.25)
hs_centers = hs_edges[0:-1]+0.125

# Create some needed columns
cs['J'] = cs.apply(lambda x:WIStools.SpecWavePower(x.wave_significant_height_m,x.wave_peak_period_seconds),axis=1)
cs['Te'] = cs.apply(lambda y:WIStools.Tp2Te(y.wave_peak_period_seconds),axis=1)
cs.dropna(inplace=True)
coastscout2.Occurrence = mhkit.wave.performance.wave_energy_flux_matrix(cs.wave_significant_height_m,
                                                    cs.wave_peak_period_seconds,
                                                    cs.J,
                                                    'frequency',
                                                    hs_centers,
                                                    tp_centers)

WIStools.plot_matrix_edges(coastscout2.Occurrence,tp_edges,hs_edges,xlabel='Tp',ylabel='Hm0',zlabel='Occurrence',show_values=False)

#%% JPD for CDIP201 (during that time period)

import sys
sys.path.insert(1, r'C:\Users\bryan\Documents\GitHub\WIStools')
import WIStools
import mhkit
cdip201 = WIStools.WISstation('000201')

tp_edges = np.arange(2.,16.5,0.5)
tp_centers = tp_edges[0:-1]+0.25
hs_edges = np.arange(0.0,3.25,0.25)
hs_centers = hs_edges[0:-1]+0.125

# Create some needed columns
ns['J'] = ns.apply(lambda x:WIStools.SpecWavePower(x.Hs,x.Tp),axis=1)
ns['Te'] = ns.apply(lambda y:WIStools.Tp2Te(y.Tp),axis=1)
ns.dropna(inplace=True)
cdip201.Occurrence = mhkit.wave.performance.wave_energy_flux_matrix(ns.Hs,
                                                    ns.Tp,
                                                    ns.J,
                                                    'frequency',
                                                    hs_centers,
                                                    tp_centers)

WIStools.plot_matrix_edges(cdip201.Occurrence,tp_edges,hs_edges,xlabel='Tp',ylabel='Hm0',zlabel='Occurrence',show_values=False)

#%% FREQUENCY BASED METRICS

#%%

period_bins = np.flipud(1/ds.variables['waveFrequency'][:]) #wave period bins, irregularly spaced over [1.724:40]

plt.plot(period_bins,np.flipud(ds.variables['waveEnergyDensity'][0,:]))





















