"""
CoastScout_Explorer_EWTEC

Python tools for exploring CoastScout/CDIP data
The notebook creates the relevant outputs for the EWTEC paper


Contact: Bryan Murray, bryan@calwave.energy

Change Log:
    + 12/29/20 [BRM]: First commit. New notebook created from elements of "CoastScout_Explorer" and "CoastScout_Explorer_Frequency"
    + 

"""

#%%
import pandas as pd
import numpy as np
import datetime
import netCDF4
import time
import calendar
import matplotlib.pyplot as plt

#%% Import and Combine Data Sets

#------------------------------------------------------
# Import CoastScout Data
# CS "report_time" field is in UTC

def report_timestamp2datetime(report_timestamp):
    import datetime
    fmt = '%Y-%m-%dT%H:%M:%S+00:00'
    return datetime.datetime.strptime(report_timestamp,fmt)

p2csdata = r'G:\My Drive\00_CalWave\01_CW_Demo\07_Partners Exchange\18_MarineLabs\EWTEC\Scripps02_FebAug2020.csv'
cs = pd.read_csv(p2csdata,
                  parse_dates=[2],
                  date_parser=report_timestamp2datetime)

# create "tepoch" column to hold posix time stamps
cs['tepoch'] = cs.report_timestamp.apply(lambda dt: int(dt.timestamp()))

#------------------------------------------------------
# Open local copy of CDIP201 netCDF "real-time" file
# it happens to cover February - August 2020
p2CDIP201_rt = r'G:\My Drive\00_CalWave\01_CW_Demo\07_Partners Exchange\18_MarineLabs\EWTEC\201p1_rt_FebAug2020.nc'
ds = netCDF4.Dataset(p2CDIP201_rt)
ds.set_always_mask(False)

# Create pandas DataFrame for CDIP201
dsTime = ds.variables['waveTime'][:]
timeall = [datetime.datetime.fromtimestamp(t) for t in dsTime] # Convert ncTime variable to datetime stamps
Hs = ds.variables['waveHs'][:]
Tp = ds.variables['waveTp'][:]
Dp = ds.variables['waveDp'][:] 

data = {'t':dsTime,'dt':timeall,'Hs':Hs,'Tp':Tp,'Dp':Dp}
cdip201 = pd.DataFrame(data)
cdip201 = cdip201.set_index('t')

# To align Scripps02 and CDIP201, we'll create a column holding posix time stamps
# we can then use the timestamps for a join

# First the CoastScout buoy
cs['thh'] = cs.tepoch
cs.thh= cs.thh.apply(lambda ti: ti - (ti%1800)) 
cs = cs.set_index('thh')

# round buoy 201 'nearshore' down to half-hour
cdip201['thh'] = cdip201.index
cdip201.thh= cdip201.thh.apply(lambda ti: ti - (ti%1800))
cdip201 = cdip201.set_index('thh')

df = pd.merge(cdip201,cs,left_index=False,right_index=True,how='inner',on=['thh','thh']) #<-- returns empty
df.dropna(inplace=True)
df['dt'] = df.tepoch.apply(lambda thh: datetime.datetime.fromtimestamp(thh))
#df.to_csv(r'path\cdip_cs.csv')

#%% Comparison Compendium
#
#
fig = plt.figure(1)
plt.suptitle('CDIP201 vs CoastScout 2, Feb-May 2020',fontsize=18,weight='bold')

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
ax1.set_title('Hs',fontsize=18,weight='bold')

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
ax2.text(10,18,'Poor agreement due to Tp binning. \n Must dig further into swell resolution codes.',color='r')
ax2.set_title('Tp',fontsize=18,weight='bold')

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
ax3.text(220,300,'Poor agreement due to Tp binning. \n Must dig further into swell resolution codes.',color='r')
ax3.set_title('Dp',fontsize=18,weight='bold')

#%% Storm Time Series
#
#

fig=plt.figure(2)
plt.suptitle('Storm March 18',fontsize=18,weight='bold')

#  Top plot: march 18 storm +/- 1 week
pHs=fig.add_subplot(211)
start_date = datetime.datetime(2020,3,11,18,0,0)
end_date = datetime.datetime(2020,3,27,18,0,0)
mask = (df.dt>start_date) & (df.dt<end_date)

pHs.plot(df.dt.loc[mask],df.Hs.loc[mask],'b-',zorder=2)
pHs.plot(df.dt.loc[mask],df.wave_significant_height_m.loc[mask],'g-',zorder=3)
pHs.legend([r'CDIP201',r'CoastScout2'])
pHs.set_ylabel('Hs, m',fontsize=18)
pHs.set_ylim(0,3)
pHs.set_yticks([0,1.,2.,3.])
pHs.grid(b=True, which='major', color='dimgrey', linestyle='--')
plt.axvline(x=datetime.datetime(2020,3,18,18,0,0),c='r',zorder=1)
plt.axvline(x=datetime.datetime(2020,3,20,18,0,0),c='r',zorder=1)

#bottom plot: zoom in on March 18 storm
pMar18=fig.add_subplot(212)

start_date = datetime.datetime(2020,3,18,18,0,0)
end_date = datetime.datetime(2020,3,20,18,0,0)
mask = (df.dt>start_date) & (df.dt<end_date)

pMar18.plot(df.dt.loc[mask],df.Hs.loc[mask],'b-',zorder=2)
pMar18.plot(df.dt.loc[mask],df.wave_significant_height_m.loc[mask],'g-',zorder=3)
pMar18.legend([r'CDIP201',r'CoastScout2'])
pMar18.set_ylabel('Hs, m',fontsize=18)
pMar18.set_ylim(0,3)
pMar18.set_yticks([0,1.,2.,3.])
pMar18.grid(b=True, which='major', color='dimgrey', linestyle='--')


#%% Frequency Based Analysis
#
#
#------------------------------------------------------
# North-West-Up (NWU)
#   Open NWU CoastScout "first 5" files and combine into two time ranges
#   20200318T181000Z UTC - 20200319T064000Z UTC
#   20200807T181000Z - 20200808T051000Z
    
NWU_f5_data_folder = r'G:\My Drive\00_CalWave\01_CW_Demo\07_Partners Exchange\18_MarineLabs\EWTEC\NWU Coordinates'

import os
NWU_f5_folders = os.listdir(NWU_f5_data_folder)
NWU_march_f5_folders = [os.path.join(NWU_f5_data_folder,folder,r'fourier_coefficients.csv') for i,folder in enumerate(NWU_f5_folders) if folder[4:6]=='03']
NWU_aug_f5_folders = [os.path.join(NWU_f5_data_folder,folder,r'fourier_coefficients.csv') for i,folder in enumerate(NWU_f5_folders) if folder[4:6]=='08']
    
NWU_marf5 = pd.concat(pd.read_csv(f,skiprows=9,header=None,index_col=False,names=['Hz','mms','a1','b1','a2','b2']) for f in NWU_march_f5_folders[23:24])
NWU_augf5 = pd.concat(pd.read_csv(f,skiprows=9,header=None,index_col=False,names=['Hz','mms','a1','b1','a2','b2']) for f in NWU_aug_f5_folders[17:18])


#------------------------------------------------------
# Spectral axis variables (of size (rows,64) or (64,))
Ed = ds.variables['waveEnergyDensity'][:] #band energy density (m*m/hz)
Fq = ds.variables['waveFrequency'][:] #band center frequency
bw = ds.variables['waveBandwidth'][:]
A1 = ds.variables['waveA1Value'][:]
B1 = ds.variables['waveB1Value'][:]
A2 = ds.variables['waveA2Value'][:]
B2 = ds.variables['waveB2Value'][:]

# Derived Series
Eb = Ed*bw # energy in contained in bands
m0 = np.sum(Eb,1) # 0th moment, sum of all energy in the spectrum for each time step. Premier among the "First 5"
H13 = 4*np.sqrt(m0) # Hs = H13


#%% March First 5
# March 18
# 1584556200 = March 18, 2020 18:30:00 = CDIP201 record 903
# 20200319T064000Z = 1584601200 = CDIP201 record 928
icdip = 928

#------------------------------------------------------
#   Plot First 5 as a function of Frequency

f, (pEd,pA1,pB1,pA2,pB2) = plt.subplots(5,1,sharex=True,figsize=(10,6),num='CDIP201 vs CoastScout, "First 5" Components')
# Set title
plt.suptitle('CDIP201 vs CoastScout "First 5" Components \n March 19, 2020 0700 GMT',fontsize=12)

#pEd.plot(ENU_marf5.Hz,ENU_marf5.mms,'g-')
#pEd.plot(NED_marf5.Hz,NED_marf5.mms,'r-')
pEd.plot(NWU_marf5.Hz,NWU_marf5.mms,'m-')
pEd.plot(Fq,Ed[icdip,:],'b-')
pEd.set_ylabel('m*m/Hz')
pEd.set_title('Energy Density')
pEd.grid(b=True, which='major', color='dimgrey', linestyle='--')
#pEd.legend(['CoastScout ENU','CoastScout NED','CDIP201'])
pEd.legend(['CoastScout NWU','CDIP201'])
pEd.set_xlim([0,1])

#pA1.plot(ENU_marf5.Hz,ENU_marf5.a1,'g-')
#pA1.plot(NED_marf5.Hz,NED_marf5.a1,'r-')
pA1.plot(NWU_marf5.Hz,NWU_marf5.a1,'m-')
pA1.plot(Fq,A1[icdip,:],'b-')
#pA1.text(0.8,0.5,'A1',fontsize=10)
pA1.set_ylabel('A1')
pA1.set_ylim([-1,1])
pA1.set_yticks([-1,0,1])
pA1.set_xlim([0,1])
pA1.grid(b=True, which='major', color='dimgrey', linestyle='--')

#pB1.plot(ENU_marf5.Hz,ENU_marf5.b1,'g-')
#pB1.plot(NED_marf5.Hz,-NED_marf5.b1,'r-')
pB1.plot(NWU_marf5.Hz,NWU_marf5.b1,'m-')
pB1.plot(Fq,B1[icdip,:],'b-')
#pB1.text(0.8,0.5,'B1',fontsize=12)
pB1.set_ylabel('B1')
pB1.set_ylim([-1,1])
pB1.set_yticks([-1,0,1])
pB1.set_xlim([0,1])
pB1.grid(b=True, which='major', color='dimgrey', linestyle='--')

#pA2.plot(ENU_marf5.Hz,ENU_marf5.a2,'g-')
#pA2.plot(NED_marf5.Hz,NED_marf5.a2,'r-')
pA2.plot(NWU_marf5.Hz,NWU_marf5.a2,'m-')
pA2.plot(Fq,A2[icdip,:],'b-')
#pA2.text(-.1,0.5,'A2',fontsize=14)
pA2.set_ylabel('A2')
pA2.set_ylim([-1,1])
pA2.set_yticks([-1,0,1])
pA2.set_xlim([0,1])
pA2.grid(b=True, which='major', color='dimgrey', linestyle='--')

#pB2.plot(ENU_marf5.Hz,ENU_marf5.b2,'g-')
#pB2.plot(NED_marf5.Hz,NED_marf5.b2,'r-')
pB2.plot(NWU_marf5.Hz,NWU_marf5.b2,'m-')
pB2.plot(Fq,B2[icdip,:],'b-')
#pB2.text(0.8,0.5,'B2',fontsize=16)
pB2.set_ylabel('B2')
pB2.set_ylim([-1,1])
pB2.set_yticks([-1,0,1])
pB2.set_xlim([0,1])
pB2.grid(b=True, which='major', color='dimgrey', linestyle='--')
pB2.set_xlabel('Frequency (Hz)')

#%% Wave Trains
#
#
p2heave=r'G:\My Drive\00_CalWave\01_CW_Demo\07_Partners Exchange\18_MarineLabs\EWTEC\march_CS2_5hz.csv'

mar5hz = pd.read_csv(p2heave,skiprows=1,header=None,names=['UTC','z'])

# NOT USED
# FFT of CoastScout hi-freq data
# import scipy.fftpack
# z=mar5hz.z.to_numpy()
# zfft = scipy.fftpack.fft(z,len(z))
# zpsd = np.abs(zfft)**2
# fftfreq = scipy.fftpack.fftfreq(len(z),0.2)
# i = fftfreq>0

# fig = plt.figure(999)
# plt.plot(1/fftfreq[i],zpsd[i],'g-')
# plt.xlim([0,25])


#%%
z = mar5hz.z[0:5000:1]

from scipy.signal import find_peaks
height = None
threshold = 0.001
distance = 5*2# distance in samples...5Hz data, 2sec wave resolution
prominence = None
peaks = find_peaks(z,height=height,threshold=threshold,distance=distance,prominence=prominence)

npeaks = find_peaks(-z,height=height,threshold=threshold,distance=distance,prominence=prominence)

fig=plt.figure(4)
plt.suptitle('Storm March 18',fontsize=18,weight='bold')
plt.plot(z)
plt.plot(peaks[0],z[peaks[0]],'g*')
plt.plot(npeaks[0],z[npeaks[0]],'r*')

#%%
wavei = list(zip(npeaks[0][:],peaks[0][:]))

waveheights=[]

for i,wi in enumerate(wavei):
    waveheights.append(z[wavei[i][1]]-z[wavei[i][0]])

bins = np.linspace(0,2,100)
plt.hist(waveheights,bins=bins)
plt.title('Wave Height Histogram (from 5Hz z timeseries)')
plt.ylabel('Wave Counts')
plt.xlabel('individual wave height')



