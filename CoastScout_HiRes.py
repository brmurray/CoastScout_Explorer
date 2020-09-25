"""
CoastScout HiRes

Python tools for playing with CoastScout high-resolution (5Hz) heave data
Focus on buoy-buoy comparsion with CDIP 201 "Scripps Nearshore" off La Jolla, CA

Contact: Bryan Murray, brmurray@mailbox.org

Change Log:
    + 9/24/2020: First commit

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

#%% Import CDIP201 archive covering Feb-Aug 2020

# Open local copy of netCDF archive file
p2archivefile = r'G:\My Drive\00_CalWave\01_CW_Demo\03_Demo_Test_Site\00_Scripps\00_MetOcean_Data\CDIP201_ScrippsNearshore\Archive_dump\CDIP201_feb-aug2020.nc'

ds = netCDF4.Dataset(p2archivefile)
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

# round down to half hour
# CDIP data usual reports time at xx:03:35 or xx:33:45
ns['thh'] = ns.index
ns.thh= ns.thh.apply(lambda ti: ti - (ti%1800))
ns = ns.set_index('thh')

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


#%% Open allCoastScout heave time series and combine into two time ranges
#   20200318T181000Z UTC - 20200319T064000Z UTC
#   20200807T181000Z - 20200808T051000Z
    
all_heave_data_folder = r'C:\Users\bryan\Documents\GitHub\CoastScout\Scripps02HeaveTimeSeries'

import os
heave_folders = os.listdir(all_heave_data_folder)
march_folders = [os.path.join(all_heave_data_folder,folder,r'timeseries.csv') for i,folder in enumerate(heave_folders) if folder[4:6]=='03']
aug_folders = [os.path.join(all_heave_data_folder,folder,r'timeseries.csv') for i,folder in enumerate(heave_folders) if folder[4:6]=='08']
    
march = pd.concat(pd.read_csv(f,skiprows=7,header=None,names=['UTC','z']) for f in march_folders[23:24])
#march['t'] = iso2posix(march.UTC)

aug = pd.concat(pd.read_csv(f,skiprows=7,header=None,names=['UTC','z']) for f in aug_folders)

#%% Open allCoastScout "first 5" files and combine into two time ranges
#   20200318T181000Z UTC - 20200319T064000Z UTC
#   20200807T181000Z - 20200808T051000Z
    
all_f5_data_folder = r'C:\Users\bryan\Documents\GitHub\CoastScout\Scripps02FourierCoefficients'

import os
f5_folders = os.listdir(all_f5_data_folder)
march_f5_folders = [os.path.join(all_f5_data_folder,folder,r'fourier_coefficients.csv') for i,folder in enumerate(f5_folders) if folder[4:6]=='03']
aug_f5_folders = [os.path.join(all_f5_data_folder,folder,r'fourier_coefficients.csv') for i,folder in enumerate(f5_folders) if folder[4:6]=='08']
    
marf5 = pd.concat(pd.read_csv(f,skiprows=9,header=None,index_col=False,names=['Hz','mms','a1','b1','a2','b2']) for f in march_f5_folders[23:24])
#augf5 = pd.concat(pd.read_csv(f,skiprows=9,header=None,index_col=False,names=['Hz','mms','a1','b1','a2','b2']) for f in aug_f5_folders[22:23])

#%% FFT of CoastScout hi-freq data
import scipy.fftpack
z=march.z.to_numpy()
zfft = scipy.fftpack.fft(z,len(z))
zpsd = np.abs(zfft)**2
fftfreq = scipy.fftpack.fftfreq(len(z),0.2)

i = fftfreq>0

plt.plot(1/fftfreq[i],zpsd[i],'g-')
plt.xlim([0,25])

#%% Compare a single half-hour record
# 1584556200 = March 18, 2020 18:30:00 = CDIP201 record 903
# 20200319T064000Z = 1584601200 = CDIP201 record 928
icdip = 928


plt.figure()
plt.stem(marf5.Hz,marf5.mms,linefmt='g-',markerfmt='g.',use_line_collection=True)
plt.stem(Fq,Ed[icdip,:],linefmt='b--',markerfmt='b.',use_line_collection=True)
#plt.plot(marf5.Hz,marf5.mms,'g-')
#plt.plot(Fq,Ed[icdip,:],'b-')


#%% Plot First 5 as a function of Frequency

f, (pEd,pA1,pB1,pA2,pB2) = plt.subplots(5,1,sharex=True,figsize=(10,6),num='CDIP201 vs CoastScout, "First 5" Components')
# Set title
plt.suptitle('CDIP201 vs CoastScout "First 5" Components \n March 19, 2020 0700 GMT',fontsize=12)

pEd.plot(marf5.Hz,marf5.mms,'g-')
pEd.plot(Fq,Ed[icdip,:],'b-')
pEd.set_ylabel('m*m/Hz')
pEd.set_title('Energy Density')
pEd.grid(b=True, which='major', color='dimgrey', linestyle='--')
pEd.legend(['CoastScout','CDIP201'])
pEd.set_xlim([0,1])

pA1.plot(marf5.Hz,marf5.a1,'g-')
pA1.plot(Fq,A1[icdip,:],'b-')
#pA1.text(0.8,0.5,'A1',fontsize=10)
pA1.set_ylabel('A1')
pA1.set_ylim([-1,1])
pA1.set_yticks([-1,0,1])
pA1.set_xlim([0,1])
pA1.grid(b=True, which='major', color='dimgrey', linestyle='--')

pB1.plot(marf5.Hz,marf5.b1,'g-')
pB1.plot(Fq,B1[icdip,:],'b-')
#pB1.text(0.8,0.5,'B1',fontsize=12)
pB1.set_ylabel('B1')
pB1.set_ylim([-1,1])
pB1.set_yticks([-1,0,1])
pB1.set_xlim([0,1])
pB1.grid(b=True, which='major', color='dimgrey', linestyle='--')

pA2.plot(marf5.Hz,marf5.a2,'g-')
pA2.plot(Fq,A2[icdip,:],'b-')
#pA2.text(-.1,0.5,'A2',fontsize=14)
pA2.set_ylabel('A2')
pA2.set_ylim([-1,1])
pA2.set_yticks([-1,0,1])
pA2.set_xlim([0,1])
pA2.grid(b=True, which='major', color='dimgrey', linestyle='--')

pB2.plot(marf5.Hz,marf5.b2,'g-')
pB2.plot(Fq,B2[icdip,:],'b-')
#pB2.text(0.8,0.5,'B2',fontsize=16)
pB2.set_ylabel('B2')
pB2.set_ylim([-1,1])
pB2.set_yticks([-1,0,1])
pB2.set_xlim([0,1])
pB2.grid(b=True, which='major', color='dimgrey', linestyle='--')
pB2.set_xlabel('Frequency (Hz)')



