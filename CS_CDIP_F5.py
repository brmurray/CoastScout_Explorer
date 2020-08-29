# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 16:02:17 2020

@author: bryan
"""

import pandas as pd
import numpy as np
import datetime
import netCDF4
import calendar
import os
import shutil
import matplotlib.pyplot as plt

#%%
# Open local copy of netCDF archive file
p2archivefile = r'G:\My Drive\00_CalWave\01_CW_Demo\03_Demo_Test_Site\00_Scripps\00_MetOcean_Data\CDIP201_ScrippsNearshore\Archive_dump\CDIP201_feb-may2020.nc'

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