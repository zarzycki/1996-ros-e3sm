#!/usr/bin/env python
# coding: utf-8

# # Creates Figure 7

import netCDF4 as nc
import numpy as np
import xarray as xr
import datetime
from datetime import timedelta
#import xesmf as xe
#import cftime

import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point

# Do we have full LENS data or just a ncks subset?
full_data=False

#loading data
if full_data:
    mon_avg = xr.open_dataset('/Users/cmz5202/NetCDF/CESMLENS/ens_mean.nc')
else:
    # ncks -d time,911 -d time,1931 -d lev,29 /Users/cmz5202/NetCDF/CESMLENS/ens_mean.nc ens_mean_T_sub.nc
    mon_avg = xr.open_dataset('./anom_data/ens_mean_T_sub.nc')

#creating array of letters for subplot labels
letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)']

## making figure to show delta technique process

projPC = ccrs.PlateCarree()

#setting up plot
fig, axs = plt.subplots(1,3,subplot_kw={'projection': projPC},figsize=(20,7),constrained_layout=True)
#choosing time indices

if full_data:
    t1 = 911  #corresponds to jan 1996
    t2 = 1931 #corresponds to jan 2081 (+4K)
    lev = 29
else:
    # Only have 2 indices and 1 level in chopped data
    t1 = 0
    t2 = 1
    lev = 0

#loop to plot three subplots of two different years and the difference of the two
for ax,i,t,L in zip(axs,[t2,t1,-1],['2081 January Surface Temperature (K)','1996 January Surface Temperature (K)','Delta: January Surface Temperature (K)'],letters):
    #plotting the difference
    if i == -1:
        cont = np.arange(-16,18,1)
        wrap_data, wrap_lon = add_cyclic_point(mon_avg['T'][t2,lev,:,:] - mon_avg['T'][t1,lev,:,:], coord=mon_avg['lon'], axis=1)
        plot = ax.contourf(wrap_lon,mon_avg['lat'],wrap_data, levels=cont,cmap='seismic')
        #adding an = between the last two subplots
        ax.text(-0.07,.48,'=',ha='center',va='center',fontsize=40,fontweight='bold',transform=ax.transAxes)
    #plotting the two years
    else:
        cont = np.arange(245,310,5)
        wrap_data, wrap_lon = add_cyclic_point(mon_avg['T'][i,lev,:,:], coord=mon_avg['lon'], axis=1)
        plot = ax.contourf(wrap_lon,mon_avg['lat'],wrap_data, levels=cont,cmap='magma',extend='min')
        #adding a subtraction sign between the first two subplots
        if i == t1:
            ax.text(-0.07,.48,'\N{MINUS SIGN}',ha='center',va='center',fontsize=40,fontweight='bold',transform=ax.transAxes)
    #adding letters to panels
    props = dict(boxstyle='square', facecolor='white', edgecolor='k')
    ax.text(.03,.95,L, fontsize=30,transform=ax.transAxes,ha='left',va='top',bbox=props)
    # Set extent
    lonW = -140
    lonE = -50
    latS = 15
    latN = 60
    ax.set_extent([lonW, lonE, latS, latN], crs=projPC)
    #adding title
    ax.set_title(t,fontsize=22,ha='center')
    #adding features to plot
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.STATES)
    #adding colorbar to bottom of subplots
    cb = plt.colorbar(plot, ax=ax, anchor=(0.5, 0.2), fraction=1.0, location='bottom')
    cb.ax.tick_params(labelsize='xx-large')
    #cb.ax.set_title('K', fontsize='xx-large')

fig.savefig('FIG_delta_example.pdf', bbox_inches='tight', pad_inches=0)
