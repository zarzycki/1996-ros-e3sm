#!/usr/bin/env python
# coding: utf-8

# # This notebook makes figures 4, 5, & 6

#loading necessary libraries
import netCDF4 as nc
import numpy as np
import xarray as xr
import datetime
from datetime import timedelta
import xesmf as xe
import pandas as pd
import shapefile as shp
import cartopy.io.shapereader as shpreader
import pytz
from pytz import timezone
from scipy import stats
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sys

### CMZ Command line args
basepath = sys.argv[1]
softpath = sys.argv[2]

#shifting mosart to be in the middle of the time period
shift_MOSART = True

#set equal to True if want to plot lines on figures 
LowRes = False
Reanalysis = True

#array of letters for subplot labels
letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)']

#loading reanalysis discharge data
rean_dis = xr.open_dataset(basepath+'/rean_glofas/1996_sub.nc')
rean_dis = rean_dis.sel(time=slice('1996-01-15', '1996-01-25'),lat=slice(45.,35.),lon=slice(-85.,-70.))

#regridding the reanalysis data to match the model data
ds_out = xr.Dataset({'lat': (['lat'], np.arange(35, 45.1, 0.125)),
                     'lon': (['lon'], np.arange(-85, -69.9, .125)),
                    })
regridder_dis = xe.Regridder(rean_dis, ds_out, 'bilinear')
    
rean_dis = regridder_dis(rean_dis)


#loading model data
procpath = basepath+'/arp-proc/'
ds_mos = xr.open_dataset(procpath+'/ds_mos_Control_ens.nc')
ds_mos_lr =  xr.open_dataset(procpath+'/Low_Res_mos_slice.nc')


#function that creates an array of time values given a start date, end date, and timestep
def datetime_range(start, end, delta):
    current = start
    if not isinstance(delta, timedelta):
        delta = timedelta(**delta)
    while current < end:
        yield current
        current += delta

#creating time arrays with datetime objects
        
start = datetime.date(1996,1,15)
end = datetime.date(1996,1,26)
time = []
for dt in datetime_range(start, end, {'days':1}):
    time.append(dt)
        
if shift_MOSART == False:
    time_m = ds_mos['time'].values
    time_m = time_m.tolist()
    Time_mos = []
    for i in range(len(ds_mos['time'])):
        Time_mos.append(datetime.datetime.strptime(time_m[i], '%Y-%m-%d-%H'))
if shift_MOSART == True:
    T_mos = []
    start_mos = datetime.datetime(1996,1,16,0)
    end_mos = datetime.datetime(1996,1,24,12)
    Time_mos = []
    for dt in datetime_range(start_mos, end_mos, {'hours':12}):
        Time_mos.append(dt)
    for i in range(len(Time_mos)):
        T_mos.append(Time_mos[i].strftime('%Y-%m-%d-%H'))

#shifting time in dataset 
if shift_MOSART == True:
    for x in [ds_mos,ds_mos_lr]:
        x['time'] = T_mos
        
Time_M = []
for i in np.arange(0,len(Time_mos),2):
    Time_M.append(Time_mos[i].strftime('%b %d'))
    if i != 16:
        Time_M.append('')


# **STATION DATA Information**
# 
# cono = USGS 01578310 SUSQUEHANNA RIVER AT CONOWINGO, MD Latitude 39°39'28.4",   Longitude 76°10'28.0"
# 
# mar = USGS 01576000 Susquehanna River at Marietta, PA Latitude 40°03'16",   Longitude 76°31'52" 
# 
# una = USGS 01500500 SUSQUEHANNA RIVER AT UNADILLA NY Latitude 42°19'17.4",   Longitude 75°18'59.0"
# 
# sun = USGS 01554000 Susquehanna River at Sunbury, PA Latitude 40°50'04",   Longitude 76°49'37" 
# 
# wil = USGS 01551500 WB Susquehanna River at Williamsport, PA Lat 41°14'10", long 76°59'49"
# 
# che = USGS 01512500 CHENANGO RIVER NEAR CHENANGO FORKS NY Lat 42°13'05", long 75°50'54"
# 
# tio = USGS 01518700 Tioga River at Tioga Junction, PA 41°57'09", long 77°06'56"
# 


#loading USGS station data
cono = pd.read_csv(softpath+"/River_obs/Cono_dis.txt", sep='\t')
mar = pd.read_csv(softpath+"/River_obs/mar_dis.txt", sep='\t')
una = pd.read_csv(softpath+"/River_obs/una_dis.txt", sep='\t')
sun = pd.read_csv(softpath+"/River_obs/sun_dis.txt", sep='\t')
wil = pd.read_csv(softpath+"/River_obs/wil_dis.txt", sep='\t')
che = pd.read_csv(softpath+"/River_obs/che_dis.txt", sep='\t')
tio = pd.read_csv(softpath+"/River_obs/tio_dis.txt", sep='\t')


#changing the time format to be more managable (string to datetime)
station = [cono,mar,sun,una,wil,che,tio]
for s in station:
    index = np.arange(0,len(s['datetime']))
    for i in index:
        s['datetime'][i] = datetime.datetime.strptime(s['datetime'][i], '%Y-%m-%d %H:%M')
        #changing from EST to UTC
        s['datetime'][i] = s['datetime'][i].astimezone(timezone('UTC'))
        
        #multiplied by .0283... to change from ft^3 to m^3/s
        s[s.columns[4]][i] = s[s.columns[4]][i]*0.0283168


#loading model lat and lon
Lat = ds_mos.coords['lat']
Lon = ds_mos.coords['lon']
Lat_lr = ds_mos_lr['lat']
Lon_lr = ds_mos_lr['lon']
Lat_r = rean_dis['lat']
Lon_r = rean_dis['lon']
# creating function to find the closest value to a certain lat or lon
def closest(Lat, K): 
      return Lat[min(range(len(Lat)), key = lambda i: abs(Lat[i]-K))]

#station lat on lons
stat_lats = [39.657889,40.054444,42.3215,40.834444,41.236111, 42.218056,41.9525]
stat_lons = [-76.174444+360 , -76.531111+360 , -75.316389+360,-76.826944+360, -76.996944+360, -75.848333+360, -77.115556+360]
stat_lonsR = [-76.174444 , -76.531111 , -75.316389,-76.826944, -76.996944, -75.848333, -77.115556]

#creating empty arrays
lats = []
lons = []
dis = []
lats_r = []
lons_r = []
dis_r = []

#finding the lats and lons that are closest to station data
for x,y,yr in zip(stat_lats, stat_lons, stat_lonsR):
    stat_latn = np.where(Lat == closest(Lat,x))
    lats.append(stat_latn)
    stat_lonn = np.where(Lon == closest(Lon,y))
    lons.append(stat_lonn)
    M_dis = ds_mos['RIVER_DISCHARGE_OVER_LAND_LIQ'][:,stat_latn[0],stat_lonn[0]]
    M_dis = M_dis.squeeze()
    dis.append(M_dis)
    
    stat_latr = np.where(Lat_r == closest(Lat_r,x))
    lats_r.append(stat_latr)
    stat_lonr = np.where(Lon_r == closest(Lon_r,yr))
    lons_r.append(stat_lonr)
    R_dis = rean_dis['dis'][:,stat_latr[0],stat_lonr[0]]
    R_dis = R_dis.squeeze()
    dis_r.append(R_dis)
    
#assigning names to the model data to the associated station name and variable
Mcono_dis,Mmar_dis,Muna_dis,Msun_dis,Mwil_dis,Mche_dis,Mtio_dis = dis
Rcono_dis,Rmar_dis,Runa_dis,Rsun_dis,Rwil_dis,Rche_dis,Rtio_dis = dis_r


#doing the same thing as above but for the low res if it is set to True
if LowRes == True:
    Lat = ds_mos.coords['lat']
    Lon = ds_mos.coords['lon']
    Lat_lr = ds_mos_lr['lat']
    Lon_lr = ds_mos_lr['lon']
    Lat_r = rean_dis['lat']
    Lon_r = rean_dis['lon']
    # creating function to find the closest value to a certain lat or lon
    def closest(Lat, K): 
          return Lat[min(range(len(Lat)), key = lambda i: abs(Lat[i]-K))]

    #station lat on lons
    stat_lats = [39.657889,40.054444,42.3215,40.834444,41.236111, 42.218056,41.9525]
    stat_lons = [-76.174444+360 , -76.531111+360 , -75.316389+360,-76.826944+360, -76.996944+360, -75.848333+360, -77.115556+360]
    stat_lonsR = [-76.174444 , -76.531111 , -75.316389,-76.826944, -76.996944, -75.848333, -77.115556]

    lats4= []
    lons4 = []
    dis4 = []
    lats_lr = []
    lons_lr = []
    dis_lr = []
    #finding the lats and lons that are closest to station data
    for x,y,yr in zip(stat_lats, stat_lons, stat_lonsR):

        stat_latlr = np.where(Lat_lr == closest(Lat_lr,x))
        lats_lr.append(stat_latlr)
        stat_lonlr = np.where(Lon_lr == closest(Lon_lr,yr))
        lons_lr.append(stat_lonlr)
        LR_dis = ds_mos_lr['RIVER_DISCHARGE_OVER_LAND_LIQ'][:,stat_latlr[0],stat_lonlr[0]]
        LR_dis = LR_dis.squeeze()
        dis_lr.append(LR_dis)
        
        
    #assigning names to the model data to the associated station name and variable
    Lcono_dis,Lmar_dis,Luna_dis,Lsun_dis,Lwil_dis,Lche_dis,Ltio_dis = dis_lr


# ## Figure 5

#plots the discharge for the reanalysis data, station data, model data, and low res data on one timeseries plot

#creating lists of the model and reanalysis data for each of the stations
s_dis = [Mcono_dis, Mmar_dis, Msun_dis,Mwil_dis,Mtio_dis,Mche_dis,Muna_dis ]  
r_dis = [Rcono_dis, Rmar_dis, Rsun_dis,Rwil_dis,Rtio_dis,Rche_dis,Runa_dis]

#list of station names
str_stat= ['Conowingo', 'Marietta', 'Sunbury', 'Williamsport', 'Tioga','Chenango','Unadilla']
#list of station data cut to where is ends at a similar time to the model
station = [cono[76:],mar[40:],sun,wil,tio[78:],che,una[39:]]

#setting up plot
fig = plt.figure(figsize=(20,13))
ax1 = plt.subplot2grid(shape=(4,4), loc=(0,0), colspan=2)
ax2 = plt.subplot2grid((4,4), (0,2), colspan=2)
ax3 = plt.subplot2grid((4,4), (1,0), colspan=2)
ax4 = plt.subplot2grid((4,4), (1,2), colspan=2)
ax5 = plt.subplot2grid((4,4), (2,0), colspan=2)
ax6 = plt.subplot2grid((4,4), (2,2), colspan=2)
ax7 = plt.subplot2grid((4,4), (3,1), colspan=2)
plt.subplots_adjust(hspace=0.35,wspace=0.35)
axs = [ax1,ax2,ax3,ax4,ax5,ax6,ax7]

#looping through stations and plotting
for s,d,r,t,ax,L in zip(station, s_dis,r_dis,str_stat,axs,letters):
    p1, = ax.plot(s['datetime'],s[s.columns[4]], color='b',  label='Station') #plotting USGS data
    p2, = ax.plot(Time_mos, d,color='r', label='Model')  #plotting model data
    if Reanalysis == True:
        p3, = ax.plot(rean_dis['time'][1:10], r[1:10], color='k', label='Reanalysis') #plotting reanalysis
        lns=[p1,p2,p3]
    if Reanalysis == False:
        lns=[p1,p2]
    if ax == ax1: #putting legend on first subplot
        ax.legend(handles=lns, loc=2,bbox_to_anchor=(0.08,.51,0.5, 0.5),fontsize='large')
    #adding subplot labels
    props = dict(boxstyle='square', facecolor='white', edgecolor='k')
    ax.text(.013,.95,L, fontsize='x-large',transform=ax.transAxes,ha='left',va='top',bbox=props)
    #adding to plot
    ax.set_ylabel('Discharge m${^3}$ s$^{-1}$', fontsize='large')
    ax.margins(x=0)
    ax.set_xticks(Time_mos)
    ax.set_xticklabels(Time_M, fontsize='large')
    ax.set_title(t, fontsize='xx-large')
    
fig.savefig('discharge_timeseries_gauges.pdf', bbox_inches='tight', pad_inches=0)


#creating empty arrays
model_max_time = []
obs_max_time = []
#arrays of model and statio data at each station
model_stat = [Mcono_dis, Mmar_dis, Msun_dis,Mwil_dis, Mtio_dis,Mche_dis,Muna_dis]
obs_stat = [cono[cono.columns[4]], mar[mar.columns[4]], sun[sun.columns[4]],wil[wil.columns[4]],
            tio[tio.columns[4]], che[che.columns[4]],una[una.columns[4]]]
obs = [cono,mar,sun,wil,tio,che,una]

#loop too find the time of maximum discharge for the station and the model 
for m,o,t in zip(model_stat, obs_stat,obs):
    k = []
    V = m
    vmax = np.nanmax(V.values)
    vmax_ll = np.where(V == vmax)
    vmax_ll = V[vmax_ll[0]]
    vmax_time = vmax_ll.coords['time'].values
    model_max_time.append(vmax_time[0])

    omax = np.nanmax(o.values)
    omax_l = np.where(o == omax)
    omax_ll = o[omax_l[0][-1]]
    omax_time = t['datetime'][omax_l[0]]
    k.append(t['datetime'][omax_l[0][-1]])
    obs_max_time.append(k)


#changing the format of the dates so they match 
Time_max = []
Model_max_time = []
for i in range(len(model_max_time)):
    Time_max.append(datetime.datetime.strptime(model_max_time[i], '%Y-%m-%d-%H'))
    Model_max_time.append(Time_max[i].replace(tzinfo=pytz.UTC))

#taking the difference between the two time of maximums 
time_diff = []
Time_diff = []
for i in range(len(Model_max_time)):
    x = Model_max_time[i] - obs_max_time[i][0]
    
    #changing to seconds 
    y = x.total_seconds() 
    #change to hours to plot
    Time_diff.append(y/3600)


# ## Figure 6

##plotting the difference in the time of maximums for the model/station for each station

stations = ['Cono','Mar','Sun','Wil','Tio','Che','Una']
#setting up figure
fig = plt.figure(figsize=(10,3.5))
ax = fig.add_subplot(111)
#plotting time difference 
ax.plot(Time_diff, linestyle='--', color='k')
ax.plot(Time_diff, 'o', color='k')
#plotting zero line
ax.plot(np.zeros(len(Time_diff)), color='b',linestyle='--',) 
#setting up labels
ax.set_xticks(range(7))
ax.set_xticklabels(stations, fontsize='large')
ax.set_ylabel('Time Difference (hrs)', fontsize='large')

fig.savefig('discharge_err_mouth_to_head.pdf', bbox_inches='tight', pad_inches=0)


#defining a function to create an outline of the SRB
def susq_outline(ax):
    sf = shp.Reader("./shapes/srb.shp")
    shape = sf.shapeRecord()
    for i in range(len(shape.shape.parts)):
        i_start = shape.shape.parts[i]
        if i==len(shape.shape.parts)-1:
            i_end = len(shape.shape.points)
        else:
            i_end = shape.shape.parts[i+1]
        x = [i[0] for i in shape.shape.points[i_start:i_end]]
        y = [i[1] for i in shape.shape.points[i_start:i_end]]
        ax.plot(x, y, transform = ccrs.AlbersEqualArea(standard_parallels=(40,43),central_latitude=39.0, central_longitude=-77.0),color='k')
        #ax.gridlines()
        ax.add_feature(cfeature.BORDERS, linewidth = 1.5, edgecolor='grey')
        ax.add_feature(cfeature.COASTLINE, linewidth = 1.5, edgecolor='grey')
        ax.add_feature(cfeature.STATES, linewidth = 1.5, edgecolor='grey')
        g = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
        g.top_labels = False
        g.right_labels = False


# ## Figure 4

##plots where the stations are located in the basin

#setting up plot
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
c = 'r'

#lat and lons of stations
stat_lats = [39.657889,40.054444,42.3215,40.834444,41.236111, 42.218056,41.9525]
stat_lons = [-76.174444+360 , -76.531111+360 , -75.316389+360,-76.826944+360, -76.996944+360, -75.848333+360, -77.115556+360]
#labels for stations
stat_lab = ['1','2','7','3','4','6','5']
#getting shapefile to plot Susq. river
shpfilename = softpath+'/shapes/ne_10m_rivers_lake_centerlines.shp'
shpfile = shpreader.Reader(shpfilename)

#This loop iterates over each river within the dataset
for rec in shpfile.records():
    name = rec.attributes['name']
    
    #This if statement plots both 'Susquehanna' and 'W. Branch Susquehanna' as blue rivers
    try :
        if 'Susquehanna' in name:
            ax.add_geometries([rec.geometry], crs = ccrs.PlateCarree(), facecolor = 'none', edgecolor = 'blue')
    except TypeError:
        del(name)
#loop plots each model point and station location with a label
for x,xx,y,yy,l in zip(stat_lats,lats,stat_lons,lons,stat_lab):
    if x == stat_lats[0]:
        ax.scatter(ds_mos['lon'][yy],ds_mos['lat'][xx], marker='*', color='k', s=60, label='Nearest Model Point')
        ax.scatter(y-360,x, marker='*', color=c, s=60, label='USGS Station')
        ax.annotate(l, (y-359.96,x+.08), weight='bold', size=14) #adding number label
    else:
        ax.scatter(ds_mos['lon'][yy],ds_mos['lat'][xx], marker='*', color='k', s=60)
    ax.scatter(y-360,x, marker='*', color=c, s=60)
    ax.annotate(l, (y-359.96,x+.08), weight='bold', size=14) 
#adding to plot
ax.set_extent([-79,-74.5,39,43])
ax.legend(loc='lower left', fontsize=13)
susq_outline(ax)

fig.savefig('map_usgs.pdf', bbox_inches='tight', pad_inches=0)

