#!/usr/bin/env python
# coding: utf-8

# # Creates Figure 1, Figure 3, and data for Table 1

# In[1]:


#importing libraries
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import xesmf as xe
import datetime 
from datetime import timedelta
import metpy
import metpy.calc as mpcalc
from metpy.units import units
from scipy import stats
import pandas as pd


# In[2]:


#creating an array of letters to put in corner of figure panels
letters = ['(a)','(b)','(c)','(d)','(e)']


# In[3]:


#loading high res model data
ds_eam = xr.open_dataset('ds_eam_Control_ens.nc')
ds_elm = xr.open_dataset('ds_elm_Control_ens.nc')

#loading low res model data
Ds_eam_lr = xr.open_dataset('Low_Res_eam_slice.nc')
Ds_elm_lr = xr.open_dataset('Low_Res_elm_slice.nc')


# In[4]:


#set Low_Res = True if want to plot low res lines on plots
Low_Res = False
#set Rean = True if you want reanalysis lines on plots
Rean = True


# In[5]:


#creating an array of time values for the xtick labels in the graphs
def datetime_range(start, end, delta):
    current = start
    if not isinstance(delta, timedelta):
        delta = timedelta(**delta)
    while current < end:
        yield current
        current += delta

start = datetime.date(1996,1,15)
end = datetime.date(1996,1,24)
time = []
for dt in datetime_range(start, end, {'days':1}):
    time.append(dt)
    
#changing time array in Month Day format for the xtick labels
time_s = []
for x in range(len(time)):
    time_s.append(time[x].strftime('%b %d'))
    
#creating time arrays for the elm and eam files (changing string to datetime)
time_m = ds_eam['time'].values
time_m = time_m.tolist()
Time_eam = []
for i in range(len(ds_eam['time'])):
    Time_eam.append(datetime.datetime.strptime(time_m[i], '%Y-%m-%d-%H'))
    
time_m2 = ds_elm['time'].values
time_m2 = time_m2.tolist()
Time_elm = []
for i in range(len(ds_elm['time'])):
    Time_elm.append(datetime.datetime.strptime(time_m2[i], '%Y-%m-%d-%H'))


# In[6]:


#loading reanalysis data
path = '/gpfs/group/cmz5202/default/arp5873/JRA_sfc/'
temperature = xr.open_dataset(path+ 'JRA.h1.1996.T2M.nc')
temperature = temperature.sel(time=slice('1996-01-15', '1996-01-22'), lat=slice(30.,50.), lon=slice(265.,300.))
u_wind = xr.open_dataset(path+ 'JRA.h1.1996.U10.nc')
u_wind = u_wind.sel(time=slice('1996-01-15', '1996-01-22'), lat=slice(30.,50.), lon=slice(265.,300.))
v_wind = xr.open_dataset(path+ 'JRA.h1.1996.V10.nc')
v_wind = v_wind.sel(time=slice('1996-01-15', '1996-01-22'), lat=slice(30.,50.), lon=slice(265.,300.))
humidity = xr.open_dataset(path+ 'JRA.h1.1996.RHREFHT.nc')
humidity = humidity.sel(time=slice('1996-01-15', '1996-01-22'), lat=slice(30.,50.), lon=slice(265.,300.))
windspeed = metpy.calc.wind_speed(u_wind['U10'], v_wind['V10'])

#regriddiing reanalysis data to match the model grid
ds_out = xr.Dataset({'lat': (['lat'], np.arange(35, 45.1, 0.125)),
                     'lon': (['lon'], np.arange(-85, -69.875, .125)),
                    })
regridder_temp = xe.Regridder(temperature, ds_out, 'bilinear')
regridder_climate = xe.Regridder(Ds_eam_lr, ds_out, 'bilinear')

Temp = temperature['T2M'] - 273  #change to C
TREFHT = regridder_temp(Temp)
WSPD = regridder_temp(windspeed)
hum = humidity['RHREFHT']
RH = regridder_temp(hum)
DP = mpcalc.dewpoint_from_relative_humidity(TREFHT*units('degC'), RH)  #calculating dew point

#regridding low res data to match high res model data
if Low_Res == True:
    ds_eam_lr = regridder_climate(Ds_eam_lr)
    ds_elm_lr = regridder_climate(Ds_elm_lr)


# In[7]:


##This cell finds the closest point to a given ASOS station and saves the data to a new variable

#loading model lat and lon
Lat = ds_eam.coords['lat']
Lon = ds_eam.coords['lon']

# creating function to find the closest value to a certain lat or lon
def closest(Lat, K): 
      return Lat[min(range(len(Lat)), key = lambda i: abs(Lat[i]-K))]

#station lat on lons
stat_lats = [41.33347,42.20856,40.21714,41.24333,40.85]
stat_lons = [ -75.5803 , -75.97972,-76.85147,-76.92167,-77.85]
#station lats in same format as the low res data
stat_lons_lr = [ -75.5803+360 , -75.97972+360,-76.85147+360,-76.92167+360,-77.85+360]

#creating empty arrays to fill in
lats = []
lons = []
temps = []
rh = []
wnd_spd = []
dwpt = [] 

#loading and creating the same things for low res data if set to true
if Low_Res == True:
    Lat_lr = ds_eam_lr.coords['lat']
    Lon_lr = ds_eam_lr.coords['lon']
    lats_lr = []
    lons_lr = []
    temps_lr = []
    rh_lr = []
    wnd_spd_lr = []
    dwpt_lr = [] 
    
#finding the lats and lons from the model that are closest to the stations and saving them in an array
for x,y,yl in zip(stat_lats, stat_lons, stat_lons_lr):
    #finding closest lats
    stat_latn = np.where(Lat == closest(Lat,x))
    lats.append(stat_latn)
    #finding closest lons
    stat_lonn = np.where(Lon == closest(Lon,y))
    lons.append(stat_lonn)
    
    #getting variable data at the above lats and lons
    M_temp = ds_eam['TREFHT'][:,stat_latn[0],stat_lonn[0]]
    M_temp = M_temp.squeeze()
    Relhum = ds_eam['RHREFHT'][:,stat_latn[0],stat_lonn[0]]
    Relhum = Relhum.squeeze()
    M_wsd = ds_eam['WSPDSRFAV'][:,stat_latn[0],stat_lonn[0]]
    M_wsp = M_wsd.squeeze()
    
    #adding unit data
    temper = np.array(M_temp) * units.degC
    rel_h = np.array(Relhum) * units.percent
    dewpoint = mpcalc.dewpoint_from_relative_humidity(temper, rel_h)
    
    #appending the data to the empty arrays
    rh.append(Relhum)
    temps.append(M_temp)
    wnd_spd.append(M_wsp)
    dwpt.append(dewpoint)
    
    #doing the same things for the low res data if true
    if Low_Res == True:
        stat_latl = np.where(Lat_lr == closest(Lat_lr,x))
        lats_lr.append(stat_latl)
        stat_lonl = np.where(Lon_lr == closest(Lon_lr,yl))
        lons_lr.append(stat_lonl)
        L_temp = ds_eam_lr['TREFHT'][:,stat_latl[0],stat_lonl[0]]
        L_temp = L_temp.squeeze()
        Relhum_lr = ds_eam_lr['RHREFHT'][:,stat_latl[0],stat_lonl[0]]
        Relhum_lr = Relhum_lr.squeeze()
        L_wsd = ds_eam_lr['WSPDSRFAV'][:,stat_latl[0],stat_lonl[0]]
        L_wsp = L_wsd.squeeze()
        temper_lr = np.array(L_temp) * units.degC
        rel_h_lr = np.array(Relhum_lr) * units.percent
        dewpoint_lr = mpcalc.dewpoint_from_relative_humidity(temper_lr, rel_h_lr)
        rh_lr.append(Relhum_lr)
        temps_lr.append(L_temp)
        wnd_spd_lr.append(L_wsp)
        dwpt_lr.append(dewpoint_lr)
        
#assigning names to the model data with the associated station name and variable
Mavp_temp,Mbgm_temp,Mcxy_temp,Mipt_temp,Munv_temp = temps
Mavp_rh,Mbgm_rh,Mcxy_rh,Mipt_rh,Munv_rh = rh
Mavp_wsd,Mbgm_wsd,Mcxy_wsd,Mipt_wsd,Munv_wsd = wnd_spd
Mavp_dwpt,Mbgm_dwpt,Mcxy_dwpt,Mipt_dwpt,Munv_dwpt = dwpt


if Low_Res== True:
    Lavp_temp,Lbgm_temp,Lcxy_temp,Lipt_temp,Lunv_temp = temps_lr
    Lavp_rh,Lbgm_rh,Lcxy_rh,Lipt_rh,Lunv_rh = rh_lr
    Lavp_wsd,Lbgm_wsd,Lcxy_wsd,Lipt_wsd,Lunv_wsd = wnd_spd_lr
    Lavp_dwpt,Lbgm_dwpt,Lcxy_dwpt,Lipt_dwpt,Lunv_dwpt = dwpt_lr


# In[8]:


##This cell does the same as above but for the reanalysis data set rather than the model data

#loading model lat and lon
Lat = ds_eam.coords['lat']
Lon = ds_eam.coords['lon']

# creating function to find the closest value to a certain lat or lon
def closest(Lat, K): 
      return Lat[min(range(len(Lat)), key = lambda i: abs(Lat[i]-K))]

#station lat on lons
stat_lats = [41.33347,42.20856,40.21714,41.24333,40.85]
stat_lons = [ -75.5803 , -75.97972,-76.85147,-76.92167,-77.85]

#creating empty array
lats = []
lons = []
Rtemps = []
Rwinds = []
Rdewps = []
Rrhs = []

#finding the lats and lons that are closest to station data
for x,y in zip(stat_lats, stat_lons):
    stat_latn = np.where(Lat == closest(Lat,x))
    lats.append(stat_latn)
    stat_lonn = np.where(Lon == closest(Lon,y))
    lons.append(stat_lonn)
    
    R_temp = TREFHT[:,stat_latn[0],stat_lonn[0]]
    R_temp = R_temp.squeeze()
    temper = np.array(R_temp) * units.degC
    Rtemps.append(R_temp)
    R_wind = WSPD[:,stat_latn[0],stat_lonn[0]]
    R_wind = R_wind.squeeze()
    Rwinds.append(R_wind)
    R_rh = RH[:,stat_latn[0],stat_lonn[0]]
    R_rh = R_rh.squeeze()
    R_rh = np.array(R_rh) * units.percent
    Rrhs.append(R_rh)
    R_dp = mpcalc.dewpoint_from_relative_humidity(temper, R_rh)
    Rdewps.append(R_dp)

#assigning names to the reanalysis data to the associated station name and variable
Ravp_temp,Rbgm_temp,Rcxy_temp,Ript_temp,Runv_temp = Rtemps
Ravp_wind,Rbgm_wind,Rcxy_wind,Ript_wind,Runv_wind = Rwinds
Ravp_dp,Rbgm_dp,Rcxy_dp,Ript_dp,Runv_dp = Rdewps
Ravp_rh,Rbgm_rh,Rcxy_rh,Ript_rh,Runv_rh = Rrhs


# ## ASOS Stations used and their lat lon coords
# - avp : Wilkes-Barre/Scranton 41.33347 -75.5803
# - bgm : Binghamton/Broome 42.20856 -75.97972
# - cxy : Harrisburg/Capital 40.21714 -76.85147
# - ipt : Williamsport 41.24333 -76.92167
# - unv : State College 40.85 -77.85

# In[9]:


#loading station data
avp = pd.read_csv("./station_obs/asos_avp.txt")
bgm = pd.read_csv("./station_obs/asos_bgm.txt")
cxy = pd.read_csv("./station_obs/asos_cxy.txt")
ipt = pd.read_csv("./station_obs/asos_ipt.txt")
unv = pd.read_csv("./station_obs/asos_unv.txt")


# In[10]:


#replacing missing data with nan and converting temp from f to C

#creating array of station data
stations = [avp,bgm,cxy,ipt,unv]

#looping through stations
for s in stations:
    #replacing missing with nans
    s['tmpf'] = s['tmpf'].replace('M', np.nan)
    s['sknt'] = s['sknt'].replace('M', np.nan)
    s['relh'] = s['relh'].replace('M', np.nan)
    #assigning units
    #s['tmpf'].attrs["units"]="deg C"
    #changing string to float
    s['tmpf']= [float(i) for i in s['tmpf']]
    s['sknt'] = [float(i) for i in s['sknt'] ]
    s['relh'] = [float(i) for i in s['relh'] ]
    #converting F to C
    s['tmpf'] = [(i-32.0) * (5.0/9.0) for i in s['tmpf']]
    s['sknt'] = s['sknt']*0.514444
    #assigning units
    temp = np.array(s['tmpf']) * units.degC
    rh = np.array(s['relh']) * units.percent
    #calculating dewpoints
    dewpoint = mpcalc.dewpoint_from_relative_humidity(temp, rh)
    s['dwpt'] = dewpoint


# In[11]:


#changing station time to datetime object to match model
for s in stations:
    index = np.arange(0,len(s['valid']))
    for i in index:
        s['valid'][i] = datetime.datetime.strptime(s['valid'][i], '%Y-%m-%d %H:%M')
       


# In[12]:


##choosing station times to match model

#creating empty datframes with the same column info as station data
avp_mt = avp[0:0]
bgm_mt = bgm[0:0]
cxy_mt = cxy[0:0]
ipt_mt = ipt[0:0]
unv_mt = unv[0:0]
stations2 = []

#filling new dataframes with the station times that match model times
for s in stations:
    for x in np.arange(0,len(Time_eam)):
        for y in np.arange(0,len(s['valid'])):
            if Time_eam[x] == s['valid'][y]:
                if s['station'][s.index[-1]] == 'AVP':
                    avp_mt.loc[y] = s.iloc[y]
                if s['station'][s.index[-1]] == 'BGM':
                    bgm_mt.loc[y] = s.iloc[y]
                if s['station'][s.index[-1]] == 'CXY':
                    cxy_mt.loc[y] = s.iloc[y]
                if s['station'][s.index[-1]] == 'IPT':
                    ipt_mt.loc[y] = s.iloc[y]
                if s['station'][s.index[-1]] == 'UNV':
                    unv_mt.loc[y] = s.iloc[y]


# In[13]:


##plots time series of station and model data for temp, RH, windspeed, and dewpoint

#creating array of variables for title and y-axis label for graph
title = ['Temperature', 'Relative Humidity', 'Wind Speed', 'Dewpoint']
ylabel = ['Temperature ($^\circ$C)', 'Relative Humidity (%)', 'Wind Speed (m s$^{-1}$)', 'Dewpoint ($^\circ$C)']

#plotting high res and low res line if True
if Low_Res == True:
    for s,t,r,w,d,tt,rr,ww,dd in zip([avp_mt,bgm_mt,ipt_mt,unv_mt,cxy_mt], [Mavp_temp,Mbgm_temp,Mipt_temp,Munv_temp,Mcxy_temp],
                         [Mavp_rh,Mbgm_rh,Mipt_rh,Munv_rh,Mcxy_rh],[Mavp_wsd,Mbgm_wsd,Mipt_wsd,Munv_wsd,Mcxy_wsd],
                         [Mavp_dwpt,Mbgm_dwpt,Mipt_dwpt,Munv_dwpt,Mcxy_dwpt],[Lavp_temp,Lbgm_temp,Lipt_temp,Lunv_temp,Lcxy_temp],
                         [Lavp_rh,Lbgm_rh,Lipt_rh,Lunv_rh,Lcxy_rh],[Lavp_wsd,Lbgm_wsd,Lipt_wsd,Lunv_wsd,Lcxy_wsd],
                         [Lavp_dwpt,Lbgm_dwpt,Lipt_dwpt,Lunv_dwpt,Lcxy_dwpt]):

        station = [s]*4
        model_temp = t
        model_rh = r
        model_wsp = w
        model_dp = d

        LR_temp = tt
        LR_rh = rr
        LR_wsp = ww
        LR_dp = dd
        fig, axs = plt.subplots(4,1,figsize=(12,15))
        for ax,s_data, m_data,lr_data,S, titl,y in zip(axs,['tmpf','relh','sknt','dwpt'],[model_temp, model_rh, model_wsp, model_dp],[LR_temp, LR_rh, LR_wsp, LR_dp], station, title,ylabel):
            x = S[s_data].interpolate(method='linear')
            ax.plot(S['valid'],x, label='Station')
            ax.plot(Time_eam,m_data, color = 'r', label='Model (Control)')
            ax.plot(Time_eam,lr_data, color = 'g', label='Low Res')
            ax.set_title(S['station'][S.index[0]] + ' Comparison to Model Data - ' + titl)
            ax.set_ylabel(y)
            ax.legend()
            plt.tight_layout()
            
#plotting model and reanalysis
else:
    #looping through model/reanalysis data at the closest model/rean point to station
    for s,t,r,w,d,Rt,Rw,Rd,Rrh in zip([avp_mt,bgm_mt,ipt_mt,unv_mt,cxy_mt], [Mavp_temp,Mbgm_temp,Mipt_temp,Munv_temp,Mcxy_temp],
        [Mavp_rh,Mbgm_rh,Mipt_rh,Munv_rh,Mcxy_rh],[Mavp_wsd,Mbgm_wsd,Mipt_wsd,Munv_wsd,Mcxy_wsd],
        [Mavp_dwpt,Mbgm_dwpt,Mipt_dwpt,Munv_dwpt,Mcxy_dwpt],[Ravp_temp,Rbgm_temp,Ript_temp,Runv_temp,Rcxy_temp],
        [Ravp_wind,Rbgm_wind,Rcxy_wind,Ript_wind,Runv_wind],[Ravp_dp,Rbgm_dp,Rcxy_dp,Ript_dp,Runv_dp],
        [Ravp_rh,Rbgm_rh,Rcxy_rh,Ript_rh,Runv_rh]):
        station = [s]*4
        model_temp = t
        model_rh = r
        model_wsp = w
        model_dp = d
        rean_temp = Rt
        rean_rh = Rrh
        rean_wsp = Rw
        rean_dp = Rd

        fig, axs = plt.subplots(4,1,figsize=(10,14))
        for ax,s_data, m_data,r_data,S, titl,y,L in zip(axs,['tmpf','relh','sknt','dwpt'],[model_temp, model_rh, model_wsp, model_dp], 
                                                      [rean_temp, rean_rh, rean_wsp, rean_dp],station, title,ylabel,letters):
            #interpolating data so the line is smooth
            x = S[s_data].interpolate(method='linear')
            #plotting station data
            ax.plot(S['valid'],x, label='Station')
            #plotting high res model data
            ax.plot(Time_eam,m_data, color = 'r', label='Model')
            #plotting reanalysis if = True
            if Rean == True:
                ax.plot(TREFHT['time'], r_data, color='k', label='Reanalysis')

            #setting up the plot
            ax.set_title(L+ ' '+titl, fontsize='xx-large')
            ax.set_ylabel(y, fontsize='x-large')
            ax.set_xticks(time)
            ax.set_xticklabels(time_s,fontsize='medium')
            ax.tick_params('y', labelsize='x-large')
            ax.tick_params('x', labelsize='large')
            #only putting the legend on the first subplot
            if s_data == 'tmpf':
                ax.legend(loc=0,fontsize='large', ncol=3)
                ax.text(0.01,.87,S['station'][S.index[0]], transform=ax.transAxes,fontweight='bold', fontsize=27)
            plt.tight_layout()
        
        thisStation=S['station'][S.index[0]]
        fig.savefig('timeseries_'+thisStation+'.pdf', bbox_inches='tight', pad_inches=0)


# ## Calculating correlation and bias of model

# In[14]:


#creating time array of datetime64 objects
time_eam = Time_eam
for t in range(len(Time_eam)):
    time_eam[t] = np.datetime64(Time_eam[t])

##choosing the model times that match the station data we previously cut so that they have the same number of values for the correlation calculations
#creating empty arrays
temps2 = []
rh2 = []
wnd_spd2 = []
dwpt2 = [] 
#looping through station and model data
for s,t,r,w,d in zip([avp_mt,bgm_mt,ipt_mt,unv_mt,cxy_mt], [Mavp_temp,Mbgm_temp,Mipt_temp,Munv_temp,Mcxy_temp],
                     [Mavp_rh,Mbgm_rh,Mipt_rh,Munv_rh,Mcxy_rh],
                     [Mavp_wsd,Mbgm_wsd,Mipt_wsd,Munv_wsd,Mcxy_wsd],
                     [Mavp_dwpt,Mbgm_dwpt,Mipt_dwpt,Munv_dwpt,Mcxy_dwpt]):
    stat_time = s['valid'].values
    for i in range(len(Time_eam)):
        #only appending model data where the two times match
        if time_eam[i] in stat_time:
            temps2.append(t[i])
            rh2.append(r[i])
            wnd_spd2.append(w[i])
            dwpt2.append(d[i])


# In[15]:


#appending the correct data at each station to a variable with that station in the name
Mavp_temp2 = temps2[0:len(avp_mt)]
Mbgm_temp2 = temps2[len(avp_mt):len(bgm_mt)+len(avp_mt)]
Mipt_temp2 = temps2[len(bgm_mt)+len(avp_mt):len(bgm_mt)+len(avp_mt)+len(ipt_mt)]
Munv_temp2 = temps2[len(bgm_mt)+len(avp_mt)+ len(ipt_mt):len(bgm_mt)+len(avp_mt)+len(ipt_mt)+len(unv_mt)]
Mcxy_temp2 = temps2[len(bgm_mt)+len(avp_mt)+ len(ipt_mt)+ len(unv_mt):len(bgm_mt)+len(avp_mt)+len(ipt_mt)+len(unv_mt)+ len(cxy_mt)]

Mavp_rh2 = rh2[0:len(avp_mt)]
Mbgm_rh2 = rh2[len(avp_mt):len(bgm_mt)+len(avp_mt)]
Mipt_rh2 = rh2[len(bgm_mt)+len(avp_mt):len(bgm_mt)+len(avp_mt)+len(ipt_mt)]
Munv_rh2 = rh2[len(bgm_mt)+len(avp_mt)+ len(ipt_mt):len(bgm_mt)+len(avp_mt)+len(ipt_mt)+len(unv_mt)]

Mcxy_rh2 = rh2[len(bgm_mt)+len(avp_mt)+ len(ipt_mt)+ len(unv_mt) :len(bgm_mt)+len(avp_mt)+len(ipt_mt)+len(unv_mt)+ len(cxy_mt)]

Mavp_wsd2 = wnd_spd2[0:len(avp_mt)]
Mbgm_wsd2 = wnd_spd2[len(avp_mt):len(bgm_mt)+len(avp_mt)]
Mipt_wsd2 = wnd_spd2[len(bgm_mt)+len(avp_mt):len(bgm_mt)+len(avp_mt)+len(ipt_mt)]
Munv_wsd2 = wnd_spd2[len(bgm_mt)+len(avp_mt)+ len(ipt_mt):len(bgm_mt)+len(avp_mt)+len(ipt_mt)+len(unv_mt)]
Mcxy_wsd2 = wnd_spd2[len(bgm_mt)+len(avp_mt)+ len(ipt_mt)+ len(unv_mt):len(bgm_mt)+len(avp_mt)+len(ipt_mt)+len(unv_mt)+ len(cxy_mt)]

Mavp_dwpt2 = tuple(dwpt2[0:len(avp_mt)])
Mbgm_dwpt2 = tuple(dwpt2[len(avp_mt):len(bgm_mt)+len(avp_mt)])
Mipt_dwpt2 = tuple(dwpt2[len(bgm_mt)+len(avp_mt):len(bgm_mt)+len(avp_mt)+len(ipt_mt)])
Munv_dwpt2 = tuple(dwpt2[len(bgm_mt)+len(avp_mt)+ len(ipt_mt):len(bgm_mt)+len(avp_mt)+len(ipt_mt)+len(unv_mt)])
Mcxy_dwpt2 = tuple(dwpt2[len(bgm_mt)+len(avp_mt)+ len(ipt_mt)+ len(unv_mt) :len(bgm_mt)+len(avp_mt)+len(ipt_mt)+len(unv_mt)+ len(cxy_mt)])


# In[16]:


##only getting the magnitude of the data (currenlty units are attached)

MAvp_dwpt2 = []
MBgm_dwpt2 = []
MIpt_dwpt2 = []
MUnv_dwpt2 = []
MCxy_dwpt2 = []
for x in range(len(Mavp_dwpt2)):
    y = Mavp_dwpt2[x].magnitude
    MAvp_dwpt2.append(y)

for x in range(len(Mbgm_dwpt2)):
    y = Mbgm_dwpt2[x].magnitude
    MBgm_dwpt2.append(y)  
    
for x in range(len(Mipt_dwpt2)):
    y = Mipt_dwpt2[x].magnitude
    MIpt_dwpt2.append(y)
    
for x in range(len(Munv_dwpt2)):
    y = Munv_dwpt2[x].magnitude
    MUnv_dwpt2.append(y)

for x in range(len(Mcxy_dwpt2)):
    y = Mcxy_dwpt2[x].magnitude
    MCxy_dwpt2.append(y)


# In[17]:


#calulates r value and bias for the model vs station data
Bias = []
R_val = []

for s,t,r,w,d,ss in zip([avp_mt,bgm_mt,ipt_mt,unv_mt,cxy_mt], [Mavp_temp2,Mbgm_temp2,Mipt_temp2,Munv_temp2,Mcxy_temp2],
                     [Mavp_rh2,Mbgm_rh2,Mipt_rh2,Munv_rh2,Mcxy_rh2],
                     [Mavp_wsd2,Mbgm_wsd2,Mipt_wsd2,Munv_wsd2,Mcxy_wsd2],[MAvp_dwpt2,MBgm_dwpt2,MIpt_dwpt2,MUnv_dwpt2,MCxy_dwpt2], ['AVP', 'BGM', 'IPT', 'UNV', 'CXY']):

    station = [s]*4
    model_temp = t
    model_rh = r
    model_wsp = w
    model_dp  = d

    for ax, s_data, m_data,S in zip(axs,['tmpf','relh','sknt', 'dwpt'],[model_temp, model_rh, model_wsp, model_dp], station):
        x = S[s_data].interpolate(method='linear')
        
        #appending that station name before the data
        if s_data == 'tmpf':
            Bias.append(str(ss))
            R_val.append(str(ss))
            
        #calculating r and p value for data
        r_val,p_val = stats.pearsonr(x, m_data)
        R_val.append(round(r_val,5))
        
        #calculating bias
        z = m_data-S[s_data]
        Bias.append(round(np.mean(z),5))


# In[18]:


#defining function to help write csvs
from csv import writer
def append_list_as_row(file_name, list_of_elem):
    # Open file in append mode
    with open(file_name, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow(list_of_elem)


# In[19]:


#creating bias and correlation csvs
f = open('model_bias.csv','w')
g =  open('model_correlation.csv','w')
vari_b = ['Bias','Temp', 'RH', 'Wind Speed', 'Dew Point']
vari_c = ['Corellation','Temp', 'RH', 'Wind Speed', 'Dew Point']
        
#only need this when first creating. the files
append_list_as_row('model_bias.csv',vari_b)
append_list_as_row('model_correlation.csv',vari_c)

#adds bias and r value data to file
b_ind = np.arange(0,len(Bias),5)
for x in b_ind:
    append_list_as_row('model_bias.csv', Bias[x:x+5])
    append_list_as_row('model_correlation.csv', R_val[x:x+5])


# ## Spatial plot of station locations

# In[20]:


#creating function to outline SRB on spatial plots
import shapefile as shp
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
        g.xlines = False
        g.ylines = False


# In[21]:


##Plots the stations and closest model point with outline of SRB

#defining colors for station and model point
cc = 'k'
c = 'r'
#setting up plot
fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()},
                         figsize=(10,6))

#getting shapefile for the Susquehanna River
shpfilename = './shapes/ne_10m_rivers_lake_centerlines.shp'
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

#putting SRB outline on plot
susq_outline(ax)
#plotting station locations
ax.scatter(-75.5803,41.33347, marker='^', color=c,s=60, label='ASOS Station')
ax.scatter(-75.97972,42.20856, marker='^', color=c,s=60)
ax.scatter(-76.85147,40.21714, marker='^', color=c,s=60)
ax.scatter(-76.92167,41.24333, marker='^', color=c, s=60)
ax.scatter(-77.85,40.85, marker='^', color=c, s=60)

#setting spatial extent of plot
ax.set_extent([-81.5,-73,38,44.5])
#choosing labels for each point 
labels= [3,5,1,2,4]
stat_label = ['Wilkes-Barre/Scranton (AVP)','Binghamton (BGM)','Harrisburg (CXY)', 'Williamsport (IPT)',
              'State College (UNV)']

#looping through closest model lat and lons and plotting model point
for x,l,s in zip(range(5),labels,stat_label):
    if x !=4:
        ax.scatter(ds_eam['lon'][lons[x]],ds_eam['lat'][lats[x]], marker='^', color=cc, s=60)
        #putting number label next to point
        ax.annotate(l, (ds_eam['lon'][lons[x]]-.09,ds_eam['lat'][lats[x]]+.09),weight='bold', size=14, label=s)
    
    #making it to where only the last point is included in legend
    if x == 4:
        ax.scatter(ds_eam['lon'][lons[x]],ds_eam['lat'][lats[x]], marker='^', color=cc, s=60, label='Nearest Model Point')
        ax.annotate(l, (ds_eam['lon'][lons[x]]-.09,ds_eam['lat'][lats[x]]+.09), weight='bold', size=14)
ax.legend(loc=2, fontsize=13)

fig.savefig('map_asos.pdf', bbox_inches='tight', pad_inches=0)


# In[ ]:




