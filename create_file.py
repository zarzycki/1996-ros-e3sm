#!/usr/bin/env python
# coding: utf-8
##Importing necessary libraries
import netCDF4 as nc
import numpy as np
import xarray as xr
from datetime import datetime
from datetime import timedelta

#defining a function which creates an array of time values which is used for the dates in the new files
def datetime_range(start, end, delta):
        current = start
        if not isinstance(delta, timedelta):
            delta = timedelta(**delta)
        while current < end:
            yield current
            current += delta

CONFIGS = ['RoS-F2010C5-ne0conus30x8-004-control', \
           'RoS-F2010C5-ne0conus30x8-004-PI', \
           'RoS-F2010C5-ne0conus30x8-004-plus4K', \
           'RoS-F2010C5-ne0conus30x8-004-plus3K', \
           'RoS-F2010C5-ne0conus30x8-004-plus2K', \
           'RoS-F2010C5-ne0conus30x8-004-plus1K', \
           'RoS-F2010C5-ne0conus30x8-005-control', \
           'RoS-F2010C5-ne0conus30x8-005-PI', \
           'RoS-F2010C5-ne0conus30x8-005-plus4K', \
           'RoS-F2010C5-ne0conus30x8-005-plus3K', \
           'RoS-F2010C5-ne0conus30x8-005-plus2K', \
           'RoS-F2010C5-ne0conus30x8-005-plus1K', \
           'RoS-F2010C5-ne0conus30x8-006-control', \
           'RoS-F2010C5-ne0conus30x8-006-PI', \
           'RoS-F2010C5-ne0conus30x8-006-plus4K', \
           'RoS-F2010C5-ne0conus30x8-006-plus3K', \
           'RoS-F2010C5-ne0conus30x8-006-plus2K', \
           'RoS-F2010C5-ne0conus30x8-006-plus1K', \
           'RoS-F2010C5-ne30-001-control', \
          ]

TITLES = ['Control-4', \
          'Pre_Ind-4', \
          'Plus4K-4', \
          'Plus3K-4', \
          'Plus2K-4', \
          'Plus1K-4', \
          'Control-5', \
          'Pre_Ind-5', \
          'Plus4K-5', \
          'Plus3K-5', \
          'Plus2K-5', \
          'Plus1K-5', \
          'Control-6', \
          'Pre_Ind-6', \
          'Plus4K-6', \
          'Plus3K-6', \
          'Plus2K-6', \
          'Plus1K-6', \
          'Low_Res', \
         ]

#choose the initialization times
date = ["1996011512","1996011600","1996011612","1996011700","1996011712"]

PATH_TO_E3SM_FILES="/gpfs/group/cmz5202/default/arp5873/E3SM/"

counter = 0
for thisconfig in CONFIGS:

    # ## Loads the model output files for EAM, ELM, and MOSART
    ### choose the configuration you want to use and uncomment it. Each one is a different model run.
    CONFIG = thisconfig  

    #change the title depending on the CONFIG chosen. Ex: Control-4, Control-5 or Plus1K-4 and so on
    title = TITLES[counter]

    #creates the file for the dates & CONFIG chosen
    for ii, DATE in enumerate(date):
        print("Appending "+DATE+"")
        ds_eam_full_single = xr.open_mfdataset(PATH_TO_E3SM_FILES+'/'+CONFIG+'/'+DATE+'/'+CONFIG+'.eam.h0.1996-*_regrid.nc')
        ds_elm_full_single = xr.open_mfdataset(PATH_TO_E3SM_FILES+'/'+CONFIG+'/'+DATE+'/'+CONFIG+'.elm.h0.1996-*_regrid.nc')
        ds_mos_full_single = xr.open_mfdataset(PATH_TO_E3SM_FILES+'/'+CONFIG+'/'+DATE+'/'+CONFIG+'.mosart.h0.1996-*_regrid.nc')
        if CONFIG == 'RoS-F2010C5-ne30-001-control':  #longitudes are in different format for low resolution file so must slice differently
            ds_eam_full_single = ds_eam_full_single.sel(time=slice('1996-01-15', '1996-01-22'), lat=slice(35.,45.), lon=slice(275.,290))
            ds_elm_full_single = ds_elm_full_single.sel(time=slice('1996-01-15', '1996-01-22'), lat=slice(35.,45.), lon=slice(275.,290))
        else: #slicing all high res files
            ds_eam_full_single = ds_eam_full_single.sel(time=slice('1996-01-15', '1996-01-22'), lat=slice(35.,45.), lon=slice(-85.,-70.))
            ds_elm_full_single = ds_elm_full_single.sel(time=slice('1996-01-15', '1996-01-22'), lat=slice(35.,45.), lon=slice(-85.,-70.))
        ds_mos_full_single = ds_mos_full_single.sel(time=slice('1996-01-15', '1996-01-25'), lat=slice(35.,45.), lon=slice(275.,290.))
        ds_eam_full_single = ds_eam_full_single.expand_dims('init')
        ds_elm_full_single = ds_elm_full_single.expand_dims('init')
        ds_mos_full_single = ds_mos_full_single.expand_dims('init')
        if ii == 0:
            ds_eam_full = ds_eam_full_single
            ds_elm_full = ds_elm_full_single
            ds_mos_full = ds_mos_full_single
        else:
            ds_eam_full = xr.concat( (ds_eam_full, ds_eam_full_single) , "init" )
            ds_elm_full = xr.concat( (ds_elm_full, ds_elm_full_single) , "init" )
            ds_mos_full = xr.concat( (ds_mos_full, ds_mos_full_single) , "init" )
        ds_eam_full_single.close
        ds_elm_full_single.close
        ds_mos_full_single.close


    # ## **Creating Soil File**

    # In[5]:


    ###creating soil file with init as a fourth dimension
    # Copy datasets to base files for ease of use with jupyter
    ds_elm = ds_elm_full

    #setting start and end dates of file to match the times of the slice
    start = datetime(1996,1,15,12)
    end = datetime(1996,1,22,23)

    #creating a time array
    time = []
    for dt in datetime_range(start, end, {'hours':6}):
        time.append(dt)
    #converting the time array into a string for the file
    Time = []
    for x in range(len(time)):
        Time.append(time[x].strftime('%Y-%m-%d-%H'))

    #creating file name
    file = title + '_init_soil.nc'
    #opening a file
    ds = nc.Dataset(file,'w', format='NETCDF4')

    #creating dimensions of file
    time = ds.createDimension('time', len(ds_elm['time']))
    init = ds.createDimension('init', len(ds_elm['init']))
    lat = ds.createDimension('lat', len(ds_elm['lat']))
    lon = ds.createDimension('lon', len(ds_elm['lon']))
    levgrnd  = ds.createDimension('levgrnd', len(ds_elm['levgrnd']))

    #creating variables for file
    times = ds.createVariable('time', 'S13', ('time',))
    lats = ds.createVariable('lat', 'f4', ('lat',))
    lons = ds.createVariable('lon', 'f4', ('lon',))
    BSW = ds.createVariable('BSW', 'f4', ('init','levgrnd', 'lat', 'lon',))
    BSW.units = 'unitless'
    TG = ds.createVariable('TG', 'f4', ('init', 'time','lat', 'lon',))
    TG.units = 'K'
    DZSOI = ds.createVariable('DZSOI', 'f4', ('init','levgrnd', 'lat', 'lon',))
    DZSOI.units = 'm'
    H2OCAN = ds.createVariable('H2OCAN', 'f4', ('init', 'time', 'lat', 'lon',))
    H2OCAN.units = 'mm'
    H2OSFC = ds.createVariable('H2OSFC', 'f4', ('init', 'time', 'lat', 'lon',))
    H2OSFC.units = 'mm'
    QTOPSOIL = ds.createVariable('QTOPSOIL', 'f4', ('init', 'time', 'lat', 'lon',))
    QTOPSOIL.units = 'mm/s'
    SUCSAT = ds.createVariable('SUCSAT', 'f4', ('init','levgrnd', 'lat', 'lon',))
    SUCSAT.units = 'mm'
    TWS = ds.createVariable('TWS', 'f4', ('init', 'time', 'lat', 'lon',))
    TWS.units = 'mm'
    WATSAT = ds.createVariable('WATSAT', 'f4', ('init', 'levgrnd', 'lat', 'lon',))
    WATSAT.units = 'mm3/mm3'
    ZWT = ds.createVariable('ZWT', 'f4', ('init', 'time', 'lat', 'lon',))
    ZWT.units = 'm'
    SOILWATER_10CM = ds.createVariable('SOILWATER_10CM', 'f4', ('init', 'time', 'lat', 'lon',))
    SOILWATER_10CM.units = 'kg/m2'
    ZSOI = ds.createVariable('ZSOI', 'f4', ('init','levgrnd', 'lat', 'lon',))
    ZSOI.units = 'm'

    #appending times from string array into the file
    for x in range(len(Time)):
        times[x] = Time[x]

    #adding data from the large file into the smaller file
    lats[:] = ds_elm['lat']
    lons[:] = ds_elm['lon']
    BSW[:] = ds_elm['BSW']
    TG[:] = ds_elm['TG']
    DZSOI[:] = ds_elm['DZSOI']
    H2OCAN[:] = ds_elm['H2OCAN']
    H2OSFC[:] = ds_elm['H2OSFC']
    QTOPSOIL[:] = ds_elm['QTOPSOIL']
    SUCSAT[:] = ds_elm['SUCSAT']
    TWS[:] = ds_elm['TWS']
    WATSAT[:] = ds_elm['WATSAT']
    ZWT[:] = ds_elm['ZWT']
    SOILWATER_10CM[:] = ds_elm['SOILWATER_10CM']
    ZSOI[:] = ds_elm['ZSOI']
    
    ds.close()


    # ## Creating H1 file for accumulated precip

    # In[6]:


    # Copy datasets to base files for ease of use with jupyter
    ds_eam_i = ds_eam_full

    ## Average dataset across different init times
    ds_eam = ds_eam_i.mean( dim=("init"), skipna=True, keep_attrs=True )

    #Converting units
    ds_eam['PRECT'] = ds_eam['PRECT'] * 8.64e7
    ds_eam['PRECT'].attrs["units"]="mm/d"

    #setting start and end dates of file to match the times of the slice
    start = datetime(1996,1,15,12)
    end = datetime(1996,1,22,23)

    #creating a time array
    time = []
    for dt in datetime_range(start, end, {'hours':3}):
        time.append(dt)
    #converting the time array into a string for the file
    Time = []
    for x in range(len(time)):
        Time.append(time[x].strftime('%Y-%m-%d-%H'))

    #creating file name and opening file
    file = title + '_eam_h1_slice.nc'
    ds = nc.Dataset(file,'w', format='NETCDF4')

    #creating the files dimensions
    time = ds.createDimension('time', len(ds_eam['time']))
    lat = ds.createDimension('lat', len(ds_eam['lat']))
    lon = ds.createDimension('lon', len(ds_eam['lon']))

    #creating variables for file
    times = ds.createVariable('time', 'S13', ('time',))
    lats = ds.createVariable('lat', 'f4', ('lat',))
    lons = ds.createVariable('lon', 'f4', ('lon',))
    PRECT = ds.createVariable('PRECT', 'f4', ('time', 'lat', 'lon',))
    PRECT.units = 'mm/d'

    #adding the data to the new file
    for x in range(len(Time)):
        times[x] = Time[x]
    lats[:] = ds_eam['lat']
    lons[:] = ds_eam['lon']
    PRECT[:] = ds_eam['PRECT']
    
    ds.close()


    # ## Creating EAM file averaged over initialization times

    # In[7]:


    # Copy datasets to base files for ease of use with jupyter
    ds_eam_i = ds_eam_full
    # Average datasets across different init times
    ds_eam = ds_eam_i.mean( dim=("init"), skipna=True, keep_attrs=True )

    #Converting units
    ds_eam['PRECT'] = ds_eam['PRECT'] * 8.64e7
    ds_eam['PRECT'].attrs["units"]="mm/d"
    ds_eam['TREFHT'] = ds_eam['TREFHT'] - 273.15
    ds_eam['TREFHT'].attrs["units"]="degC"

    #creating an array with the time values for the file
    start = datetime(1996,1,15,12)
    end = datetime(1996,1,22,23)

    #creating a time array
    time = []
    for dt in datetime_range(start, end, {'hours':3}):
        time.append(dt)
    
    #converting the time array into a string for the file
    Time = []
    for x in range(len(time)):
        Time.append(time[x].strftime('%Y-%m-%d-%H'))

    #creating file name and opening file
    file = title + '_eam_slice.nc'
    ds = nc.Dataset(file,'w', format='NETCDF4')

    #creating the files dimensions
    time = ds.createDimension('time', len(ds_eam['time']))
    lat = ds.createDimension('lat', len(ds_eam['lat']))
    lon = ds.createDimension('lon', len(ds_eam['lon']))

    #creating variables for file
    times = ds.createVariable('time', 'S13', ('time',))
    lats = ds.createVariable('lat', 'f4', ('lat',))
    lons = ds.createVariable('lon', 'f4', ('lon',))
    PRECT = ds.createVariable('PRECT', 'f4', ('time', 'lat', 'lon',))
    PRECT.units = 'mm/d'
    TREFHT = ds.createVariable('TREFHT', 'f4', ('time', 'lat', 'lon',))
    TREFHT.units = 'C'
    RHREFHT = ds.createVariable('RHREFHT', 'f4', ('time', 'lat', 'lon',))
    RHREFHT.units = 'fraction'
    WSPDSRFAV = ds.createVariable('WSPDSRFAV', 'f4', ('time', 'lat', 'lon',))
    WSPDSRFAV.units = 'm/s'
    CAPE = ds.createVariable('CAPE', 'f4', ('time', 'lat', 'lon',))
    CAPE.units = 'J/Kg'
    PSL = ds.createVariable('PSL', 'f4', ('time', 'lat', 'lon',))
    PSL.units = 'Pa'
    TMQ = ds.createVariable('TMQ', 'f4', ('time', 'lat', 'lon',))
    TMQ.units = 'kg/m^2'
    OMEGA850 = ds.createVariable('OMEGA850', 'f4', ('time', 'lat', 'lon',))
    OMEGA850.units = 'Pa/s'
    OMEGA500 = ds.createVariable('OMEGA500', 'f4', ('time', 'lat', 'lon',))
    OMEGA500.units = 'Pa/s'

    #adding the data to the new file
    for x in range(len(Time)):
        times[x] = Time[x]
    lats[:] = ds_eam['lat']
    lons[:] = ds_eam['lon']
    PRECT[:] = ds_eam['PRECT']
    TREFHT[:] = ds_eam['TREFHT']
    RHREFHT[:] = ds_eam['RHREFHT']
    WSPDSRFAV[:] = ds_eam['WSPDSRFAV']
    CAPE[:] = ds_eam['CAPE']
    PSL[:] = ds_eam['PSL']
    TMQ[:] = ds_eam['TMQ']
    OMEGA850[:] = ds_eam['OMEGA850']
    OMEGA500[:] = ds_eam['OMEGA500']
    
    ds.close()


    # ## Creating EAM file NOT averaged over initialization times

    # In[8]:


    # Copy datasets to base files for ease of use with jupyter
    ds_eam = ds_eam_full

    #setting start and end dates of file to match the times of the slice
    start = datetime(1996,1,15,12)
    end = datetime(1996,1,22,23)

    #Converting units
    ds_eam['PRECT'] = ds_eam['PRECT'] * 8.64e7
    ds_eam['PRECT'].attrs["units"]="mm/d"
    ds_eam['TREFHT'] = ds_eam['TREFHT'] - 273.15
    ds_eam['TREFHT'].attrs["units"]="degC"

    #creating a time array
    time = []
    for dt in datetime_range(start, end, {'hours':3}):
        time.append(dt)
    Time = []
    for x in range(len(time)):
        Time.append(time[x].strftime('%Y-%m-%d-%H'))

    #creating file name and opening file
    file = title + '_init_eam_slice.nc'
    ds = nc.Dataset(file,'w', format='NETCDF4')
    
    #creating the files dimensions
    time = ds.createDimension('time', len(ds_eam['time']))
    init = ds.createDimension('init', len(ds_eam['init']))
    lat = ds.createDimension('lat', len(ds_eam['lat']))
    lon = ds.createDimension('lon', len(ds_eam['lon']))

    #creating the variables for file
    times = ds.createVariable('time', 'S13', ('time',))
    lats = ds.createVariable('lat', 'f4', ('lat',))
    lons = ds.createVariable('lon', 'f4', ('lon',))
    PRECT = ds.createVariable('PRECT', 'f4', ('init', 'time', 'lat', 'lon',))
    PRECT.units = 'mm/d'
    TREFHT = ds.createVariable('TREFHT', 'f4', ('init', 'time', 'lat', 'lon',))
    TREFHT.units = 'C'
    RHREFHT = ds.createVariable('RHREFHT', 'f4', ('init', 'time', 'lat', 'lon',))
    RHREFHT.units = 'fraction'
    WSPDSRFAV = ds.createVariable('WSPDSRFAV', 'f4', ('init', 'time', 'lat', 'lon',))
    WSPDSRFAV.units = 'm/s'
    CAPE = ds.createVariable('CAPE', 'f4', ('init', 'time', 'lat', 'lon',))
    CAPE.units = 'J/Kg'
    PSL = ds.createVariable('PSL', 'f4', ('init', 'time', 'lat', 'lon',))
    PSL.units = 'Pa'
    TMQ = ds.createVariable('TMQ', 'f4', ('init', 'time', 'lat', 'lon',))
    TMQ.units = 'kg/m^2'
    OMEGA850 = ds.createVariable('OMEGA850', 'f4', ('init', 'time', 'lat', 'lon',))
    OMEGA850.units = 'Pa/s'
    OMEGA500 = ds.createVariable('OMEGA500', 'f4', ('init', 'time', 'lat', 'lon',))
    OMEGA500.units = 'Pa/s'

    #adding data to file
    for x in range(len(Time)):
        times[x] = Time[x]
    lats[:] = ds_eam['lat']
    lons[:] = ds_eam['lon']
    PRECT[:] = ds_eam['PRECT']
    TREFHT[:] = ds_eam['TREFHT']
    RHREFHT[:] = ds_eam['RHREFHT']
    WSPDSRFAV[:] = ds_eam['WSPDSRFAV']
    CAPE[:] = ds_eam['CAPE']
    PSL[:] = ds_eam['PSL']
    TMQ[:] = ds_eam['TMQ']
    OMEGA850[:] = ds_eam['OMEGA850']
    OMEGA500[:] = ds_eam['OMEGA500']
    
    ds.close()


    # ## Creating ELM file averaged over initialization times

    # In[9]:


    import warnings
    warnings.filterwarnings('ignore')

    # Copy datasets to base files for ease of use with jupyter
    ds_elm_i = ds_elm_full
    # Average datasets across different init times
    ds_elm = ds_elm_i.mean( dim=("init"), skipna=True, keep_attrs=True )

    #converting units 
    ds_elm['QRUNOFF'] = ds_elm['QRUNOFF'] * 86400
    ds_elm['QRUNOFF'].attrs["units"]="mm/d"
    ds_elm['QSNOMELT'] = ds_elm['QSNOMELT'] * 86400
    ds_elm['QSNOMELT'].attrs["units"]="mm/d"
    ds_elm['QOVER'] = ds_elm['QOVER'] * 86400
    ds_elm['QOVER'].attrs["units"]="mm/d"
    ds_elm['QFLOOD'] = ds_elm['QFLOOD'] * 86400
    ds_elm['QFLOOD'].attrs["units"]="mm/d"
    ds_elm['QH2OSFC'] = ds_elm['QH2OSFC'] * 86400
    ds_elm['QH2OSFC'].attrs["units"]="mm/d"
    ds_elm['SNOWDP'] = ds_elm['SNOWDP'] * 1000
    ds_elm['SNOWDP'].attrs["units"]="mm"
    ds_elm['TSA'] = ds_elm['TSA'] - 273.15

    #setting start and end dates of file to match the times of the slice
    start = datetime(1996,1,15,12)
    end = datetime(1996,1,22,23)

    #creating a time array
    time = []
    for dt in datetime_range(start, end, {'hours':6}):
        time.append(dt)
    #converting the time array into a string for the file
    Time = []
    for x in range(len(time)):
        Time.append(time[x].strftime('%Y-%m-%d-%H'))

    #creating file name and opening file
    file = title + '_elm_slice.nc'
    ds = nc.Dataset(file,'w', format='NETCDF4')

    #creating dimensions of file
    time = ds.createDimension('time', len(ds_elm['time']))
    lat = ds.createDimension('lat', len(ds_elm['lat']))
    lon = ds.createDimension('lon', len(ds_elm['lon']))

    #creating variables for the file
    times = ds.createVariable('time', 'S13', ('time',))
    lats = ds.createVariable('lat', 'f4', ('lat',))
    lons = ds.createVariable('lon', 'f4', ('lon',))
    QRUNOFF = ds.createVariable('QRUNOFF', 'f4', ('time', 'lat', 'lon',))
    QRUNOFF.units = 'mm/d'
    QOVER = ds.createVariable('QOVER', 'f4', ('time', 'lat', 'lon',))
    QOVER.units = 'mm/d'
    QFLOOD = ds.createVariable('QFLOOD', 'f4', ('time', 'lat', 'lon',))
    QFLOOD.units = 'mm/d'
    QH2OSFC = ds.createVariable('QH2OSFC', 'f4', ('time', 'lat', 'lon',))
    QH2OSFC.units = 'mm/d'
    QSNOMELT = ds.createVariable('QSNOMELT', 'f4', ('time', 'lat', 'lon',))
    QSNOMELT.units = 'mm/d'
    H2OSNO = ds.createVariable('H2OSNO', 'f4', ('time', 'lat', 'lon',))
    H2OSNO.units = 'mm'
    LWup = ds.createVariable('LWup', 'f4', ('time', 'lat', 'lon',))
    LWup.units = 'W/m^2'
    LWdown = ds.createVariable('LWdown', 'f4', ('time', 'lat', 'lon',))
    LWdown.units = 'W/m^2'
    SWup = ds.createVariable('SWup', 'f4', ('time', 'lat', 'lon',))
    SWup.units = 'W/m^2'
    SWdown = ds.createVariable('SWdown', 'f4', ('time', 'lat', 'lon',))
    SWdown.units = 'W/m^2'
    Qh = ds.createVariable('Qh', 'f4', ('time', 'lat', 'lon',))
    Qh.units = 'W/m^2'
    Qle = ds.createVariable('Qle', 'f4', ('time', 'lat', 'lon',))
    Qle.units = 'W/m^2'
    H2OSFC = ds.createVariable('H2OSFC', 'f4', ('time', 'lat', 'lon',))
    H2OSFC.units = 'mm'
    Rnet = ds.createVariable('Rnet', 'f4', ('time', 'lat', 'lon',))
    Rnet.units = 'W/m^2'
    TSA = ds.createVariable('TSA', 'f4', ('time', 'lat', 'lon',))
    TSA.units = 'C'
    TWS = ds.createVariable('TWS', 'f4', ('time', 'lat', 'lon',))
    TWS.units = 'mm'
    H2OSNO_TOP = ds.createVariable('H2OSNO_TOP', 'f4', ('time', 'lat', 'lon',))
    H2OSNO_TOP.units = 'kg/m^2'
    SNOWDP = ds.createVariable('SNOWDP', 'f4', ('time', 'lat', 'lon',))
    SNOWDP.units = 'mm'
    FGR = ds.createVariable('FGR', 'f4', ('time', 'lat', 'lon',))
    FGR.units = 'W/m^2'
    FSM = ds.createVariable('FSM', 'f4', ('time', 'lat', 'lon',))
    FSM.units = 'W/m^2'
    SNOWLIQ = ds.createVariable('SNOWLIQ', 'f4', ('time', 'lat', 'lon',))
    SNOWLIQ.units = 'kg/m^2'

    #adding data to the file
    for x in range(len(Time)):
        times[x] = Time[x]
    lats[:] = ds_elm['lat']
    lons[:] = ds_elm['lon']
    QRUNOFF[:] = ds_elm['QRUNOFF']
    QOVER[:] = ds_elm['QOVER']
    QSNOMELT[:] = ds_elm['QSNOMELT']
    QH2OSFC[:] = ds_elm['QH2OSFC']
    QFLOOD[:] = ds_elm['QFLOOD']
    H2OSNO[:] = ds_elm['H2OSNO']
    LWup[:] = ds_elm['LWup']
    LWdown[:] = ds_elm['LWdown']
    SWup[:] = ds_elm['SWup']
    SWdown[:] = ds_elm['SWdown']
    Qh[:] = ds_elm['Qh']
    Qle[:] = ds_elm['Qle']
    H2OSFC[:] = ds_elm['H2OSFC']
    Rnet[:] = ds_elm['Rnet']
    TSA[:] = ds_elm['TSA']
    TWS[:] = ds_elm['TWS']
    FGR[:] = ds_elm['FGR']
    FSM[:] = ds_elm['FSM']
    H2OSNO_TOP[:] = ds_elm['H2OSNO_TOP']
    SNOWLIQ[:] = ds_elm['SNOWLIQ']
    SNOWDP[:] = ds_elm['SNOWDP']
    
    ds.close()


    # ## Creating ELM file NOT averaged over initialization times

    # In[10]:


    # Copy datasets to base files for ease of use with jupyter
    ds_elm = ds_elm_full

    #setting start and end dates of file to match the times of the slice
    start = datetime(1996,1,15,12)
    end = datetime(1996,1,22,23)

    #converting units 
    ds_elm['QRUNOFF'] = ds_elm['QRUNOFF'] * 86400
    ds_elm['QRUNOFF'].attrs["units"]="mm/d"
    ds_elm['QSNOMELT'] = ds_elm['QSNOMELT'] * 86400
    ds_elm['QSNOMELT'].attrs["units"]="mm/d"
    ds_elm['QOVER'] = ds_elm['QOVER'] * 86400
    ds_elm['QOVER'].attrs["units"]="mm/d"
    ds_elm['QFLOOD'] = ds_elm['QFLOOD'] * 86400
    ds_elm['QFLOOD'].attrs["units"]="mm/d"
    ds_elm['QH2OSFC'] = ds_elm['QH2OSFC'] * 86400
    ds_elm['QH2OSFC'].attrs["units"]="mm/d"
    ds_elm['SNOWDP'] = ds_elm['SNOWDP'] * 1000
    ds_elm['SNOWDP'].attrs["units"]="mm"
    ds_elm['TSA'] = ds_elm['TSA'] - 273.15

    #creating a time array
    time = []
    for dt in datetime_range(start, end, {'hours':6}):
        time.append(dt)
    Time = []
    for x in range(len(time)):
        Time.append(time[x].strftime('%Y-%m-%d-%H'))

    #creating file name and opening file
    file = title + '_init_elm_slice.nc'
    ds = nc.Dataset(file,'w', format='NETCDF4')

    #creating dimensions of file
    time = ds.createDimension('time', len(ds_elm['time']))
    init = ds.createDimension('init', len(ds_elm['init']))
    lat = ds.createDimension('lat', len(ds_elm['lat']))
    lon = ds.createDimension('lon', len(ds_elm['lon']))

    #Creating variables for file
    times = ds.createVariable('time', 'S13', ('time',))
    lats = ds.createVariable('lat', 'f4', ('lat',))
    lons = ds.createVariable('lon', 'f4', ('lon',))
    QRUNOFF = ds.createVariable('QRUNOFF', 'f4', ('init','time', 'lat', 'lon',))
    QRUNOFF.units = 'mm/d'
    QOVER = ds.createVariable('QOVER', 'f4', ('init','time', 'lat', 'lon',))
    QOVER.units = 'mm/d'
    QFLOOD = ds.createVariable('QFLOOD', 'f4', ('init','time', 'lat', 'lon',))
    QFLOOD.units = 'mm/d'
    QH2OSFC = ds.createVariable('QH2OSFC', 'f4', ('init','time', 'lat', 'lon',))
    QH2OSFC.units = 'mm/d'
    QSNOMELT = ds.createVariable('QSNOMELT', 'f4', ('init','time', 'lat', 'lon',))
    QSNOMELT.units = 'mm/d'
    H2OSNO = ds.createVariable('H2OSNO', 'f4', ('init','time', 'lat', 'lon',))
    H2OSNO.units = 'mm'
    LWup = ds.createVariable('LWup', 'f4', ('init','time', 'lat', 'lon',))
    LWup.units = 'W/m^2'
    LWdown = ds.createVariable('LWdown', 'f4', ('init','time', 'lat', 'lon',))
    LWdown.units = 'W/m^2'
    SWup = ds.createVariable('SWup', 'f4', ('init','time', 'lat', 'lon',))
    SWup.units = 'W/m^2'
    SWdown = ds.createVariable('SWdown', 'f4', ('init','time', 'lat', 'lon',))
    SWdown.units = 'W/m^2'
    Qh = ds.createVariable('Qh', 'f4', ('init','time', 'lat', 'lon',))
    Qh.units = 'W/m^2'
    Qle = ds.createVariable('Qle', 'f4', ('init','time', 'lat', 'lon',))
    Qle.units = 'W/m^2'
    H2OSFC = ds.createVariable('H2OSFC', 'f4', ('init','time', 'lat', 'lon',))
    H2OSFC.units = 'mm'
    Rnet = ds.createVariable('Rnet', 'f4', ('init','time', 'lat', 'lon',))
    Rnet.units = 'W/m^2'
    TSA = ds.createVariable('TSA', 'f4', ('init','time', 'lat', 'lon',))
    TSA.units = 'C'
    TWS = ds.createVariable('TWS', 'f4', ('init','time', 'lat', 'lon',))
    TWS.units = 'mm'
    H2OSNO_TOP = ds.createVariable('H2OSNO_TOP', 'f4', ('init','time', 'lat', 'lon',))
    H2OSNO_TOP.units = 'kg/m^2'
    SNOWDP = ds.createVariable('SNOWDP', 'f4', ('init','time', 'lat', 'lon',))
    SNOWDP.units = 'mm'
    FGR = ds.createVariable('FGR', 'f4', ('init','time', 'lat', 'lon',))
    FGR.units = 'W/m^2'
    FSM = ds.createVariable('FSM', 'f4', ('init','time', 'lat', 'lon',))
    FSM.units = 'W/m^2'
    SNOWLIQ = ds.createVariable('SNOWLIQ', 'f4', ('init','time', 'lat', 'lon',))
    SNOWLIQ.units = 'kg/m^2'

    #adding data to file
    for x in range(len(Time)):
        times[x] = Time[x]
    lats[:] = ds_elm['lat']
    lons[:] = ds_elm['lon']
    QRUNOFF[:] = ds_elm['QRUNOFF']
    QOVER[:] = ds_elm['QOVER']
    QSNOMELT[:] = ds_elm['QSNOMELT']
    QH2OSFC[:] = ds_elm['QH2OSFC']
    QFLOOD[:] = ds_elm['QFLOOD']
    H2OSNO[:] = ds_elm['H2OSNO']
    LWup[:] = ds_elm['LWup']
    LWdown[:] = ds_elm['LWdown']
    SWup[:] = ds_elm['SWup']
    SWdown[:] = ds_elm['SWdown']
    Qh[:] = ds_elm['Qh']
    Qle[:] = ds_elm['Qle']
    H2OSFC[:] = ds_elm['H2OSFC']
    Rnet[:] = ds_elm['Rnet']
    TSA[:] = ds_elm['TSA']
    TWS[:] = ds_elm['TWS']
    FGR[:] = ds_elm['FGR']
    FSM[:] = ds_elm['FSM']
    H2OSNO_TOP[:] = ds_elm['H2OSNO_TOP']
    SNOWLIQ[:] = ds_elm['SNOWLIQ']
    SNOWDP[:] = ds_elm['SNOWDP']
    
    ds.close()

    # ## Creates MOSART file (averaged over initialization time)

    #averaging file over initialization time
    ds_mos_i = ds_mos_full
    ds_mos = ds_mos_i.mean( dim=("init"), skipna=True, keep_attrs=True )

    #setting start and end dates of file to match the times of the slice
    start = datetime(1996,1,16,12)
    end = datetime(1996,1,24,13)

    #creating a time array
    time = []
    for dt in datetime_range(start, end, {'hours':12}):
        time.append(dt)
    Time = []
    for x in range(len(time)):
        Time.append(time[x].strftime('%Y-%m-%d-%H'))

    #creating file name and opening it
    file = title + '_mos_slice.nc'
    ds = nc.Dataset(file,'w', format='NETCDF4')

    #creating dimensions of file
    time = ds.createDimension('time', len(ds_mos['time']))
    lat = ds.createDimension('lat', len(ds_mos['lat']))
    lon = ds.createDimension('lon', len(ds_mos['lon']))

    #Creating variables for file
    times = ds.createVariable('time', 'S13', ('time',))
    lats = ds.createVariable('lat', 'f4', ('lat',))
    lons = ds.createVariable('lon', 'f4', ('lon',))
    RIVER_DISCHARGE_OVER_LAND_LIQ = ds.createVariable('RIVER_DISCHARGE_OVER_LAND_LIQ', 'f4', ('time', 'lat', 'lon',))
    RIVER_DISCHARGE_OVER_LAND_LIQ.units = 'm^3/s'
    RIVER_DISCHARGE_TO_OCEAN_LIQ = ds.createVariable('RIVER_DISCHARGE_TO_OCEAN_LIQ', 'f4', ('time', 'lat', 'lon',))
    RIVER_DISCHARGE_TO_OCEAN_LIQ.units = 'C'

    #adding data to file
    for x in range(len(Time)):
        times[x] = Time[x]
    lats[:] = ds_mos['lat']
    lons[:] = ds_mos['lon']
    RIVER_DISCHARGE_OVER_LAND_LIQ[:] = ds_mos['RIVER_DISCHARGE_OVER_LAND_LIQ']
    RIVER_DISCHARGE_TO_OCEAN_LIQ[:] = ds_mos['RIVER_DISCHARGE_TO_OCEAN_LIQ']

    ds.close()


    # ## Creates MOSART file (NOT averaged over initialization time)

    # In[12]:


    #Copy dataset to base file (which is not averaged over init)
    ds_mos = ds_mos_full

    #creating files where the initialization times are not averaged but are instead a fourth dimension
    start = datetime(1996,1,16,12)
    end = datetime(1996,1,24,13)

    #creating a time array
    time = []
    for dt in datetime_range(start, end, {'hours':12}):
        time.append(dt)
    Time = []
    for x in range(len(time)):
        Time.append(time[x].strftime('%Y-%m-%d-%H'))

    #creating file name and opening it
    file = title + '_init_mos_slice.nc'
    ds = nc.Dataset(file,'w', format='NETCDF4')

    #creating dimensions of file
    time = ds.createDimension('time', len(ds_mos['time']))
    init = ds.createDimension('init', len(ds_mos['init']))
    lat = ds.createDimension('lat', len(ds_mos['lat']))
    lon = ds.createDimension('lon', len(ds_mos['lon']))

    #creating variables for file
    times = ds.createVariable('time', 'S13', ('time',))
    lats = ds.createVariable('lat', 'f4', ('lat',))
    lons = ds.createVariable('lon', 'f4', ('lon',))
    RIVER_DISCHARGE_OVER_LAND_LIQ = ds.createVariable('RIVER_DISCHARGE_OVER_LAND_LIQ', 'f4', ('init','time', 'lat', 'lon',))
    RIVER_DISCHARGE_OVER_LAND_LIQ.units = 'm^3/s'
    RIVER_DISCHARGE_TO_OCEAN_LIQ = ds.createVariable('RIVER_DISCHARGE_TO_OCEAN_LIQ', 'f4', ('init','time', 'lat', 'lon',))
    RIVER_DISCHARGE_TO_OCEAN_LIQ.units = 'C'

    #adding data to file
    for x in range(len(Time)):
        times[x] = Time[x]
    lats[:] = ds_mos['lat']
    lons[:] = ds_mos['lon']
    RIVER_DISCHARGE_OVER_LAND_LIQ[:] = ds_mos['RIVER_DISCHARGE_OVER_LAND_LIQ']
    RIVER_DISCHARGE_TO_OCEAN_LIQ[:] = ds_mos['RIVER_DISCHARGE_TO_OCEAN_LIQ']

    ds.close()

    # Increment counter
    counter += 1