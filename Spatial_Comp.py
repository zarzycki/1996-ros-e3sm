#!/usr/bin/env python
# coding: utf-8

# # Creates Figure 2, Figure 8, Figure 9, Figure 10, Figure 16 

# importing necessary libraries
import numpy as np
import netCDF4
import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xesmf as xe
import datetime
from datetime import timedelta
import shapefile as shp
import geopandas as gpd

#set to true to use shapefile to average over SRB set to false to use rough basin average
SRB_AVG = True

#ELM and MOSART variables are time averaged, so to shift the timeseries plots to fall at the middle 
#of the time period rather than the end of the time period set equal to True rather than False
shift_ELM = True
shift_MOSART = True

#creating letter array to label subplots
letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)']

#loading ensemble files
ds_eam_Control = xr.open_dataset('ds_eam_Control_ens.nc')
ds_elm_Control = xr.open_dataset('ds_elm_Control_ens.nc')
ds_mos_Control = xr.open_dataset('ds_mos_Control_ens.nc')

ds_eam_Plus_1K = xr.open_dataset('ds_eam_Plus_1K_ens.nc')
ds_elm_Plus_1K = xr.open_dataset('ds_elm_Plus_1K_ens.nc')
ds_mos_Plus_1K = xr.open_dataset('ds_mos_Plus_1K_ens.nc')

ds_eam_Plus_2K = xr.open_dataset('ds_eam_Plus_2K_ens.nc')
ds_elm_Plus_2K = xr.open_dataset('ds_elm_Plus_2K_ens.nc')
ds_mos_Plus_2K = xr.open_dataset('ds_mos_Plus_2K_ens.nc')

ds_eam_Plus_3K = xr.open_dataset('ds_eam_Plus_3K_ens.nc')
ds_elm_Plus_3K = xr.open_dataset('ds_elm_Plus_3K_ens.nc')
ds_mos_Plus_3K = xr.open_dataset('ds_mos_Plus_3K_ens.nc')

ds_eam_Plus_4K = xr.open_dataset('ds_eam_Plus_4K_ens.nc')
ds_elm_Plus_4K = xr.open_dataset('ds_elm_Plus_4K_ens.nc')
ds_mos_Plus_4K = xr.open_dataset('ds_mos_Plus_4K_ens.nc')

ds_eam_Pre_Ind = xr.open_dataset('ds_eam_Pre_Ind_ens.nc')
ds_elm_Pre_Ind = xr.open_dataset('ds_elm_Pre_Ind_ens.nc')
ds_mos_Pre_Ind = xr.open_dataset('ds_mos_Pre_Ind_ens.nc')

#loading ensemble files
ds_eam_Controli = xr.open_dataset('ds_eam_Controli_ens.nc')
ds_elm_Controli = xr.open_dataset('ds_elm_Controli_ens.nc')
ds_mos_Controli = xr.open_dataset('ds_mos_Controli_ens.nc')

ds_eam_Plus_1Ki = xr.open_dataset('ds_eam_Plus_1Ki_ens.nc')
ds_elm_Plus_1Ki = xr.open_dataset('ds_elm_Plus_1Ki_ens.nc')
ds_mos_Plus_1Ki = xr.open_dataset('ds_mos_Plus_1Ki_ens.nc')

ds_eam_Plus_2Ki = xr.open_dataset('ds_eam_Plus_2Ki_ens.nc')
ds_elm_Plus_2Ki = xr.open_dataset('ds_elm_Plus_2Ki_ens.nc')
ds_mos_Plus_2Ki = xr.open_dataset('ds_mos_Plus_2Ki_ens.nc')

ds_eam_Plus_3Ki = xr.open_dataset('ds_eam_Plus_3Ki_ens.nc')
ds_elm_Plus_3Ki = xr.open_dataset('ds_elm_Plus_3Ki_ens.nc')
ds_mos_Plus_3Ki = xr.open_dataset('ds_mos_Plus_3Ki_ens.nc')

ds_eam_Plus_4Ki = xr.open_dataset('ds_eam_Plus_4Ki_ens.nc')
ds_elm_Plus_4Ki = xr.open_dataset('ds_elm_Plus_4Ki_ens.nc')
ds_mos_Plus_4Ki = xr.open_dataset('ds_mos_Plus_4Ki_ens.nc')

ds_eam_Pre_Indi = xr.open_dataset('ds_eam_Pre_Indi_ens.nc')
ds_elm_Pre_Indi = xr.open_dataset('ds_elm_Pre_Indi_ens.nc')
ds_mos_Pre_Indi = xr.open_dataset('ds_mos_Pre_Indi_ens.nc')

#function that creates an array of time values given a start date, end date, and timestep
def datetime_range(start, end, delta):
    current = start
    if not isinstance(delta, timedelta):
        delta = timedelta(**delta)
    while current < end:
        yield current
        current += delta

#creating time arrays with datetime objects
time_m = ds_eam_Controli['time'].values
time_m = time_m.tolist()
Time_eam = []
for i in range(len(ds_eam_Controli['time'])):
    Time_eam.append(datetime.datetime.strptime(time_m[i], '%Y-%m-%d-%H'))

if shift_ELM == False:
    time_m2 = ds_elm_Controli['time'].values
    time_m2 = time_m2.tolist()
    Time_elm = []
    for i in range(len(ds_elm_Controli['time'])):
        Time_elm.append(datetime.datetime.strptime(time_m2[i], '%Y-%m-%d-%H'))
if shift_MOSART == False:
    time_m3 = ds_mos_Controli['time'].values
    time_m3 = time_m3.tolist()
    Time_mos = []
    for i in range(len(ds_mos_Controli['time'])):
        Time_mos.append(datetime.datetime.strptime(time_m3[i], '%Y-%m-%d-%H'))

#shifting elm and mos times if set to true
if shift_ELM == True:
    start_elm = datetime.datetime(1996,1,15,9)
    end_elm = datetime.datetime(1996,1,22,21)
    Time_elm = []
    for dt in datetime_range(start_elm, end_elm, {'hours':6}):
        Time_elm.append(dt)
    
if shift_MOSART == True:
    start_mos = datetime.datetime(1996,1,16,0)
    end_mos = datetime.datetime(1996,1,24,12)
    Time_mos = []
    for dt in datetime_range(start_mos, end_mos, {'hours':12}):
        Time_mos.append(dt)
        
#Changing the time data in the files if the shift is set to True
ELM_runs = [ds_elm_Pre_Indi, ds_elm_Controli, ds_elm_Plus_1Ki, ds_elm_Plus_2Ki, ds_elm_Plus_3Ki,ds_elm_Plus_4Ki]
MOS_runs = [ds_mos_Pre_Indi, ds_mos_Controli, ds_mos_Plus_1Ki, ds_mos_Plus_2Ki, ds_mos_Plus_3Ki,ds_mos_Plus_4Ki]
if shift_ELM == True:
    for x in ELM_runs:
        x['time'] = Time_elm
if shift_MOSART == True:
    for x in MOS_runs:
        x['time'] = Time_mos

#regridding mosart data to match high res control
ds_out = xr.Dataset({'lat': (['lat'], np.arange(35, 45.1, 0.125)),
                     'lon': (['lon'], np.arange(-85, -69.875, .125)),
                    })

regridder_mos = xe.Regridder(ds_mos_Controli, ds_out, 'bilinear')
regridder_mos2 = xe.Regridder(ds_mos_Control, ds_out, 'bilinear')
ds_mos_Controli = regridder_mos(ds_mos_Controli)
ds_mos_Plus_1Ki = regridder_mos(ds_mos_Plus_1Ki)
ds_mos_Plus_2Ki = regridder_mos(ds_mos_Plus_2Ki)
ds_mos_Plus_3Ki = regridder_mos(ds_mos_Plus_3Ki)
ds_mos_Plus_4Ki = regridder_mos(ds_mos_Plus_4Ki)
ds_mos_Pre_Indi = regridder_mos(ds_mos_Pre_Indi)

ds_mos_Control = regridder_mos2(ds_mos_Control)
ds_mos_Plus_1K = regridder_mos2(ds_mos_Plus_1K)
ds_mos_Plus_2K = regridder_mos2(ds_mos_Plus_2K)
ds_mos_Plus_3K = regridder_mos2(ds_mos_Plus_3K)
ds_mos_Plus_4K = regridder_mos2(ds_mos_Plus_4K)
ds_mos_Pre_Ind = regridder_mos2(ds_mos_Pre_Ind)


#this cell creates the necessary functions to average a variable over the actual SRB rather than the square above

from rasterio import features
from affine import Affine

def transform_from_latlon(lat, lon):
    """ input 1D array of lat / lon and output an Affine transformation
    """
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    trans = Affine.translation(lon[0], lat[0])
    scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
    return trans * scale

def rasterize(shapes, coords, latitude='lat', longitude='lon',
              fill=np.nan, **kwargs):
    """Rasterize a list of (geometry, fill_value) tuples onto the given
    xray coordinates. This only works for 1d latitude and longitude
    arrays.

    usage:
    -----
    1. read shapefile to geopandas.GeoDataFrame
          `states = gpd.read_file(shp_dir+shp_file)`
    2. encode the different shapefiles that capture those lat-lons as different
        numbers i.e. 0.0, 1.0 ... and otherwise np.nan
          `shapes = (zip(states.geometry, range(len(states))))`
    3. Assign this to a new coord in your original xarray.DataArray
          `ds['states'] = rasterize(shapes, ds.coords, longitude='X', latitude='Y')`

    arguments:
    ---------
    : **kwargs (dict): passed to `rasterio.rasterize` function

    attrs:
    -----
    :transform (affine.Affine): how to translate from latlon to ...?
    :raster (numpy.ndarray): use rasterio.features.rasterize fill the values
      outside the .shp file with np.nan
    :spatial_coords (dict): dictionary of {"X":xr.DataArray, "Y":xr.DataArray()}
      with "X", "Y" as keys, and xr.DataArray as values

    returns:
    -------
    :(xr.DataArray): DataArray with `values` of nan for points outside shapefile
      and coords `Y` = latitude, 'X' = longitude.


    """
    transform = transform_from_latlon(coords[latitude], coords[longitude])
    out_shape = (len(coords[latitude]), len(coords[longitude]))
    raster = features.rasterize(shapes, out_shape=out_shape,
                                fill=fill, transform=transform,
                                dtype=float, **kwargs)
    spatial_coords = {latitude: coords[latitude], longitude: coords[longitude]}
    return xr.DataArray(raster, coords=spatial_coords, dims=(latitude, longitude))

def add_shape_coord_from_data_array(xr_da, shp_path, coord_name):
    """ Create a new coord for the xr_da indicating whether or not it 
         is inside the shapefile

        Creates a new coord - "coord_name" which will have integer values
         used to subset xr_da for plotting / analysis/

        Usage:
        -----
        precip_da = add_shape_coord_from_data_array(precip_da, "awash.shp", "awash")
        awash_da = precip_da.where(precip_da.awash==0, other=np.nan) 
    """
    # 1. read in shapefile
    shp_gpd = gpd.read_file(shp_path)
    shp_gpd = shp_gpd.to_crs(epsg=4326)


    # 2. create a list of tuples (shapely.geometry, id)
    #    this allows for many different polygons within a .shp file (e.g. States of US)
    shapes = [(shape, n) for n, shape in enumerate(shp_gpd.geometry)]

    # 3. create a new coord in the xr_da which will be set to the id in `shapes`
    xr_da[coord_name] = rasterize(shapes, xr_da.coords, 
                               longitude='lon', latitude='lat')

    return xr_da

#masking out all data outside of SRB if set to TRUE
if SRB_AVG == True:
    ds_eam_Control = add_shape_coord_from_data_array(ds_eam_Control, "./shapes/srb.shp", 'coord' )
    ds_eam_Control = ds_eam_Control.where(ds_eam_Control.coord==0, other=np.nan)
    ds_elm_Control = add_shape_coord_from_data_array(ds_elm_Control, "./shapes/srb.shp", 'coord' )
    ds_elm_Control = ds_elm_Control.where(ds_elm_Control.coord==0, other=np.nan)
    ds_mos_Control = add_shape_coord_from_data_array(ds_mos_Control, "./shapes/srb.shp", 'coord' )
    ds_mos_Control = ds_mos_Control.where(ds_mos_Control.coord==0, other=np.nan)

    ds_eam_Plus_1K = add_shape_coord_from_data_array(ds_eam_Plus_1K, "./shapes/srb.shp", 'coord' )
    ds_eam_Plus_1K = ds_eam_Plus_1K.where(ds_eam_Plus_1K.coord==0, other=np.nan)
    ds_elm_Plus_1K = add_shape_coord_from_data_array(ds_elm_Plus_1K, "./shapes/srb.shp", 'coord' )
    ds_elm_Plus_1K = ds_elm_Plus_1K.where(ds_elm_Plus_1K.coord==0, other=np.nan)
    ds_mos_Plus_1K = add_shape_coord_from_data_array(ds_mos_Plus_1K, "./shapes/srb.shp", 'coord' )
    ds_mos_Plus_1K = ds_mos_Plus_1K.where(ds_mos_Plus_1K.coord==0, other=np.nan)

    ds_eam_Plus_2K = add_shape_coord_from_data_array(ds_eam_Plus_2K, "./shapes/srb.shp", 'coord' )
    ds_eam_Plus_2K = ds_eam_Plus_2K.where(ds_eam_Plus_2K.coord==0, other=np.nan)
    ds_elm_Plus_2K = add_shape_coord_from_data_array(ds_elm_Plus_2K, "./shapes/srb.shp", 'coord' )
    ds_elm_Plus_2K = ds_elm_Plus_2K.where(ds_elm_Plus_2K.coord==0, other=np.nan)
    ds_mos_Plus_2K = add_shape_coord_from_data_array(ds_mos_Plus_2K, "./shapes/srb.shp", 'coord' )
    ds_mos_Plus_2K = ds_mos_Plus_2K.where(ds_mos_Plus_2K.coord==0, other=np.nan)

    ds_eam_Plus_3K = add_shape_coord_from_data_array(ds_eam_Plus_3K, "./shapes/srb.shp", 'coord' )
    ds_eam_Plus_3K = ds_eam_Plus_3K.where(ds_eam_Plus_3K.coord==0, other=np.nan)
    ds_elm_Plus_3K = add_shape_coord_from_data_array(ds_elm_Plus_3K, "./shapes/srb.shp", 'coord' )
    ds_elm_Plus_3K = ds_elm_Plus_3K.where(ds_elm_Plus_3K.coord==0, other=np.nan)
    ds_mos_Plus_3K = add_shape_coord_from_data_array(ds_mos_Plus_3K, "./shapes/srb.shp", 'coord' )
    ds_mos_Plus_3K = ds_mos_Plus_3K.where(ds_mos_Plus_3K.coord==0, other=np.nan)

    ds_eam_Plus_4K = add_shape_coord_from_data_array(ds_eam_Plus_4K, "./shapes/srb.shp", 'coord' )
    ds_eam_Plus_4K = ds_eam_Plus_4K.where(ds_eam_Plus_4K.coord==0, other=np.nan)
    ds_elm_Plus_4K = add_shape_coord_from_data_array(ds_elm_Plus_4K, "./shapes/srb.shp", 'coord' )
    ds_elm_Plus_4K = ds_elm_Plus_4K.where(ds_elm_Plus_4K.coord==0, other=np.nan)
    ds_mos_Plus_4K = add_shape_coord_from_data_array(ds_mos_Plus_4K, "./shapes/srb.shp", 'coord' )
    ds_mos_Plus_4K = ds_mos_Plus_4K.where(ds_mos_Plus_4K.coord==0, other=np.nan)

    ds_eam_Pre_Ind = add_shape_coord_from_data_array(ds_eam_Pre_Ind, "./shapes/srb.shp", 'coord' )
    ds_eam_Pre_Ind = ds_eam_Pre_Ind.where(ds_eam_Pre_Ind.coord==0, other=np.nan)
    ds_elm_Pre_Ind = add_shape_coord_from_data_array(ds_elm_Pre_Ind, "./shapes/srb.shp", 'coord' )
    ds_elm_Pre_Ind = ds_elm_Pre_Ind.where(ds_elm_Pre_Ind.coord==0, other=np.nan)
    ds_mos_Pre_Ind = add_shape_coord_from_data_array(ds_mos_Pre_Ind, "./shapes/srb.shp", 'coord' )
    ds_mos_Pre_Ind = ds_mos_Pre_Ind.where(ds_mos_Pre_Ind.coord==0, other=np.nan)

    ds_eam_Controli = add_shape_coord_from_data_array(ds_eam_Controli, "./shapes/srb.shp", 'coord' )
    ds_eam_Controli = ds_eam_Controli.where(ds_eam_Controli.coord==0, other=np.nan)
    ds_elm_Controli = add_shape_coord_from_data_array(ds_elm_Controli, "./shapes/srb.shp", 'coord' )
    ds_elm_Controli = ds_elm_Controli.where(ds_elm_Controli.coord==0, other=np.nan)
    ds_mos_Controli = add_shape_coord_from_data_array(ds_mos_Controli, "./shapes/srb.shp", 'coord' )
    ds_mos_Controli = ds_mos_Controli.where(ds_mos_Controli.coord==0, other=np.nan)

    ds_eam_Plus_1Ki = add_shape_coord_from_data_array(ds_eam_Plus_1Ki, "./shapes/srb.shp", 'coord' )
    ds_eam_Plus_1Ki = ds_eam_Plus_1Ki.where(ds_eam_Plus_1Ki.coord==0, other=np.nan)
    ds_elm_Plus_1Ki = add_shape_coord_from_data_array(ds_elm_Plus_1Ki, "./shapes/srb.shp", 'coord' )
    ds_elm_Plus_1Ki = ds_elm_Plus_1Ki.where(ds_elm_Plus_1Ki.coord==0, other=np.nan)
    ds_mos_Plus_1Ki = add_shape_coord_from_data_array(ds_mos_Plus_1Ki, "./shapes/srb.shp", 'coord' )
    ds_mos_Plus_1Ki = ds_mos_Plus_1Ki.where(ds_mos_Plus_1Ki.coord==0, other=np.nan)

    ds_eam_Plus_2Ki = add_shape_coord_from_data_array(ds_eam_Plus_2Ki, "./shapes/srb.shp", 'coord' )
    ds_eam_Plus_2Ki = ds_eam_Plus_2Ki.where(ds_eam_Plus_2Ki.coord==0, other=np.nan)
    ds_elm_Plus_2Ki = add_shape_coord_from_data_array(ds_elm_Plus_2Ki, "./shapes/srb.shp", 'coord' )
    ds_elm_Plus_2Ki = ds_elm_Plus_2Ki.where(ds_elm_Plus_2Ki.coord==0, other=np.nan)
    ds_mos_Plus_2Ki = add_shape_coord_from_data_array(ds_mos_Plus_2Ki, "./shapes/srb.shp", 'coord' )
    ds_mos_Plus_2Ki = ds_mos_Plus_2Ki.where(ds_mos_Plus_2Ki.coord==0, other=np.nan)

    ds_eam_Plus_3Ki = add_shape_coord_from_data_array(ds_eam_Plus_3Ki, "./shapes/srb.shp", 'coord' )
    ds_eam_Plus_3Ki = ds_eam_Plus_3Ki.where(ds_eam_Plus_3Ki.coord==0, other=np.nan)
    ds_elm_Plus_3Ki = add_shape_coord_from_data_array(ds_elm_Plus_3Ki, "./shapes/srb.shp", 'coord' )
    ds_elm_Plus_3Ki = ds_elm_Plus_3Ki.where(ds_elm_Plus_3Ki.coord==0, other=np.nan)
    ds_mos_Plus_3Ki = add_shape_coord_from_data_array(ds_mos_Plus_3Ki, "./shapes/srb.shp", 'coord' )
    ds_mos_Plus_3Ki = ds_mos_Plus_3Ki.where(ds_mos_Plus_3Ki.coord==0, other=np.nan)

    ds_eam_Plus_4Ki = add_shape_coord_from_data_array(ds_eam_Plus_4Ki, "./shapes/srb.shp", 'coord' )
    ds_eam_Plus_4Ki = ds_eam_Plus_4Ki.where(ds_eam_Plus_4Ki.coord==0, other=np.nan)
    ds_elm_Plus_4Ki = add_shape_coord_from_data_array(ds_elm_Plus_4Ki, "./shapes/srb.shp", 'coord' )
    ds_elm_Plus_4Ki = ds_elm_Plus_4Ki.where(ds_elm_Plus_4Ki.coord==0, other=np.nan)
    ds_mos_Plus_4Ki = add_shape_coord_from_data_array(ds_mos_Plus_4Ki, "./shapes/srb.shp", 'coord' )
    ds_mos_Plus_4Ki = ds_mos_Plus_4Ki.where(ds_mos_Plus_4Ki.coord==0, other=np.nan)

    ds_eam_Pre_Indi = add_shape_coord_from_data_array(ds_eam_Pre_Indi, "./shapes/srb.shp", 'coord' )
    ds_eam_Pre_Indi = ds_eam_Pre_Indi.where(ds_eam_Pre_Indi.coord==0, other=np.nan)
    ds_elm_Pre_Indi = add_shape_coord_from_data_array(ds_elm_Pre_Indi, "./shapes/srb.shp", 'coord' )
    ds_elm_Pre_Indi = ds_elm_Pre_Indi.where(ds_elm_Pre_Indi.coord==0, other=np.nan)
    ds_mos_Pre_Indi = add_shape_coord_from_data_array(ds_mos_Pre_Indi, "./shapes/srb.shp", 'coord' )
    ds_mos_Pre_Indi = ds_mos_Pre_Indi.where(ds_mos_Pre_Indi.coord==0, other=np.nan)

#creating an array of time values for the xtick labels in the graphs
def datetime_range(start, end, delta):
    current = start
    if not isinstance(delta, timedelta):
        delta = timedelta(**delta)
    while current < end:
        yield current
        current += delta

start = datetime.datetime(1996,1,15,6)
end = datetime.datetime(1996,1,24,3)
start2 = datetime.date(1996,1,16)
end2 = datetime.date(1996,1,26)
end3 = datetime.date(1996,1,25)
time = []
time3 = []
time2 = []
for dt in datetime_range(start, end, {'hours':6}):
    time.append(dt)
for dt in datetime_range(start2, end3, {'days':1}):
    time3.append('')
    time3.append('')
    time3.append('')
    time3.append(dt.strftime('%b %d'))
for dt in datetime_range(start2, end2, {'days':1}):
    time2.append(dt)


# ## Figure 2

##plots a timeseries of the 1996 event

#calculating the mean of each variable over the basin at each time period
if SRB_AVG == False:
    mask = (
    (ds_elm_Control.coords["lat"] > x)
    & (ds_elm_Control.coords["lat"] < x+1)
    & (ds_elm_Control.coords["lon"] > y)
    & (ds_elm_Control.coords["lon"] < y+1)
        )

    ds_mos1 = ds_mos_Control.where(mask)
    q = ds_elm_Control['QRUNOFF'].where(mask)
    p = ds_eam_Control['PRECT'].where(mask)
    t = ds_eam_Control['TREFHT'].where(mask)
    s = ds_elm_Control['H2OSNO'].where(mask)
    q_mean = q.mean( dim=("lat","lon"), skipna=True )
    p_mean = p.mean( dim=("lat","lon"), skipna=True )
    t_mean = t.mean( dim=("lat","lon"), skipna=True )
    s_mean = s.mean( dim=("lat","lon"), skipna=True )
    ds_mos_mean = ds_mos1.mean( dim=("lat","lon"), skipna=True )
    Time1 = Time_elm
    Time2 = Time_eam

if SRB_AVG == True:    

    q = ds_elm_Control['QRUNOFF']
    p = ds_eam_Control['PRECT']
    t = ds_eam_Control['TREFHT']
    s = ds_elm_Control['H2OSNO']
    q_mean = q.mean( dim=("lat","lon"), skipna=True )
    p_mean = p.mean( dim=("lat","lon"), skipna=True )
    t_mean = t.mean( dim=("lat","lon"), skipna=True )
    s_mean = s.mean( dim=("lat","lon"), skipna=True )
    ds_mos_mean = ds_mos_Control.mean( dim=("lat","lon"), skipna=True )
    Time1 = Time_elm
    Time2 = Time_eam
            
    ## Setup plot
    fig = plt.figure(figsize=(14,5))
    host1 = fig.add_subplot(111)

    ## Create extra y-axes that shares x-axis with "host"
    yaxis1 = host1.twinx()
    yaxis2 = host1.twinx()
    yaxis3 = host1.twinx()

    ## Set axis limits
    host1.set_ylim([-5, 140])
    yaxis1.set_ylim([-5, 250])
    yaxis2.set_ylim([-25,10])
    yaxis3.set_ylim([-5,3000])

    # Offset the right spine and show it
    yaxis2.spines["right"].set_position(("axes", 1.1))
    yaxis2.spines["right"].set_visible(True)
    yaxis3.spines["right"].set_position(("axes", 1.2))
    yaxis3.spines["right"].set_visible(True)

    ## Set labels for axes
    host1.set_ylabel("Snow Water Equivalent (mm)", color='c', fontsize='x-large')
    yaxis1.set_ylabel("Runoff & Precipitation (mm day$^{-1}$)", color='b',fontsize='x-large')
    yaxis2.set_ylabel("2-m Temperature ($^\circ$C)", color='r',fontsize='x-large')
    yaxis3.set_ylabel('Discharge m${^3}$ s$^{-1}$', color='m',fontsize='x-large')

    ## Select colors
    color1 = 'c'
    color2 = 'b'
    color3 = 'g'
    color4 = 'r'
    color5 = 'm'

    ## Plot each line

    p1, = host1.plot(Time1, s_mean, '-', color=color1,label="SWE")
    p2, = yaxis1.plot(Time1, q_mean, '-', color=color2, label="Runoff")
    p3, = yaxis1.plot(Time2, p_mean, '-', color=color3, label="Precipitation")
    p4, = yaxis2.plot(Time2, t_mean, '-', color=color4, label="Temperature")
    p6, = yaxis2.plot(Time2, np.zeros(len(Time2)), color='silver',linestyle='--', label='0$^\circ$C line')
    p5, = yaxis3.plot(Time_mos[:-2], ds_mos_mean['RIVER_DISCHARGE_OVER_LAND_LIQ'][:-2], color=color5,label='River Discharge')

    ## Add legend, specify loc with 2x tuple of bottom left corner of box in 0->1 coords
    lns = [p1, p2, p3, p4, p5, p6]
    host1.legend(handles=lns, fontsize='medium', ncol=2)
    host1.set_xticks(time)
    host1.set_xticklabels(time3)
    host1.margins(x=0)
    host1.tick_params('both', labelsize='large')
    yaxis1.tick_params('y', labelsize='large')
    yaxis2.tick_params('y', labelsize='large')
    yaxis3.tick_params('y', labelsize='large')
    fig.tight_layout(pad=3.0)
    
fig.savefig('event_hydrograph.pdf', bbox_inches='tight', pad_inches=0)







## Fig ????
## CMZ load some stuff

swe_obs_ds = xr.open_dataset("/Users/cmz5202/Software/1996-ros-e3sm/ncl/gen-nc-files-for-obs/obs_swe.nc")
swe_obs = swe_obs_ds.SWE_avg
swe_obs_time = swe_obs_ds.time

swe_obs2_ds = xr.open_dataset("/Users/cmz5202/Software/1996-ros-e3sm/ncl/gen-nc-files-for-obs/ERA5_swe.nc")
swe_obs2 = swe_obs2_ds.SWE_avg
swe_obs2_time = swe_obs2_ds.time

swe_obs3_ds = xr.open_dataset("/Users/cmz5202/Software/1996-ros-e3sm/ncl/gen-nc-files-for-obs/JRA_swe.nc")
swe_obs3 = swe_obs3_ds.SWE_avg
swe_obs3_time = swe_obs3_ds.time

swe_obs4_ds = xr.open_dataset("/Users/cmz5202/Software/1996-ros-e3sm/ncl/gen-nc-files-for-obs/NLDAS_swe.nc")
swe_obs4 = swe_obs4_ds.SWE_avg
swe_obs4_time = swe_obs4_ds.time

precip_mod_ds = xr.open_dataset("/Users/cmz5202/Software/1996-ros-e3sm/ncl/gen-nc-files-for-obs/mod_precip.nc")
precip_mod = precip_mod_ds.SWE_avg
precip_mod_time = precip_mod_ds.time

precip_obs_ds = xr.open_dataset("/Users/cmz5202/Software/1996-ros-e3sm/ncl/gen-nc-files-for-obs/obs_precip.nc")
precip_obs = precip_obs_ds.SWE_avg
precip_obs_time = precip_obs_ds.time

tmax_mod_ds = xr.open_dataset("/Users/cmz5202/Software/1996-ros-e3sm/ncl/gen-nc-files-for-obs/mod_tmax.nc")
tmax_mod = tmax_mod_ds.SWE_avg
tmax_mod = tmax_mod - 273.15
tmax_mod_time = tmax_mod_ds.time

tmax_obs_ds = xr.open_dataset("/Users/cmz5202/Software/1996-ros-e3sm/ncl/gen-nc-files-for-obs/obs_tmax.nc")
tmax_obs = tmax_obs_ds.SWE_avg
tmax_obs = tmax_obs - 273.15
tmax_obs_time = tmax_obs_ds.time

tmin_mod_ds = xr.open_dataset("/Users/cmz5202/Software/1996-ros-e3sm/ncl/gen-nc-files-for-obs/mod_tmin.nc")
tmin_mod = tmin_mod_ds.SWE_avg
tmin_mod = tmin_mod - 273.15
tmin_mod_time = tmin_mod_ds.time

tmin_obs_ds = xr.open_dataset("/Users/cmz5202/Software/1996-ros-e3sm/ncl/gen-nc-files-for-obs/obs_tmin.nc")
tmin_obs = tmin_obs_ds.SWE_avg
tmin_obs = tmin_obs - 273.15
tmin_obs_time = tmin_obs_ds.time

#calculating the mean of each variable over the basin at each time period
if SRB_AVG == False:
    mask = (
    (ds_elm_Control.coords["lat"] > x)
    & (ds_elm_Control.coords["lat"] < x+1)
    & (ds_elm_Control.coords["lon"] > y)
    & (ds_elm_Control.coords["lon"] < y+1)
        )

    ds_mos1 = ds_mos_Control.where(mask)
    q = ds_elm_Control['QRUNOFF'].where(mask)
    p = ds_eam_Control['PRECT'].where(mask)
    t = ds_eam_Control['TREFHT'].where(mask)
    s = ds_elm_Control['H2OSNO'].where(mask)
    q_mean = q.mean( dim=("lat","lon"), skipna=True )
    p_mean = p.mean( dim=("lat","lon"), skipna=True )
    t_mean = t.mean( dim=("lat","lon"), skipna=True )
    s_mean = s.mean( dim=("lat","lon"), skipna=True )
    ds_mos_mean = ds_mos1.mean( dim=("lat","lon"), skipna=True )
    Time1 = Time_elm
    Time2 = Time_eam

if SRB_AVG == True:    

    q = ds_elm_Control['QRUNOFF']
    p = ds_eam_Control['PRECT']
    t = ds_eam_Control['TREFHT']
    s = ds_elm_Control['H2OSNO']
    q_mean = q.mean( dim=("lat","lon"), skipna=True )
    p_mean = p.mean( dim=("lat","lon"), skipna=True )
    t_mean = t.mean( dim=("lat","lon"), skipna=True )
    s_mean = s.mean( dim=("lat","lon"), skipna=True )
    ds_mos_mean = ds_mos_Control.mean( dim=("lat","lon"), skipna=True )
    Time1 = Time_elm
    Time2 = Time_eam
            
    ## Setup plot
    fig = plt.figure(figsize=(10.75,5))
    host1 = fig.add_subplot(111)

    ## Create extra y-axes that shares x-axis with "host"
    yaxis1 = host1.twinx()
    yaxis2 = host1.twinx()

    ## Set axis limits
    host1.set_ylim([-5, 160])
    yaxis1.set_ylim([-5, 40])
    yaxis2.set_ylim([-25,25])
    
    #host1.set_xlim([datetime.datetime(1996, 1, 15), datetime.datetime(1996, 1, 22)])
    #yaxis1.set_xlim([datetime.datetime(1996, 1, 15), datetime.datetime(1996, 1, 22)])
    #yaxis2.set_xlim([datetime.datetime(1996, 1, 15), datetime.datetime(1996, 1, 22)])
    
    # Offset the right spine and show it
    yaxis2.spines["right"].set_position(("axes", 1.1))
    yaxis2.spines["right"].set_visible(True)

    ## Set labels for axes
    host1.set_ylabel("Snow Water Equivalent (mm)", color='c', fontsize='x-large')
    yaxis1.set_ylabel("Precipitation (mm day$^{-1}$)", color='g',fontsize='x-large')
    yaxis2.set_ylabel("2-m Temperature ($^\circ$C)", color='r',fontsize='x-large')

    ## Select colors
    color1 = 'c'
    color2 = 'b'
    color3 = 'g'
    color4 = 'darksalmon'
    color5 = 'm'
    color6 = 'darkred' # added for tmin

    ## Plot each line
    
    p1, = host1.plot(Time1, s_mean, '-', color=color1,label="E3SM SWE", linewidth=2.7)
    p2, = host1.plot(swe_obs_time, swe_obs, '--', color=color1,label="UA SWE", linewidth=2)
    p10, = host1.plot(swe_obs2_time, swe_obs2, '-.', color=color1,label="ERA5 SWE", linewidth=2)
    p11, = host1.plot(swe_obs3_time, swe_obs3, ':', color=color1,label="JRA SWE", linewidth=2)
    p12, = host1.plot(swe_obs4_time, swe_obs4, linestyle=(0, (3, 1, 1, 1, 1, 1)), color=color1,label="NLDAS SWE", linewidth=2.5)

    # temps
    p3, = yaxis2.plot(tmax_obs_time, tmax_obs, '--', color=color4,label="CPC Max. Temp.", linewidth=2.0)
    p4, = yaxis2.plot(tmax_mod_time, tmax_mod, '-', color=color4,label="E3SM Max. Temp.", linewidth=2.7)
    p5, = yaxis2.plot(tmin_obs_time, tmin_obs, '--', color=color6,label="CPC Min. Temp.", linewidth=1.5)
    p6, = yaxis2.plot(tmin_mod_time, tmin_mod, '-', color=color6,label="E3SM Min. Temp.", linewidth=2.2)
        
    p7, = yaxis1.plot(precip_obs_time, precip_obs, '--', color=color3,label="CPC Precip.", linewidth=1.5)
    p8, = yaxis1.plot(precip_mod_time, precip_mod, '-', color=color3,label="E3SM Precip.", linewidth=2.2)
    
    p9, = yaxis2.plot(Time2, np.zeros(len(Time2)), color='silver',linestyle=':', label='0$^\circ$C line')

    ## Add legend, specify loc with 2x tuple of bottom left corner of box in 0->1 coords
    lns = [p6, p5, p4, p3, p8, p7, p1, p2, p10, p11, p12 ]
    host1.legend(handles=lns, fontsize='small', ncol=2,handlelength=3)
    host1.set_xticks(time)
    host1.set_xticklabels(time3)
    host1.margins(x=0)
    host1.tick_params('both', labelsize='large')
    yaxis1.tick_params('y', labelsize='large')
    yaxis2.tick_params('y', labelsize='large')
    fig.tight_layout(pad=3.0)
    
    #fig.xlim(xmin=datetime.datetime(1996, 1, 15),xmax=datetime.datetime(1996, 1, 22))
    
fig.savefig('comp_event_hydrograph.pdf', bbox_inches='tight', pad_inches=0)









# ## Figure 8

### This cell plots a timeseries of a chosen variables for the 6 different counterfactual simulations with
### the 25th to 75th percentiles shaded

#set shading = True if want to plot the 25th and 75th percentiles on timeseries set =False if don't want shading
shading = True

#setting basin outline incase SRB = False
minlat = 40
maxlat = 42.5
minlon = -78.5
maxlon = -76
#creating array of variables
eam_var = ['TREFHT','PRECT', 'CAPE', 'OMEGA500', 'OMEGA850', 'PSL', 'TMQ']
elm_var = ['QRUNOFF','H2OSNO', 'H2OSFC', 'TSA', 'TWS', 'FGR']
mos_var = ['RIVER_DISCHARGE_OVER_LAND_LIQ']
variables = ['TREFHT','PRECT','H2OSNO','QRUNOFF','RIVER_DISCHARGE_OVER_LAND_LIQ']
#array of simulations
EAM_runs = [ds_eam_Pre_Indi, ds_eam_Controli, ds_eam_Plus_1Ki, ds_eam_Plus_2Ki, ds_eam_Plus_3Ki,ds_eam_Plus_4Ki]
ELM_runs = [ds_elm_Pre_Indi, ds_elm_Controli, ds_elm_Plus_1Ki, ds_elm_Plus_2Ki, ds_elm_Plus_3Ki,ds_elm_Plus_4Ki]
MOS_runs = [ds_mos_Pre_Indi, ds_mos_Controli, ds_mos_Plus_1Ki, ds_mos_Plus_2Ki, ds_mos_Plus_3Ki,ds_mos_Plus_4Ki]
#creating arrays for loops
runs = ['Pre_Ind', 'Control', 'Plus_1K', 'Plus_2K', 'Plus_3K','Plus_4K']
colors = ['c','b','g','r','m','y']
#craeting empty arrays
max_snow = []
min_snow = []
delta_swe = []
#setting up plot
fig, axs = plt.subplots(5,1,figsize=(15,22))

#loops through variables and plots them
for ax,vari,L in zip(axs,variables,letters):
    #creates labels for the plots based on variable chosen above
    if vari == 'TREFHT':
        titl = '2-m Temperature'
        label = '2-m Temperature ($^\circ$C)'
    if vari == 'QRUNOFF':
        titl = 'Runoff'
        label = 'Runoff (mm day$^{-1}$)'
    if vari == 'PRECT':
        titl = 'Precipitation'
        label = 'Precipitation (mm day$^{-1}$)'
    if vari == 'H2OSNO':
        titl = 'Snow water equivalent'
        label = 'Snow water equivalent (mm)'
    if vari == 'RIVER_DISCHARGE_OVER_LAND_LIQ':
        titl = 'River Discharge'
        label = 'River Discharge (m${^3}$ s$^{-1}$)'
        
    #averaging the variables over the basin and plotting it
    for a,l,m,c,r in zip(EAM_runs,ELM_runs,MOS_runs,colors,runs):
        if SRB_AVG == False:
            mask = (
            (ds_elm_Controli.coords["lat"] >= minlat)
            & (ds_elm_Controli.coords["lat"] <= maxlat)
            & (ds_elm_Controli.coords["lon"] >= minlon)
            & (ds_elm_Controli.coords["lon"] <= maxlon))

            if vari in elm_var:
                v = l[vari].where(mask)
                Time = Time_elm
            if vari in eam_var:
                v = a[vari].where(mask)
                Time = Time_eam
            if vari in mos_var:
                v = m[vari].where(mask)
                Time = Time_mos
        else: 
            if vari in elm_var:
                v = l[vari]
                Time = Time_elm
            if vari in eam_var:
                v = a[vari]
                Time = Time_eam
            if vari in mos_var:
                v = m[vari]
                Time = Time_mos

        v_mean = v.mean( dim=("lat","lon"), skipna=True )
        v_mean2 = v.mean(dim=("lat","lon", 'init'),skipna=True )
        
        #calculating delta SWE if vari is H2OSNO
        if vari == 'H2OSNO':
                max_snow.append(v_mean2.max())
                min_snow.append(v_mean2.min())
                delta_swe.append(v_mean2.max()-v_mean2.min())
                
        #creating arrays for the 25th and 75th percentiles
        if shading == True:
            twenty_fifth = []
            seventy_fifth = []
            for i in range(len(Time)):
                twenty_fifth.append(np.nanpercentile(v_mean[:,i],25))
                seventy_fifth.append(np.nanpercentile(v_mean[:,i],75))
            #plotting the shading
            ax.fill_between(Time[0:], twenty_fifth[0:], seventy_fifth[0:], alpha=0.2, color=c)
            shading_titl = ' (25th-75th percentile shaded)'
        else:
            shading_titl = ''

        #plotting the lines
        if vari == 'H2OSNO':
            ax.set_ylim([-5,150])
        ax.plot(Time[0:], v_mean2[0:], '-', color=c, label=r)
    
    #adding stuff to plot
    ax.margins(x=0)
    ax.text(.02,.97,L+' Basin Averaged ' + titl , fontsize='xx-large',transform=ax.transAxes,ha='left',va='top')#,bbox=props)
    ax.set_xticks(time)
    ax.tick_params('both', labelsize='x-large')
    ax.tick_params('x',direction='inout',length=18)
    ax.set_ylabel(label, fontsize='x-large')
    ax.grid(color='gainsboro',linewidth=.5)
    
    #adding legend to precip subplot
    if vari == 'PRECT':
        ax.legend(loc=1, fontsize='xx-large')
        
    #adding xtick labels to bottom subplot only
    if vari == 'RIVER_DISCHARGE_OVER_LAND_LIQ':
        ax.set_xticklabels(time3, fontsize='large')
    else:
        ax.set_xticklabels('')
    #making subplots touch (no space)
    plt.subplots_adjust(hspace=10**-6)

fig.savefig('timeseries_counterfactuals.pdf', bbox_inches='tight', pad_inches=0)


##calculates statistcs (mean and/or delta) for each variable and obtains the spread of ensemble points

#arrays of simulations
EAM_runs = [ds_eam_Pre_Indi, ds_eam_Controli, ds_eam_Plus_1Ki, ds_eam_Plus_2Ki, ds_eam_Plus_3Ki,ds_eam_Plus_4Ki]
ELM_runs = [ds_elm_Pre_Indi, ds_elm_Controli, ds_elm_Plus_1Ki, ds_elm_Plus_2Ki, ds_elm_Plus_3Ki,ds_elm_Plus_4Ki]
MOS_runs = [ds_mos_Pre_Indi, ds_mos_Controli, ds_mos_Plus_1Ki, ds_mos_Plus_2Ki, ds_mos_Plus_3Ki,ds_mos_Plus_4Ki]

#variables
eam_v = ['PRECT', 'TREFHT']
elm_v = ['H2OSNO', 'QRUNOFF']
mos_v = ['RIVER_DISCHARGE_OVER_LAND_LIQ']
varis = ['PRECT', 'H2OSNO', 'TREFHT', 'QRUNOFF','RIVER_DISCHARGE_OVER_LAND_LIQ']

#creating empty arrays for loop to fill in
sno_diff = []
avg_temp = []
avg_pre = []
avg_run = []
avg_dis = []

sno_diffi = []
avg_tempi = []
avg_prei = []
avg_runi = []
avg_disi = []

#looping through variables and calcluating statistics to add to empty arrays
for va in varis:
    if va in eam_v:
        for l in EAM_runs:
            if SRB_AVG == False:
                mask = (
                (ds_elm_Controli.coords["lat"] >= minlat)
                & (ds_elm_Controli.coords["lat"] <= maxlat)
                & (ds_elm_Controli.coords["lon"] >= minlon)
                & (ds_elm_Controli.coords["lon"] <= maxlon))

                v = l[va].where(mask)
                v_mean = v.mean( dim=("lat","lon"), skipna=True )  
                if va == 'PRECT':
                    avg_pre.append(v_mean.mean())
                    for x in range(len(ds_elm_Controli['init'])):
                        avg_prei.append(v_mean[x,:].mean(skipna=True))

                if va == 'TREFHT':
                    avg_temp.append(v_mean.mean())
                    for x in range(len(ds_elm_Controli['init'])):
                        avg_tempi.append(v_mean[x,:].mean())
            else:
                #averaging across basin (lat/lon)
                v = l[va]
                v_mean = v.mean( dim=("lat","lon"), skipna=True )  
                
                #calculating average of variable across the basin and init time
                if va == 'PRECT':
                    avg_pre.append(v_mean.mean())
                    #calculating average of variable for each init time
                    for x in range(len(ds_elm_Controli['init'])):
                        avg_prei.append(v_mean[x,:].mean(skipna=True))

                if va == 'TREFHT':
                    avg_temp.append(v_mean.mean())
                    for x in range(len(ds_elm_Controli['init'])):
                        avg_tempi.append(v_mean[x,:].mean())
    
    #doing the same thing as above for elm variables 
    if va in elm_v:
        for l in ELM_runs:
            if SRB_AVG ==False:
                mask = (
                (ds_elm_Controli.coords["lat"] >= minlat)
                & (ds_elm_Controli.coords["lat"] <= maxlat)
                & (ds_elm_Controli.coords["lon"] >= minlon)
                & (ds_elm_Controli.coords["lon"] <= maxlon))

                v = l[va].where(mask)
                v_mean = v.mean( dim=("lat","lon"), skipna=True )
                v_mean2 = v_mean.mean(dim='init', skipna=True)
                if va == 'H2OSNO':
                    sno_diff.append(v_mean2.max() - v_mean2.min())
                    for x in range(len(ds_elm_Controli['init'])):
                        sno_diffi.append(v_mean[x,:].max()-v_mean[x,:].min())
        
                if va == 'QRUNOFF':
                    avg_run.append(v_mean.mean())
                    for x in range(len(ds_elm_Controli['init'])):
                        avg_runi.append(v_mean[x,:].mean())
                        
            else:
                v = l[va]
                v_mean = v.mean( dim=("lat","lon"), skipna=True )
                v_mean2 = v_mean.mean(dim='init', skipna=True) #calculating second mean acorss all init times
                if va == 'H2OSNO':
                    #calculating delta SWE for ensemble avg
                    sno_diff.append(v_mean2.max() - v_mean2.min())
                    for x in range(len(ds_elm_Controli['init'])):
                        #calculating delta SWE for each ensemble member
                        sno_diffi.append(v_mean[x,:].max()-v_mean[x,:].min())
                if va == 'QRUNOFF':
                    #calculating average runoff for ensemble
                    avg_run.append(v_mean.mean())
                    for x in range(len(ds_elm_Controli['init'])):
                        #calculating average runoff for each ensemble member
                        avg_runi.append(v_mean[x,:].mean())

                        
    if va in mos_v:
        for l in MOS_runs:
            if SRB_AVG ==False:
                mask = (
                (ds_mos_Controli.coords["lat"] >= minlat)
                & (ds_mos_Controli.coords["lat"] <= maxlat)
                & (ds_mos_Controli.coords["lon"] >= minlon)
                & (ds_mos_Controli.coords["lon"] <= maxlon))

                v = l[va].where(mask)
                v_mean = v.mean( dim=("lat","lon"), skipna=True )
                if va == 'RIVER_DISCHARGE_OVER_LAND_LIQ':
                    avg_dis.append(v_mean.mean())
                    for x in range(len(ds_mos_Controli['init'])):
                        avg_disi.append(v_mean[x,:].mean())
                        
            else:
                v = l[va]
                v_mean = v.mean( dim=("lat","lon"), skipna=True )
                if va == 'RIVER_DISCHARGE_OVER_LAND_LIQ':
                    avg_dis.append(v_mean.mean())
                    for x in range(len(ds_mos_Controli['init'])):
                        avg_disi.append(v_mean[x,:].mean())


# ## Figure 9

##plots the avg temp, avg pre, avg, runoff, and delta swe for each of the 6 simulations showing the trends for future climates

#creating spacing for xticks (chose .3 bc 1996 was roughly .3C above pre industrial)
temps = [0, .3 ,1, 2, 3, 4]
#creating labels for figure
xtik = ['PI', 'Control', '+1K', '+2K', '+3K', '+4K']

#setting up plot
fig = plt.figure(figsize=(13,5))
host1 = fig.add_subplot(111)

## Create extra y-axes that shares x-axis with "host"
yaxis1 = host1.twinx()
yaxis2 = host1.twinx()
yaxis3 = host1.twinx()
yaxis4 = host1.twinx()

# Offsetting the spines and and show it
yaxis4.spines["left"].set_position(("axes", -0.1))
yaxis4.spines["left"].set_visible(True)   
yaxis4.yaxis.set_label_position('left')
yaxis4.yaxis.set_ticks_position('left')
yaxis2.spines["right"].set_position(("axes", 1.1))
yaxis2.spines["right"].set_visible(True)
yaxis3.spines["right"].set_position(("axes", 1.2))
yaxis3.spines["right"].set_visible(True)

#plotting variables
p1, = host1.plot(temps,sno_diff, linestyle='dashed', marker='*', markersize=9, color='b', label = '$\Delta$'+' SWE')
p2, = yaxis1.plot(temps,avg_run, linestyle='dashed', marker='^',markersize=9,color='g', label = 'Runoff')
p5, = yaxis4.plot(temps,avg_temp, linestyle='dashed', marker='o',markersize=9,color='r', label = 'Temperature' )
p3, = yaxis2.plot(temps,avg_dis, linestyle='dashed', marker='D',markersize=9,color='grey', label = 'River Discharge' )
p4, =  yaxis3.plot(temps,avg_pre, linestyle='dashed', marker='x',markersize=9,color = 'c', label =  'Precipitation')
lns = [p1, p2, p3, p4, p5]

## Set labels for axes
host1.set_ylabel("SWE difference (mm)",fontsize=15, color='b');
yaxis1.set_ylabel("Average Runoff (mm day$^{-1}$)",fontsize=15, color='g');
yaxis4.set_ylabel("Average 2-m temperature ($^\circ$C)",fontsize=15, color='r');
yaxis2.set_ylabel("Average River Discharge (m${^3}$ s$^{-1}$)",fontsize=15, color='grey');
yaxis3.set_ylabel('Average Precipitation (mm day$^{-1}$)',fontsize=15, color='c');

#adding stuff to plot
plt.legend(handles=lns,loc=8, bbox_to_anchor=(0.42, 0., 0.5, 0.5), fontsize='large' );
host1.set_xticks(temps)
host1.set_xticklabels(xtik, fontsize='xx-large');
host1.tick_params('both', labelsize='x-large')
yaxis1.tick_params('y', labelsize='x-large')
yaxis2.tick_params('y', labelsize='x-large')
yaxis3.tick_params('y', labelsize='x-large')
yaxis4.tick_params('y', labelsize='x-large')

fig.savefig('avgchange_vs_warming_line.pdf', bbox_inches='tight', pad_inches=0)


# ## Figure 10

#plots the spread of the ensemble for each of the variables

#creating array in order to plot points below 
temps2 = np.repeat(temps,15)
#array of non ensemble averaged variable stats
init_dots=[avg_tempi,avg_prei,sno_diffi,avg_runi, avg_disi]  
#array of ensemble averaged variable stats
var = [avg_temp,avg_pre,sno_diff,avg_run, avg_dis]

#labels for plot
titl=['Basin Average Temperature', 'Basin Average Precipitation', 'Delta SWE', 'Basin Average Runoff', 'Basin Average River Discharge']
ylab = ['Temperature ($^\circ$C)','Precipitation (mm day$^{-1}$)','Delta SWE (mm)', 'Runoff (mm day$^{-1}$)','River Discharge (m${^3}$ s$^{-1}$)']
colors = ['r', 'c','b', 'g', 'grey']

#setting up plot
fig = plt.figure(figsize=(20,13))
ax1 = plt.subplot2grid(shape=(3,4), loc=(0,0), colspan=2)
ax2 = plt.subplot2grid((3,4), (0,2), colspan=2)
ax3 = plt.subplot2grid((3,4), (1,0), colspan=2)
ax4 = plt.subplot2grid((3,4), (1,2), colspan=2)
ax5 = plt.subplot2grid((3,4), (2,1), colspan=2)
plt.subplots_adjust(hspace=0.3,wspace=0.3)
axs = [ax1,ax2,ax3,ax4,ax5]

#looping through variables and plotting them in subplots
for v,i,c,t,y,ax,L in zip(var,init_dots,colors,titl,ylab,axs,letters): #axs.ravel()):
    ax.scatter(temps2, i,marker='o', color=c) #plotting non ensemble averaged point
    ax.plot(temps, v, color=c) #plotting ensemble average
    #calculating 25-75 ensemble spread
    twenty_fifth = []
    seventy_fifth = []
    for x in np.arange(0,90,15):
        twenty_fifth.append(np.nanpercentile(i[x:x+15],25))
        seventy_fifth.append(np.nanpercentile(i[x:x+15],75))
        
    #plotting the shading
    ax.fill_between(temps, twenty_fifth[0:], seventy_fifth[0:], alpha=0.2, color=c)
    shading_titl = ' (25th-75th percentile shaded)'
    
    #adding stuff to plot
    ax.set_ylabel(y, fontsize='x-large')
    ax.set_title(L+ ' '+t, fontsize='xx-large')
    ax.set_xticks(temps)
    ax.set_xticklabels(xtik, fontsize='x-large');
    plt.tight_layout()

fig.savefig('enschange_vs_warming_line.pdf', bbox_inches='tight', pad_inches=0)


#calculating the ensemble averaged delta SWE by averaging across the delta swe for each ensemble member for a simulation
snow_diff2 = []
for i in range(6):
    snow_diff2.append(np.mean(sno_diffi[i*15:(i*15)+15]))


# ## Figure 16

##plots comparison of two averaging techniques for delta swe - part a

#setting up arrays for xtixk labels
temps = [0, .3 ,1, 2, 3, 4]
xtik = ['PI', 'Control', '+1K', '+2K', '+3K', '+4K']

#setting up plot
fig = plt.figure(figsize=(12,4))
ax = fig.add_subplot(111)

#plotting the two different averaging techniques for delta swe
p1, = ax.plot(temps,sno_diff, linestyle='dashed', marker='*', markersize=9, color='b', label = 'Average 1')
p2, = ax.plot(temps,snow_diff2, linestyle='dashed', marker='^',markersize=9,color='k', label = 'Average 2')
lns = [p1, p2]

#Set labels for axes
ax.set_ylabel("SWE difference (mm)",fontsize='large');
ax.set_xticks(temps)
ax.set_xticklabels(xtik);
ax.tick_params('both', labelsize='large')
#plot legent
plt.legend(handles=lns,loc=1, fontsize='large' );

#putting subplot labels on
props = dict(boxstyle='square', facecolor='white', edgecolor='k')
ax.text(.01,.97,letters[0], fontsize='xx-large',transform=ax.transAxes,ha='left',va='top',bbox=props)

fig.savefig('dSWE_avg_comp.pdf', bbox_inches='tight', pad_inches=0)


##plots comparison of two averaging techniques for delta swe - parts b/c

#setting up arrays of variables for loops
init_dots=[sno_diffi,sno_diffi]  
var = [sno_diff, snow_diff2]
#arrays for labeling
titl=['Delta SWE - Average 1','Delta SWE - Average 2']
ylab = ['Delta SWE (mm)','Delta SWE (mm)']
colors = ['b','k']
#setting up figure
fig = plt.figure(figsize=(17,8.5))
ax1 = plt.subplot2grid(shape=(3,4), loc=(0,0), colspan=2)
ax2 = plt.subplot2grid((3,4), (0,2), colspan=2)
axs = [ax1,ax2]
host1.set_xticks(temps)
host1.set_xticklabels(xtik);
host1.tick_params('both', labelsize='x-large')

#looping through variables and plotting line/dots/shading
for v,i,c,t,y,ax,j in zip(var,init_dots,colors,titl,ylab,axs,[1,2]): #axs.ravel()):
    ax.scatter(temps2, i,marker='o', color=c)
    ax.plot(temps, v, color=c)
    #calculating ensemble spread
    twenty_fifth = []
    seventy_fifth = []
    for x in np.arange(0,90,15):
        twenty_fifth.append(np.nanpercentile(i[x:x+15],25))
        seventy_fifth.append(np.nanpercentile(i[x:x+15],75))
    #plotting the shading
    ax.fill_between(temps, twenty_fifth[0:], seventy_fifth[0:], alpha=0.2, color=c)
    shading_titl = ' (25th-75th percentile shaded)'
    #adding labels
    ax.set_ylabel(y, fontsize='x-large')
    ax.set_xticks(temps)
    ax.set_xticklabels(xtik, fontsize='x-large');
    #adding subplot label
    props = dict(boxstyle='square', facecolor='white', edgecolor='k')
    ax.text(.01,.97,letters[j], fontsize='xx-large',transform=ax.transAxes,ha='left',va='top',bbox=props)
    plt.tight_layout()

fig.savefig('dSWE_avg_comp_2.pdf', bbox_inches='tight', pad_inches=0)
