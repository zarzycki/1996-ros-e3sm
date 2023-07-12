#!/usr/bin/env python
# coding: utf-8

# # Creates Figure 11, Figure 12, Figure 13, Figure 15, Figure 14, Figure 17, & Figure 18

#importing libraries
import numpy as np
import netCDF4
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xesmf as xe
import datetime
from datetime import timedelta
import shapefile as shp
import geopandas as gpd
import sys

### CMZ Command line args
basepath = sys.argv[1]


#set to true to use shapefile to average over SRB only, set to false to use rough basin average
SRB_AVG = True
#setting up rough outline of SRB if SRB_AVG set to false
if SRB_AVG == False:
    minlat = 40
    maxlat = 42.5
    minlon = -78.5
    maxlon = -76

#ELM and MOSART variables are time averaged, so to shift the timeseries plots to fall at the middle
#of the time period rather than the end of the time period set equal to True rather than False
shift_ELM = True
shift_MOSART = True

#setting up array of letters to put on corner of subplots
letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)','(p)','(q)','(r)','(s)','(t)','(u)',
          '(v)','(w)','(x)','(y)']

procpath = basepath+'/proc_arp/'
#loading low res files
ds_eamLR = xr.open_dataset(procpath+'Low_Res_eam_slice.nc')
ds_elmLR = xr.open_dataset(procpath+'Low_Res_elm_slice.nc')
ds_mosLR = xr.open_dataset(procpath+'Low_Res_mos_slice.nc')
#loading the non ensemble Control file
ds_eam6 = xr.open_dataset(procpath+'Control-6_eam_slice.nc')
ds_elm6 = xr.open_dataset(procpath+'Control-6_elm_slice.nc')
ds_mos6 = xr.open_dataset(procpath+'Control-6_mos_slice.nc')

#loading ensemble files
ds_eam_Control = xr.open_dataset(procpath+'ds_eam_Control_ens.nc')
ds_elm_Control = xr.open_dataset(procpath+'ds_elm_Control_ens.nc')
ds_mos_Control = xr.open_dataset(procpath+'ds_mos_Control_ens.nc')
Control_h1 = xr.open_dataset(procpath+'ds_eam_Control_h1_ens.nc')

ds_eam_Plus_1K = xr.open_dataset(procpath+'ds_eam_Plus_1K_ens.nc')
ds_elm_Plus_1K = xr.open_dataset(procpath+'ds_elm_Plus_1K_ens.nc')
ds_mos_Plus_1K = xr.open_dataset(procpath+'ds_mos_Plus_1K_ens.nc')
Plus_1K_h1 = xr.open_dataset(procpath+'ds_eam_Plus_1K_h1_ens.nc')

ds_eam_Plus_2K = xr.open_dataset(procpath+'ds_eam_Plus_2K_ens.nc')
ds_elm_Plus_2K = xr.open_dataset(procpath+'ds_elm_Plus_2K_ens.nc')
ds_mos_Plus_2K = xr.open_dataset(procpath+'ds_mos_Plus_2K_ens.nc')
Plus_2K_h1 = xr.open_dataset(procpath+'ds_eam_Plus_2K_h1_ens.nc')

ds_eam_Plus_3K = xr.open_dataset(procpath+'ds_eam_Plus_3K_ens.nc')
ds_elm_Plus_3K = xr.open_dataset(procpath+'ds_elm_Plus_3K_ens.nc')
ds_mos_Plus_3K = xr.open_dataset(procpath+'ds_mos_Plus_3K_ens.nc')
Plus_3K_h1 = xr.open_dataset(procpath+'ds_eam_Plus_3K_h1_ens.nc')

ds_eam_Plus_4K = xr.open_dataset(procpath+'ds_eam_Plus_4K_ens.nc')
ds_elm_Plus_4K = xr.open_dataset(procpath+'ds_elm_Plus_4K_ens.nc')
ds_mos_Plus_4K = xr.open_dataset(procpath+'ds_mos_Plus_4K_ens.nc')
Plus_4K_h1 = xr.open_dataset(procpath+'ds_eam_Plus_4K_h1_ens.nc')

ds_eam_Pre_Ind = xr.open_dataset(procpath+'ds_eam_Pre_Ind_ens.nc')
ds_elm_Pre_Ind = xr.open_dataset(procpath+'ds_elm_Pre_Ind_ens.nc')
ds_mos_Pre_Ind = xr.open_dataset(procpath+'ds_mos_Pre_Ind_ens.nc')
Pre_Ind_h1 = xr.open_dataset(procpath+'ds_eam_Pre_Ind_h1_ens.nc')

#loading ensemble files not averaged over initializatioin times
ds_eam_Controli = xr.open_dataset(procpath+'ds_eam_Controli_ens.nc')
ds_elm_Controli = xr.open_dataset(procpath+'ds_elm_Controli_ens.nc')
ds_mos_Controli = xr.open_dataset(procpath+'ds_mos_Controli_ens.nc')

ds_eam_Plus_1Ki = xr.open_dataset(procpath+'ds_eam_Plus_1Ki_ens.nc')
ds_elm_Plus_1Ki = xr.open_dataset(procpath+'ds_elm_Plus_1Ki_ens.nc')
ds_mos_Plus_1Ki = xr.open_dataset(procpath+'ds_mos_Plus_1Ki_ens.nc')

ds_eam_Plus_2Ki = xr.open_dataset(procpath+'ds_eam_Plus_2Ki_ens.nc')
ds_elm_Plus_2Ki = xr.open_dataset(procpath+'ds_elm_Plus_2Ki_ens.nc')
ds_mos_Plus_2Ki = xr.open_dataset(procpath+'ds_mos_Plus_2Ki_ens.nc')

ds_eam_Plus_3Ki = xr.open_dataset(procpath+'ds_eam_Plus_3Ki_ens.nc')
ds_elm_Plus_3Ki = xr.open_dataset(procpath+'ds_elm_Plus_3Ki_ens.nc')
ds_mos_Plus_3Ki = xr.open_dataset(procpath+'ds_mos_Plus_3Ki_ens.nc')

ds_eam_Plus_4Ki = xr.open_dataset(procpath+'ds_eam_Plus_4Ki_ens.nc')
ds_elm_Plus_4Ki = xr.open_dataset(procpath+'ds_elm_Plus_4Ki_ens.nc')
ds_mos_Plus_4Ki = xr.open_dataset(procpath+'ds_mos_Plus_4Ki_ens.nc')

ds_eam_Pre_Indi = xr.open_dataset(procpath+'ds_eam_Pre_Indi_ens.nc')
ds_elm_Pre_Indi = xr.open_dataset(procpath+'ds_elm_Pre_Indi_ens.nc')
ds_mos_Pre_Indi = xr.open_dataset(procpath+'ds_mos_Pre_Indi_ens.nc')

#Loading reanalysis data
path = basepath+'/rean_JRA/'
runoff = xr.open_dataset(path+ 'JRA.h1.1996.ROF.nc')
runoff = runoff.sel(time=slice('1996-01-15', '1996-01-22'), lat=slice(30.,50.), lon=slice(265.,300.))
precip = xr.open_dataset(path+ 'JRA.h1.1996.PRECT.nc')
precip = precip.sel(time=slice('1996-01-15', '1996-01-22'), lat=slice(30.,50.), lon=slice(265.,300.))
temperature = xr.open_dataset(path+ 'JRA.h1.1996.T2M.nc')
temperature = temperature.sel(time=slice('1996-01-15', '1996-01-22'), lat=slice(30.,50.), lon=slice(265.,300.))
swe = xr.open_dataset(path+ 'JRA.h1.1996.SWE.nc')
swe = swe.sel(time=slice('1996-01-15', '1996-01-22'), lat=slice(30.,50.), lon=slice(265.,300.))

hydropath = basepath+'/rean_glofas/'
rean_dis = xr.open_dataset(hydropath+'/1996_sub.nc')
rean_dis = rean_dis.sel(time=slice('1996-01-15', '1996-01-25'),lat=slice(45.,35.),lon=slice(-85.,-70.))

#regridding data to match model grid
ds_out = xr.Dataset({'lat': (['lat'], np.arange(35, 45.1, 0.125)),
                     'lon': (['lon'], np.arange(-85, -69.875, .125)),
                    })
ds_out2 = xr.Dataset({'lat': (['lat'], np.arange(35.0625, 45, 0.125)),
                     'lon': (['lon'], np.arange(-84.9375, -70, .125)),
                    })
regridder_swe = xe.Regridder(swe, ds_out, 'bilinear')
regridder_temp = xe.Regridder(temperature, ds_out, 'bilinear')
regridder_pre = xe.Regridder(precip, ds_out, 'bilinear')
regridder_runoff = xe.Regridder(runoff, ds_out, 'bilinear')
regridder_dis = xe.Regridder(rean_dis, ds_out, 'bilinear')
regridder_lr = xe.Regridder(ds_eamLR, ds_out, 'bilinear')
regridder_lr_mos = xe.Regridder(ds_mosLR, ds_out, 'bilinear')

SWE = swe['SWE']
H2OSNO = regridder_swe(SWE)
Temp = temperature['T2M'] - 273  #change to C
TREFHT = regridder_temp(Temp)
Precip = precip['PRECT']*8.64e7 #was m/s now mm/day
PRECT = regridder_pre(Precip)
Runoff = runoff['ROF']*8.64e7   #was m/s now mm/day
QRUNOFF = regridder_runoff(Runoff)
Discharge = rean_dis['dis']
DISCHARGE = regridder_dis(Discharge)

#regridding mosart file to match elm/eam files
regridder_mos = xe.Regridder(ds_mos_Control, ds_out2, 'bilinear')
ds_mos_Control = regridder_mos(ds_mos_Control)
ds_mos_Plus_1K = regridder_mos(ds_mos_Plus_1K)
ds_mos_Plus_2K = regridder_mos(ds_mos_Plus_2K)
ds_mos_Plus_3K = regridder_mos(ds_mos_Plus_3K)
ds_mos_Plus_4K = regridder_mos(ds_mos_Plus_4K)
ds_mos_Pre_Ind = regridder_mos(ds_mos_Pre_Ind)

regridder_mosi = xe.Regridder(ds_mos_Controli, ds_out2, 'bilinear')
ds_mos_Controli = regridder_mosi(ds_mos_Controli)
ds_mos_Plus_1Ki = regridder_mosi(ds_mos_Plus_1Ki)
ds_mos_Plus_2Ki = regridder_mosi(ds_mos_Plus_2Ki)
ds_mos_Plus_3Ki = regridder_mosi(ds_mos_Plus_3Ki)
ds_mos_Plus_4Ki = regridder_mosi(ds_mos_Plus_4Ki)
ds_mos_Pre_Indi = regridder_mosi(ds_mos_Pre_Indi)

#regridding low res files
ds_eamLR = regridder_lr(ds_eamLR)
ds_elmLR = regridder_lr(ds_elmLR)
ds_mosLR = regridder_lr_mos(ds_mosLR)
#regridding control mosart file
ds_mos6 = regridder_mos(ds_mos6)

#this block creates the necessary functions to average a variable over the shapefile

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

##masks out all data that is not in the SRB shapefile
if SRB_AVG == True:
    ds_eam_Control = add_shape_coord_from_data_array(ds_eam_Control, "./shapes/srb.shp", 'coord' )
    ds_eam_Control = ds_eam_Control.where(ds_eam_Control.coord==0, other=np.nan)
    ds_elm_Control = add_shape_coord_from_data_array(ds_elm_Control, "./shapes/srb.shp", 'coord' )
    ds_elm_Control = ds_elm_Control.where(ds_elm_Control.coord==0, other=np.nan)
    ds_mos_Control = add_shape_coord_from_data_array(ds_mos_Control, "./shapes/srb.shp", 'coord' )
    ds_mos_Control = ds_mos_Control.where(ds_mos_Control.coord==0, other=np.nan)
    Control_h1 = add_shape_coord_from_data_array(Control_h1, "./shapes/srb.shp", 'coord' )
    Control_h1 = Control_h1.where(Control_h1.coord==0, other=np.nan)

    ds_eamLR = add_shape_coord_from_data_array(ds_eamLR, "./shapes/srb.shp", 'coord' )
    ds_eamLR = ds_eamLR.where(ds_eamLR.coord==0, other=np.nan)
    ds_elmLR = add_shape_coord_from_data_array(ds_elmLR, "./shapes/srb.shp", 'coord' )
    ds_elmLR= ds_elmLR.where(ds_elmLR.coord==0, other=np.nan)
    ds_mosLR = add_shape_coord_from_data_array(ds_mosLR, "./shapes/srb.shp", 'coord' )
    ds_mosLR = ds_mosLR.where(ds_mosLR.coord==0, other=np.nan)

    ds_eam6 = add_shape_coord_from_data_array(ds_eam6, "./shapes/srb.shp", 'coord' )
    ds_eam6 = ds_eam6.where(ds_eam6.coord==0, other=np.nan)
    ds_elm6 = add_shape_coord_from_data_array(ds_elm6, "./shapes/srb.shp", 'coord' )
    ds_elm6= ds_elm6.where(ds_elm6.coord==0, other=np.nan)
    ds_mos6 = add_shape_coord_from_data_array(ds_mos6, "./shapes/srb.shp", 'coord' )
    ds_mos6 = ds_mos6.where(ds_mos6.coord==0, other=np.nan)

    ds_eam_Plus_1K = add_shape_coord_from_data_array(ds_eam_Plus_1K, "./shapes/srb.shp", 'coord' )
    ds_eam_Plus_1K = ds_eam_Plus_1K.where(ds_eam_Plus_1K.coord==0, other=np.nan)
    ds_elm_Plus_1K = add_shape_coord_from_data_array(ds_elm_Plus_1K, "./shapes/srb.shp", 'coord' )
    ds_elm_Plus_1K = ds_elm_Plus_1K.where(ds_elm_Plus_1K.coord==0, other=np.nan)
    ds_mos_Plus_1K = add_shape_coord_from_data_array(ds_mos_Plus_1K, "./shapes/srb.shp", 'coord' )
    ds_mos_Plus_1K = ds_mos_Plus_1K.where(ds_mos_Plus_1K.coord==0, other=np.nan)
    Plus_1K_h1 = add_shape_coord_from_data_array(Plus_1K_h1, "./shapes/srb.shp", 'coord' )
    Plus_1K_h1 = Plus_1K_h1.where(Plus_1K_h1.coord==0, other=np.nan)

    ds_eam_Plus_2K = add_shape_coord_from_data_array(ds_eam_Plus_2K, "./shapes/srb.shp", 'coord' )
    ds_eam_Plus_2K = ds_eam_Plus_2K.where(ds_eam_Plus_2K.coord==0, other=np.nan)
    ds_elm_Plus_2K = add_shape_coord_from_data_array(ds_elm_Plus_2K, "./shapes/srb.shp", 'coord' )
    ds_elm_Plus_2K = ds_elm_Plus_2K.where(ds_elm_Plus_2K.coord==0, other=np.nan)
    ds_mos_Plus_2K = add_shape_coord_from_data_array(ds_mos_Plus_2K, "./shapes/srb.shp", 'coord' )
    ds_mos_Plus_2K = ds_mos_Plus_2K.where(ds_mos_Plus_2K.coord==0, other=np.nan)
    Plus_2K_h1 = add_shape_coord_from_data_array(Plus_2K_h1, "./shapes/srb.shp", 'coord' )
    Plus_2K_h1 = Plus_2K_h1.where(Plus_2K_h1.coord==0, other=np.nan)

    ds_eam_Plus_3K = add_shape_coord_from_data_array(ds_eam_Plus_3K, "./shapes/srb.shp", 'coord' )
    ds_eam_Plus_3K = ds_eam_Plus_3K.where(ds_eam_Plus_3K.coord==0, other=np.nan)
    ds_elm_Plus_3K = add_shape_coord_from_data_array(ds_elm_Plus_3K, "./shapes/srb.shp", 'coord' )
    ds_elm_Plus_3K = ds_elm_Plus_3K.where(ds_elm_Plus_3K.coord==0, other=np.nan)
    ds_mos_Plus_3K = add_shape_coord_from_data_array(ds_mos_Plus_3K, "./shapes/srb.shp", 'coord' )
    ds_mos_Plus_3K = ds_mos_Plus_3K.where(ds_mos_Plus_3K.coord==0, other=np.nan)
    Plus_3K_h1 = add_shape_coord_from_data_array(Plus_3K_h1, "./shapes/srb.shp", 'coord' )
    Plus_3K_h1 = Plus_3K_h1.where(Plus_3K_h1.coord==0, other=np.nan)

    ds_eam_Plus_4K = add_shape_coord_from_data_array(ds_eam_Plus_4K, "./shapes/srb.shp", 'coord' )
    ds_eam_Plus_4K = ds_eam_Plus_4K.where(ds_eam_Plus_4K.coord==0, other=np.nan)
    ds_elm_Plus_4K = add_shape_coord_from_data_array(ds_elm_Plus_4K, "./shapes/srb.shp", 'coord' )
    ds_elm_Plus_4K = ds_elm_Plus_4K.where(ds_elm_Plus_4K.coord==0, other=np.nan)
    ds_mos_Plus_4K = add_shape_coord_from_data_array(ds_mos_Plus_4K, "./shapes/srb.shp", 'coord' )
    ds_mos_Plus_4K = ds_mos_Plus_4K.where(ds_mos_Plus_4K.coord==0, other=np.nan)
    Plus_4K_h1 = add_shape_coord_from_data_array(Plus_4K_h1, "./shapes/srb.shp", 'coord' )
    Plus_4K_h1 = Plus_4K_h1.where(Plus_4K_h1.coord==0, other=np.nan)

    ds_eam_Pre_Ind = add_shape_coord_from_data_array(ds_eam_Pre_Ind, "./shapes/srb.shp", 'coord' )
    ds_eam_Pre_Ind = ds_eam_Pre_Ind.where(ds_eam_Pre_Ind.coord==0, other=np.nan)
    ds_elm_Pre_Ind = add_shape_coord_from_data_array(ds_elm_Pre_Ind, "./shapes/srb.shp", 'coord' )
    ds_elm_Pre_Ind = ds_elm_Pre_Ind.where(ds_elm_Pre_Ind.coord==0, other=np.nan)
    ds_mos_Pre_Ind = add_shape_coord_from_data_array(ds_mos_Pre_Ind, "./shapes/srb.shp", 'coord' )
    ds_mos_Pre_Ind = ds_mos_Pre_Ind.where(ds_mos_Pre_Ind.coord==0, other=np.nan)
    Pre_Ind_h1 = add_shape_coord_from_data_array(Pre_Ind_h1, "./shapes/srb.shp", 'coord' )
    Pre_Ind_h1 = Pre_Ind_h1.where(Pre_Ind_h1.coord==0, other=np.nan)

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

    QRUNOFF = add_shape_coord_from_data_array(QRUNOFF, "./shapes/srb.shp", 'coord' )
    QRUNOFF = QRUNOFF.where(QRUNOFF.coord==0, other=np.nan)

    PRECT = add_shape_coord_from_data_array(PRECT, "./shapes/srb.shp", 'coord' )
    PRECT = PRECT.where(PRECT.coord==0, other=np.nan)

    TREFHT = add_shape_coord_from_data_array(TREFHT, "./shapes/srb.shp", 'coord' )
    TREFHT = TREFHT.where(TREFHT.coord==0, other=np.nan)

    H2OSNO = add_shape_coord_from_data_array(H2OSNO, "./shapes/srb.shp", 'coord' )
    H2OSNO = H2OSNO.where(H2OSNO.coord==0, other=np.nan)

    DISCHARGE = add_shape_coord_from_data_array(DISCHARGE, "./shapes/srb.shp", 'coord' )
    DISCHARGE = DISCHARGE.where(DISCHARGE.coord==0, other=np.nan)

#function that creates an array of time values given a start date, end date, and timestep
def datetime_range(start, end, delta):
    current = start
    if not isinstance(delta, timedelta):
        delta = timedelta(**delta)
    while current < end:
        yield current
        current += delta

#creating time arrays of datetime objects
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

#creating arrays for the xticks
start = datetime.date(1996,1,15)
end = datetime.date(1996,1,25)
time = []
for dt in datetime_range(start, end, {'days':1}):
    time.append(dt)

start2 = datetime.datetime(1996,1,15,12)
end2 = datetime.datetime(1996,1,24)
time2 = []
for dt in datetime_range(start2, end2, {'hours':12}):
    time2.append(dt)

start3 = datetime.date(1996,1,16)
end3 = datetime.date(1996,1,24)
time3 = []
for dt in datetime_range(start3, end3, {'days':1}):
    time3.append('')
    time3.append(dt.strftime('%b %d'))
time3.append('')

#creating different arrays if the timing shift is set to True
if shift_ELM == True:
    T_elm = []
    start_elm = datetime.datetime(1996,1,15,9)
    end_elm = datetime.datetime(1996,1,22,21)
    Time_elm = []
    for dt in datetime_range(start_elm, end_elm, {'hours':6}):
        Time_elm.append(dt)
    for i in range(len(Time_elm)):
        T_elm.append(Time_elm[i].strftime('%Y-%m-%d-%H'))
if shift_MOSART == True:
    T_mos = []
    start_mos = datetime.datetime(1996,1,16,0)
    end_mos = datetime.datetime(1996,1,24,12)
    Time_mos = []
    for dt in datetime_range(start_mos, end_mos, {'hours':12}):
        Time_mos.append(dt)
    for i in range(len(Time_mos)):
        T_mos.append(Time_mos[i].strftime('%Y-%m-%d-%H'))

#Changing the time data in the files if the shift is set to True
elm_runs = [ds_elm_Pre_Ind, ds_elm_Control, ds_elm_Plus_1K, ds_elm_Plus_2K, ds_elm_Plus_3K,ds_elm_Plus_4K]
mos_runs = [ds_mos_Pre_Ind, ds_mos_Control, ds_mos_Plus_1K, ds_mos_Plus_2K, ds_mos_Plus_3K,ds_mos_Plus_4K]
ELM_runs = [ds_elm_Pre_Indi, ds_elm_Controli, ds_elm_Plus_1Ki, ds_elm_Plus_2Ki, ds_elm_Plus_3Ki,ds_elm_Plus_4Ki]
MOS_runs = [ds_mos_Pre_Indi, ds_mos_Controli, ds_mos_Plus_1Ki, ds_mos_Plus_2Ki, ds_mos_Plus_3Ki,ds_mos_Plus_4Ki]
Elm_runs = [ds_elm6,ds_elmLR]
Mos_runs = [ds_mos6,ds_mosLR]
if shift_ELM == True:
    for x in ELM_runs+elm_runs+Elm_runs:
        x['time'] = T_elm
if shift_MOSART == True:
    for x in MOS_runs+mos_runs+Mos_runs:
        x['time'] = T_mos

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
        g.bottom_labels = False
        g.left_labels = False
        g.xlines = False
        g.ylines = False


# ## Figure 12
##Plotting difference between counterfactual simulations and Control for river discharge

#setting up figure
fig, axs = plt.subplots(6,4,subplot_kw={'projection': ccrs.PlateCarree()},figsize=(23,24),constrained_layout=True)

#picking time indices to plot
times = [6,8,10,12]
#creating own colorbar
colors = ['w','#fde725','#a0da39','#4ac16d','#1fa187','#277f8e','#365c8d','#46327e','#440154' ]
#creating array for ylabels
runs = ['Pre-Ind', 'Control','+1K', '+2K', '+3K','+4K']

#looping through simulations
for i,r,title in zip(range(6),[ds_mos_Pre_Ind,ds_mos_Control,ds_mos_Plus_1K,ds_mos_Plus_2K,ds_mos_Plus_3K,ds_mos_Plus_4K],runs):
    for j,t in zip(range(4),times):
        #creating box for subplot letter label
        props = dict(boxstyle='square', facecolor='white', edgecolor='k')
        axs[i][j].text(.02,.97,letters[j + 4*i], fontsize='xx-large',transform=axs[i][j].transAxes,ha='left',va='top',bbox=props)
        #plotting the control in a different colorbar
        if title == 'Control':
            cont=[0,250,500,1000,2000,4000,7000,12000,18000]
            riv_dis_plot = axs[i][j].contourf(ds_mos_Control['lon'],ds_mos_Control['lat'],r['RIVER_DISCHARGE_OVER_LAND_LIQ'][t,:,:]
                                    ,levels=cont, colors=colors, extend='max')
        #plotting all other simulations as a difference plot
        else:
            cont = [-8000,-4000,-2000,-750,-250,-50,50,250,750,2000,4000,8000]
            riv_dis_plot = axs[i][j].contourf(ds_mos_Control['lon'],ds_mos_Control['lat'],r['RIVER_DISCHARGE_OVER_LAND_LIQ'][t,:,:] - ds_mos_Control['RIVER_DISCHARGE_OVER_LAND_LIQ'][t,:,:]
                                        ,levels=cont, cmap='seismic', extend='both')

        #adding extra to plot
        susq_outline(axs[i][j])
        axs[i][j].set_ylabel('+1',fontsize='x-large')
        axs[i][j].set_extent([-79.2,-74.5,39.4,43.1])
        #only adding time above first row
        if i == 0:
            axs[i][j].text(.5,1.1,Time_mos[t].strftime('%Y-%m-%d %HZ'),ha='center',va='center',fontsize=26,fontweight='bold',
                           transform=axs[i][j].transAxes)
        #only adding lat and lon labels to last row and first column
        if i == 5:
            g = axs[i][j].gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
            g.top_labels = False
            g.right_labels = False
            g.left_labels = False
            g.xlines = False
            g.ylines = False
            g.xlabel_style = {'size': 15}
        if j == 0:
            g = axs[i][j].gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
            g.top_labels = False
            g.right_labels = False
            g.bottom_labels = False
            g.xlines = False
            g.ylines = False
            g.ylabel_style = {'size': 15}

    #adding text to left of figures
    if title == 'Control':
        axs[i][0].text(-.35,.5,title,transform=axs[i][0].transAxes,fontsize=27,fontweight='bold',ha='center',va='center')
    else:
        axs[i][0].text(-.35,.5,title + ' -\nControl',transform=axs[i][0].transAxes,fontsize=27,fontweight='bold',ha='center',va='center',ma='center')

    cb = plt.colorbar(riv_dis_plot, ax=axs[i][3], aspect = 15,location='right')
    cb.ax.tick_params(labelsize='x-large')
    cb.set_label('m${^3}$ s$^{-1}$', fontsize='xx-large')

fig.savefig('map_streamflow_warming_diff.pdf', bbox_inches='tight', pad_inches=0)


#Finding the time of the maximum and the maximum value of each variable for each simulation

#choosing variables
vari = ['TREFHT','PRECT','H2OSNO', 'QRUNOFF', 'RIVER_DISCHARGE_OVER_LAND_LIQ']
titl= ['Temp','Precip','Snow','Runoff', 'Discharge']
#array of simulations
elm = [ds_elm_Pre_Ind, ds_elm_Control, ds_elm_Plus_1K, ds_elm_Plus_2K, ds_elm_Plus_3K, ds_elm_Plus_4K]
eam = [ds_eam_Pre_Ind, ds_eam_Control, ds_eam_Plus_1K, ds_eam_Plus_2K, ds_eam_Plus_3K, ds_eam_Plus_4K]
mos = [ds_mos_Pre_Ind, ds_mos_Control, ds_mos_Plus_1K, ds_mos_Plus_2K, ds_mos_Plus_3K, ds_mos_Plus_4K]
titl2 = ['Pre Ind', 'Control', 'Plus 1K', 'Plus 2K', 'Plus 3K', 'Plus 4K']

#empty arrays
max_time = []
max_s = []
min_s = []
max_r = []
max_d = []
max_p = []
max_t = []
del_swe = []

#change these depending on what snow variable you want to plot: Percent used in paper
SNOMAX = False
PERCENT = True

for v,t in zip(vari,titl):
    #looing at elm variables
    if v == 'H2OSNO' or v == 'QRUNOFF':
        #looping through elm simulations
        for d,tt in zip(elm, titl2):
            if SRB_AVG == False:
                mask = (
                (ds_elm_Control.coords["lat"] >= minlat)
                & (ds_elm_Control.coords["lat"] <= maxlat)
                & (ds_elm_Control.coords["lon"] >= minlon)
                & (ds_elm_Control.coords["lon"] <= maxlon))

                V = d[v].where(mask)
                V = V.mean( dim=("lat","lon"), skipna=True )
            else:
                #avergaing over basin before calculating time
                V = d[v].mean( dim=("lat","lon"), skipna=True )

            #finding max/min in variable
            vmax = np.nanmax(V.values)
            vmin = np.nanmin(V.values)
            #finding time of max
            vmax_ll = np.where(V == vmax)
            vmax_ll = V[vmax_ll[0]]
            vmax_time = vmax_ll.coords['time'].values

            if v == 'H2OSNO':
                #finding time where there is still 75% snow (25% snowmelt)
                if PERCENT == True:
                    threshold = vmax*0.75
                    index = []
                    for ii in range(len(V)):
                        if V[ii] < threshold:
                            index.append(ii)
                    v_perct = V[index[0]]
                    v_perct_time = v_perct.coords['time'].values
                    v_perct_time = v_perct_time.tolist()
                    max_time.append(v_perct_time)
                    max_s.append(v_perct)

                #finding time of the maximum snow amount
                if SNOMAX == True:
                    max_s.append(vmax)
                    max_time.append(vmax_time[0])
                min_s.append(vmin)
                #calculating delta SWE
                del_swe.append(vmax-vmin)

            #appending values for runoff
            if v == 'QRUNOFF':
                max_r.append(vmax)
                max_time.append(vmax_time[0])

    #looking at eam variables
    if v == 'PRECT' or v == 'TREFHT':
        #looping through eam simulations
        for d,tt in zip(eam, titl2):
            if SRB_AVG == False:
                mask = (
                (ds_eam_Control.coords["lat"] >= minlat)
                & (ds_eam_Control.coords["lat"] <= maxlat)
                & (ds_eam_Control.coords["lon"] >= minlon)
                & (ds_eam_Control.coords["lon"] <= maxlon))

                V = d[v].where(mask)
                V = V.mean( dim=("lat","lon"), skipna=True )
            else:
                #averaging over basin first
                V = d[v].mean( dim=("lat","lon"), skipna=True )

            #calculating max value
            vmax = np.nanmax(V.values)
            #finding time of max value
            vmax_ll = np.where(V == vmax)
            vmax_ll = V[vmax_ll[0]]
            vmax_time = vmax_ll.coords['time'].values
            max_time.append(vmax_time[0])
            if v == 'PRECT':
                max_p.append(vmax)
            if v == 'TREFHT':
                max_t.append(vmax)

    #looking at mosart variable
    if v == 'RIVER_DISCHARGE_OVER_LAND_LIQ':
        #looping through mosart simulations
        for d,tt in zip(mos, titl2):
            if SRB_AVG == False:
                minlon = 281.5
                maxlon = 284
                mask = (
                (ds_mos_Control.coords["lat"] >= minlat)
                & (ds_mos_Control.coords["lat"] <= maxlat)
                & (ds_mos_Control.coords["lon"] >= minlon)
                & (ds_mos_Control.coords["lon"] <= maxlon))

                V = d[v].where(mask)
                V = V.mean( dim=("lat","lon"), skipna=True )
            else:
                #averaging over basin first
                V = d[v].mean( dim=("lat","lon"), skipna=True )

            #calc max value
            vmax = np.nanmax(V.values)
            #finding time of max value
            vmax_ll = np.where(V == vmax)
            vmax_ll = V[vmax_ll[0]]
            vmax_time = vmax_ll.coords['time'].values
            max_time.append(vmax_time[0])
            max_d.append(vmax)

#doing same as above cell but calculating a time for each of the 15 enesemble members rather than as an esmeble mean like above

minlat = 40
maxlat = 42.5
minlon = -78.5
maxlon = -76
vari = ['TREFHT','PRECT','H2OSNO', 'QRUNOFF', 'RIVER_DISCHARGE_OVER_LAND_LIQ']
titl= ['Temp','Precip','Snow','Runoff', 'Discharge']
elm = [ds_elm_Pre_Indi, ds_elm_Controli, ds_elm_Plus_1Ki, ds_elm_Plus_2Ki, ds_elm_Plus_3Ki, ds_elm_Plus_4Ki]
eam = [ds_eam_Pre_Indi, ds_eam_Controli, ds_eam_Plus_1Ki, ds_eam_Plus_2Ki, ds_eam_Plus_3Ki, ds_eam_Plus_4Ki]
mos = [ds_mos_Pre_Indi, ds_mos_Controli, ds_mos_Plus_1Ki, ds_mos_Plus_2Ki, ds_mos_Plus_3Ki, ds_mos_Plus_4Ki]
titl2 = ['Pre Ind', 'Control', 'Plus 1K', 'Plus 2K', 'Plus 3K', 'Plus 4K']
max_time3 = []
max_si = []
max_ri = []
max_di = []
max_pi = []
max_ti = []

for v,t in zip(vari,titl):
    if v == 'H2OSNO' or v == 'QRUNOFF':
        for d,tt in zip(elm, titl2):
            for i in range(15):
                if SRB_AVG == False:
                    mask = (
                    (ds_elm_Controli.coords["lat"] >= minlat)
                    & (ds_elm_Controli.coords["lat"] <= maxlat)
                    & (ds_elm_Controli.coords["lon"] >= minlon)
                    & (ds_elm_Controli.coords["lon"] <= maxlon))

                    V = d[v].where(mask)
                    V = V.mean( dim=("lat","lon"), skipna=True )
                    V = V[i,:]
                else:
                    V = d[v].mean( dim=("lat","lon"), skipna=True )
                    V = V[i,:]

                vmax = np.nanmax(V.values)
                vmax_ll = np.where(V == vmax)
                vmax_ll = V[vmax_ll[0]]
                vmax_time = vmax_ll.coords['time'].values
                if v == 'H2OSNO':
                    if PERCENT == True:
                        threshold = vmax*0.75
                        index = []
                        for ii in range(len(V)):
                            if V[ii] < threshold:
                                index.append(ii)
                        v_perct = V[index[0]]
                        v_perct_time = v_perct.coords['time'].values
                        v_perct_time = v_perct_time.tolist()
                        max_time3.append(v_perct_time)
                        max_si.append(v_perct)
                    if SNOMAX == True:
                        max_si.append(vmax)
                        max_time3.append(vmax_time[0])
                if v == 'QRUNOFF':
                    v_m_t = datetime.datetime.strptime(vmax_time[0], '%Y-%m-%d-%H')
                    threshold = datetime.datetime(1996, 1, 21, 0, 0)
                    if v_m_t > threshold:
                        x = V.values
                        X = np.delete(x,np.where(V == vmax)[0][0])
                        vmax2 = np.nanmax(X)
                        vmax_2 = np.where(V == vmax2)
                        vmax_2 = V[vmax_2[0]]
                        vmax_time = vmax_2.coords['time'].values
                    max_time3.append(vmax_time[0])
                    max_ri.append(vmax)

    if v == 'PRECT' or v == 'TREFHT':
        for d,tt in zip(eam, titl2):
            for i in range(15):
                if SRB_AVG == False:
                    mask = (
                    (ds_eam_Controli.coords["lat"] >= minlat)
                    & (ds_eam_Controli.coords["lat"] <= maxlat)
                    & (ds_eam_Controli.coords["lon"] >= minlon)
                    & (ds_eam_Controli.coords["lon"] <= maxlon))

                    V = d[v].where(mask)
                    V = V.mean( dim=("lat","lon"), skipna=True )
                    V = V[i,:]
                else:
                    V = d[v].mean( dim=("lat","lon"), skipna=True )
                    V = V[i,:]

                vmax = np.nanmax(V.values)
                vmax_ll = np.where(V == vmax)
                vmax_ll = V[vmax_ll[0]]
                vmax_time = vmax_ll.coords['time'].values
                max_time3.append(vmax_time[0])
                if v == 'PRECT':
                    max_pi.append(vmax)
                if v == 'TREFHT':
                    max_ti.append(vmax)

    if v == 'RIVER_DISCHARGE_OVER_LAND_LIQ':
        for d,tt in zip(mos, titl2):
            for i in range(15):
                if SRB_AVG == False:
                    minlon = 281.5
                    maxlon = 284
                    mask = (
                    (ds_mos_Controli.coords["lat"] >= minlat)
                    & (ds_mos_Controli.coords["lat"] <= maxlat)
                    & (ds_mos_Controli.coords["lon"] >= minlon)
                    & (ds_mos_Controli.coords["lon"] <= maxlon))

                    V = d[v].where(mask)
                    V = V.mean( dim=("lat","lon"), skipna=True )
                    V = V[i,:]
                else:
                    V = d[v].mean( dim=("lat","lon"), skipna=True )
                    V = V[i,:]

                vmax = np.nanmax(V.values)
                vmax_ll = np.where(V == vmax)
                vmax_ll = V[vmax_ll[0]]
                vmax_time = vmax_ll.coords['time'].values
                max_time3.append(vmax_time[0])
                max_di.append(vmax)

#creating function to create csvs
from csv import writer
def append_list_as_row(file_name, list_of_elem):
    # Open file in append mode
    with open(file_name, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow(list_of_elem)

#creating csv of time of the maximums
f = open('Time_of_Maximums.csv','w')
vari = ['vari','Pre Ind', 'Control', 'Plus 1K', 'Plus 2K', 'Plus 3K', 'Plus 4K']
append_list_as_row('Time_of_Maximums.csv',vari)

b_ind = np.arange(0,len(max_time),7)
for x in b_ind:
    append_list_as_row('Time_of_Maximums.csv', max_time[x:x+7])


#turning the max time into a datetime object for ensemble mean
max_date = []
for i in range(len(max_time)):
    max_date.append(datetime.datetime.strptime(max_time[i], '%Y-%m-%d-%H'))

#turning the max time into a datetime object for ensemble members
max_datei = []
for i in range(len(max_time3)):
    max_datei.append(datetime.datetime.strptime(max_time3[i], '%Y-%m-%d-%H'))

#calculating the average max time for the individual ensemble members to get a different version of ensemble mean max time
avgTime = []
for x in np.arange(0,450,15):
    avgTime.append(datetime.datetime.strftime(datetime.datetime.fromtimestamp(sum(map(datetime.datetime.timestamp,max_datei[x:x+15]))/len(max_datei[x:x+15])),'%Y-%m-%d-%H'))

#turning above into datetime object
avgTimei = []
for j in range(len(avgTime)):
    avgTimei.append(datetime.datetime.strptime(avgTime[j], '%Y-%m-%d-%H'))


#creating array to set xticks
start = datetime.datetime(1996,1,17)
end = datetime.datetime(1996,1,24)
time_plot2 = []
for dt in datetime_range(start, end, {'hours':6}):
    time_plot2.append(dt)

#creating array for xtick labels
time_s = []
for x in np.arange(0,len(time_plot2),4):
    time_s.append(time_plot2[x].strftime('%b %d'))
    for y in range(3):
        time_s.append('')


# ## Figure 13/15
##Plots time of maximum for each variable to see trend in timing of variables as the simulations warm

#Set AVG_2 = True if want to make figure 17 of thesis (plots both averaging techniques) set = False for figure 13
#CMZ added ii loop so both AVG_2 = True and False upon single call of script...

for ii in range(0, 2):

  if ii == 0:
    AVG_2 = False
  else:
    AVG_2 = True

  print(ii, AVG_2)

  #creating array for y axis values
  temps = [0, .3 ,1, 2, 3, 4]
  temps2 = temps*5
  temps3 = np.repeat(temps,15)
  #creating labels
  ytik = ['PI', 'Control', '+1K', '+2K', '+3K', '+4K']
  vari = ['Maximum Temperature','Maximum Precipitation','25% Snow melt ', 'Maximum Runoff', 'Maximum River Discharge']
  #picking colors for plots
  colors = ['r','b','g','k','m']
  #setting up plot
  fig, axs = plt.subplots(5,1,figsize=(9,13))
  for x,y,c,v,ax,L in zip(np.arange(0,len(max_datei),90),np.arange(0,30,6), colors,vari,axs,letters):
      #plotting ensemble members
      ax.scatter(max_datei[x:x+90], temps3, color=c, alpha = 0.1)
      #plotting ensemble mean
      p1 = ax.plot(max_date[y:y+6], temps2[y:y+6], color=c,label='Average 1')

      #plotting ensemble mean using average 2 technique
      if AVG_2 == True:
          p2 = ax.plot(avgTimei[y:y+6], temps2[y:y+6], color=c, linestyle='dashed', label='Average 2')
      ax.margins(x=0)

      #putting the subplot title in a different location for discharge than the other variables
      if v == 'Maximum River Discharge':
          ax.text(.57,.95,L+' Time of '+v, fontsize='x-large',transform=ax.transAxes,ha='right',va='top')
          ax.set_xticks(time_plot2[:-2])
          ax.set_xticklabels(time_s[:-2],fontsize='large')
      #placing subplot title
      else:
          ax.text(.99,.97,L+' Time of '+v, fontsize='x-large',transform=ax.transAxes,ha='right',va='top')
          ax.set_xticks(time_plot2[:-2])
          ax.set_xticklabels([])

      #adding to plot
      ax.set_yticks(temps)
      ax.set_yticklabels(ytik, fontsize='medium')
      ax.grid(color='gainsboro',linewidth=.5)
      #only adding legend to first plot if plotting both averages
      if AVG_2 == True:
          if v == 'Maximum Temperature' :
              props = dict(boxstyle='square', facecolor='white', edgecolor='grey')
              ax.text(.015,.95,'Solid Line: Average 1 \nDashed Line: Average 2', fontsize='large',transform=ax.transAxes,ha='left',va='top',bbox=props)
      #making plots close togther
      plt.subplots_adjust(hspace=10**-6)

  if ii == 0:
    fig.savefig('timing_shifts_line.pdf', bbox_inches='tight', pad_inches=0)
  else:
    fig.savefig('orig_avg_comp_timing.pdf', bbox_inches='tight', pad_inches=0)


# ## Figure 11

##plots 4 variables spatially for the different simulations

#setting up figure
fig, axs = plt.subplots(6,4,subplot_kw={'projection': ccrs.PlateCarree()},figsize=(17,21),constrained_layout=True)
#creating titles
runs = ['Pre-Ind', 'Control','+1K', '+2K', '+3K','+4K']
var = ['Max Temperature','Accumulated Precip.','Delta SWE', 'Max Runoff']

#loop that goes through the 6 simulations and plots the variable in question for each sim
for i,l,a,h,title,L in zip(range(6),[ds_elm_Pre_Ind,ds_elm_Control,ds_elm_Plus_1K,ds_elm_Plus_2K,ds_elm_Plus_3K,ds_elm_Plus_4K],
                     [ds_eam_Pre_Ind,ds_eam_Control,ds_eam_Plus_1K,ds_eam_Plus_2K,ds_eam_Plus_3K,ds_eam_Plus_4K],
                       [Pre_Ind_h1,Control_h1, Plus_1K_h1,Plus_2K_h1,Plus_3K_h1,Plus_4K_h1 ],runs,letters):

    #calculating max runoff and max temp for whole time period
    run_max = np.amax(l['QRUNOFF'], axis=0)
    temp_max = np.amax(a['TREFHT'], axis=0)
    #calculating accumulated precip
    pre_acc = np.sum(h['PRECT']*(3/24), axis=0)   #changing mm/day to mm by multiplying each rate by the time interval (3 hrs/)
    #calculating delta swe
    swe = np.amax(l['H2OSNO'], axis=0) - np.amin(l['H2OSNO'], axis=0)

    #setting contour values
    cont_pacc = np.arange(25,81,5)
    cont_tmax = np.arange(0,21,2)
    cont_swe = np.arange(40,121,10)
    cont_maxr = np.arange(0,421,30)
    #plotting variables
    pa_plot = axs[i][1].contourf(ds_elm_Control['lon'],ds_elm_Control['lat'],pre_acc, levels=cont_pacc,cmap = 'viridis_r')
    tm_plot = axs[i][0].contourf(ds_elm_Control['lon'],ds_elm_Control['lat'],temp_max,levels = cont_tmax,cmap = 'YlOrRd')
    swe_plot = axs[i][2].contourf(ds_elm_Control['lon'],ds_elm_Control['lat'],swe, levels=cont_swe,cmap = 'BuPu')
    rm_plot = axs[i][3].contourf(ds_elm_Control['lon'],ds_elm_Control['lat'],run_max,levels=cont_maxr, cmap = 'inferno_r')#'GnBu')
    #adding text to left hand side of each row
    axs[i][0].text(-.4,.5,title,transform=axs[i][0].transAxes,fontsize=20,fontweight='bold',ha='center',va='center')
    #looping through subplots and adding letter to each subplot
    for j,v in zip(range(4),var):
        props = dict(boxstyle='square', facecolor='white', edgecolor='k')
        axs[i][j].text(.02,.97,letters[j + 4*i], fontsize='xx-large',transform=axs[i][j].transAxes,ha='left',va='top',bbox=props)
        #adding a title to column for the first row only
        if i == 0:
            axs[i][j].text(.5,1.18,v,ha='center',va='center',fontsize=20,fontweight='bold',transform=axs[i][j].transAxes)
        #adding susquehanna river basin outline
        susq_outline(axs[i][j])
        axs[i][j].set_extent([-79.2,-74.5,39.4,43.1])

    #adding longitude labels to last row
    if i == 5:
        for j in range(4):
            g = axs[i][j].gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
            g.top_labels = False
            g.right_labels = False
            g.left_labels = False
            g.xlines = False
            g.ylines = False
            g.xlabel_style = {'size': 13}
    #adding latitude labels for first column
    gl = axs[i][0].gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
    gl.xlines = False
    gl.ylines = False
    gl.ylabel_style = {'size': 13}
#adding colorbars to bottom of each column
cb = plt.colorbar(pa_plot, ax=[axs[5,1]], shrink=.95,aspect=17,location='bottom')
cb.ax.tick_params(labelsize='x-large')
cb.ax.set_title('mm', fontsize='xx-large')
cb1 = plt.colorbar(tm_plot,  ax=[axs[5, 0]],shrink=.95,aspect=17,location='bottom')
cb1.ax.tick_params(labelsize='x-large')
cb1.ax.set_title('$^\circ$C', fontsize='xx-large')
cb2 = plt.colorbar(swe_plot, ax=[axs[5,2]],shrink=.95,aspect=17,location='bottom')
cb2.ax.tick_params(labelsize='x-large')
cb2.ax.set_title('mm', fontsize='xx-large')
cb3 = plt.colorbar(rm_plot, ax=[axs[5,3]],shrink=.95,aspect=17,location='bottom')
cb3.ax.tick_params(labelsize='x-large')
cb3.ax.set_title('mm day$^{-1}$', fontsize='xx-large');

fig.savefig('spatial_maps_bulk_quants.pdf', bbox_inches='tight', pad_inches=0)


# ## Figure 18

##plots initial SWE values across basin for different simulations

#setting up figure
fig, axs = plt.subplots(1,6,subplot_kw={'projection': ccrs.PlateCarree()},figsize=(20,4),constrained_layout=True)
#creating labels
runs = ['Pre-Ind', 'Control','+1K', '+2K', '+3K','+4K']
var = ['Max Temperature']
#looping through elm simulations
for ax,l,title,L in zip(axs,[ds_elm_Pre_Ind,ds_elm_Control,ds_elm_Plus_1K,ds_elm_Plus_2K,ds_elm_Plus_3K,ds_elm_Plus_4K],
                    runs,letters):
    #calculating max SWE value
    SWE = np.amax(l['H2OSNO'], axis=0)
    #establishing contours
    cont_SWE = np.arange(0,161,20)
    #plotting SWE
    SWE_plot = ax.contourf(ds_elm_Control['lon'],ds_elm_Control['lat'],SWE,levels=cont_SWE,cmap = 'plasma_r')
    #adding run title above each box
    props = dict(boxstyle='square', facecolor='white', edgecolor='k')
    ax.text(.5,1.12,title,ha='center',va='center',fontsize='xx-large',fontweight='bold',transform=ax.transAxes)
    #adding SRB outline
    susq_outline(ax)
    ax.set_extent([-79.2,-74.5,39.4,43.1])
    #adding variable title to left of subplots
    if title == 'Pre-Ind':
        ax.text(-.3,.5,'Max SWE',transform=ax.transAxes,fontsize='xx-large',fontweight='bold',ha='center',va='center')
    #adding colorbar to right of figure
    if title == '+4K':
        cb1 = plt.colorbar(SWE_plot,  ax=ax,shrink=.65,aspect=60, location='right')
        cb1.ax.tick_params(labelsize='x-large')
        cb1.set_label('mm', fontsize='xx-large')

fig.savefig('spatial_maps_maxSWE.pdf', bbox_inches='tight', pad_inches=0)


# ## Figure 17

##Plotting comparison of Low res, high res, and reanalysis

#creating own colorbar
colors = ['w','#fde725','#a0da39','#4ac16d','#1fa187','#277f8e','#365c8d','#46327e','#440154' ]
#setting up plot
fig, axs = plt.subplots(2,3,subplot_kw={'projection': ccrs.PlateCarree()},figsize=(12,5),constrained_layout=True)
#looping through two variables and plotting them soatially for the 3 different sims
for i,v,vari in zip(range(2),['QRUNOFF','RIVER_DISCHARGE_OVER_LAND_LIQ'],['Runoff','River \nDischarge']):
    if v == 'RIVER_DISCHARGE_OVER_LAND_LIQ':
        ds = ds_mos6
        dsLR = ds_mosLR
        rean = DISCHARGE
    else:
        ds = ds_elm6
        dsLR = ds_elmLR
        rean = QRUNOFF
    #looping through the 3 different sims for the variable
    for j,r,t in zip(range(3), [ds,dsLR,rean], ['14km Resolution','110km Resolution','Reanalysis']):
        #adding letters to each subplot
        props = dict(boxstyle='square', facecolor='white', edgecolor='k')
        axs[i][j].text(.02,.97,letters[j + 3*i], fontsize='x-large',transform=axs[i][j].transAxes,ha='left',va='top',bbox=props)
        #adding text above the first row
        if i == 0:
            axs[i][j].text(.5,1.1,t,ha='center',va='center',fontsize=16,fontweight='bold',
                           transform=axs[i][j].transAxes)
        #adding text to left of the two rows
        axs[i][0].text(-.45,.5,vari,transform=axs[i][0].transAxes,fontsize=16,fontweight='bold',ha='center',va='center')
        #plotting discharge
        if v == 'RIVER_DISCHARGE_OVER_LAND_LIQ':
            cont=[0,250,500,750,1500,3000,6000,11000,21000]
            #plotting with two different time indicies bc reanalysis is different than model simulations
            if t == 'Reanalysis':
                plot = axs[i][j].contourf(r['lon'],r['lat'],r[7,:,:],levels=cont, colors = colors)
            else:
                plot = axs[i][j].contourf(r['lon'],r['lat'],r[v][12,:,:],levels=cont, colors = colors)
            susq_outline(axs[i][j])
            axs[i][j].set_extent([-79.2,-74.5,39.5,43])
        #plotting runoff
        else:
            cont = [0,50,100,150,200,250,300,350]

            #plotting with two different time indicies bc reanalysis is different than model simulations
            if t == 'Reanalysis':
                plot = axs[i][j].contourf(ds['lon'],ds['lat'],r[19,:,:],levels=cont, cmap = 'inferno_r')
            else:
                plot = axs[i][j].contourf(ds['lon'],ds['lat'],r[v][17,:,:],levels=cont, cmap = 'inferno_r')
            susq_outline(axs[i][j])
            axs[i][j].set_extent([-79.2,-74.5,39.4,43.1])
        #adding longitude labels to bottom row
        if i == 1:
            g = axs[i][j].gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
            g.top_labels = False
            g.right_labels = False
            g.left_labels = False
            g.xlines = False
            g.ylines = False
            g.xlabel_style = {'size': 12}
        #adding latitude labels to left column
        if j == 0:
            g = axs[i][j].gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
            g.top_labels = False
            g.right_labels = False
            g.bottom_labels = False
            g.xlines = False
            g.ylines = False
            g.ylabel_style = {'size': 12}
    #adding colorbar to right of row
    cb = plt.colorbar(plot, ax = axs[i][2], shrink=.95,aspect= 20,label = 'mm/day',location='right')
    cb.ax.tick_params(labelsize='large')
    cb.set_label('mm day$^{-1}$', fontsize='x-large')

fig.savefig('spatial_maps_resolution.pdf', bbox_inches='tight', pad_inches=0)


# ## Figure 14

##compares timeseries of the 1996 event for the high res and low res model runs

#taking average of the variable over the basin
if SRB_AVG == False:
    mask = (
    (ds_elm_Control.coords["lat"] > x)
    & (ds_elm_Control.coords["lat"] < x+1)
    & (ds_elm_Control.coords["lon"] > y)
    & (ds_elm_Control.coords["lon"] < y+1)
        )

    ds_mos = ds_mos_Control.where(mask)
    ds_elm = ds_elm_Control.where(mask)
    ds_eam = ds_eam_Control.where(mask)
    ds_eam_mean = ds_eam.mean( dim=("lat","lon"), skipna=True )
    ds_elm_mean = ds_elm.mean( dim=("lat","lon"), skipna=True )
    ds_mos_mean = ds_mos.mean( dim=("lat","lon"), skipna=True )
if SRB_AVG == True:
    ds_eam_mean = ds_eam_Control.mean( dim=("lat","lon"), skipna=True )
    ds_elm_mean = ds_elm_Control.mean( dim=("lat","lon"), skipna=True )
    ds_mos_mean = ds_mos_Control.mean( dim=("lat","lon"), skipna=True )
    cl_eam_mean = ds_eamLR.mean( dim=("lat","lon"), skipna=True )
    cl_elm_mean = ds_elmLR.mean( dim=("lat","lon"), skipna=True )
    cl_mos_mean = ds_mosLR.mean( dim=("lat","lon"), skipna=True )

    ## Setup plot
    fig = plt.figure(figsize=(14,5))
    host1 = fig.add_subplot(111)

    ## Create extra y-axes that shares x-axis with "host"
    yaxis1 = host1.twinx()
    yaxis2 = host1.twinx()

    ## Set axis limits
    host1.set_ylim([-5, 130])
    yaxis1.set_ylim([-5, 180])
    yaxis2.set_ylim([-30,20])


    # Offset the right spine of par2. and show it
    yaxis2.spines["right"].set_position(("axes", 1.1))
    yaxis2.spines["right"].set_visible(True)

    ## Set labels for axes
    host1.set_ylabel("Snow Water Equivalent (mm)", color='c', fontsize='x-large')
    yaxis1.set_ylabel("Runoff & Precipitation (mm day$^{-1}$)", color='b',fontsize='x-large')
    yaxis2.set_ylabel("2-m Temperature ($^\circ$C)", color='r',fontsize='x-large')

    ## Select colors
    color1 = 'c'
    color2 = 'b'
    color3 = 'g'
    color4 = 'r'

    ## Plot each line
    p1, = host1.plot(Time_elm, ds_elm_mean['H2OSNO'], '-', color=color1,label="HR - SWE")
    p2, = yaxis1.plot(Time_elm, ds_elm_mean['QRUNOFF'], '-', color=color2, label="HR - Runoff")
    p3, = yaxis1.plot(Time_eam, ds_eam_mean['PRECT'], '-', color=color3, label="HR - Precip")
    p4, = yaxis2.plot(Time_eam, ds_eam_mean['TREFHT'], '-', color=color4, label="HR - Temp")

    p5, = host1.plot(Time_elm, cl_elm_mean['H2OSNO'], '--', color=color1,label="LR - SWE")
    p6, = yaxis1.plot(Time_elm, cl_elm_mean['QRUNOFF'], '--', color=color2, label="LR - Runoff")
    p7, = yaxis1.plot(Time_eam, cl_eam_mean['PRECT'], '--', color=color3, label="LR - Precip")
    p8, = yaxis2.plot(Time_eam, cl_eam_mean['TREFHT'], '--', color=color4, label="LR - Temp")

    ## Add legend ad labels
    lns = [p1, p2, p3, p4,p5,p6,p7,p8]
    host1.legend(handles=lns, loc=(0.68,0.69), ncol=2, fontsize=12)
    host1.set_xticks(time2)
    host1.set_xticklabels(time3, fontsize='large')
    host1.margins(x=0)
    host1.tick_params('both', labelsize='large')
    yaxis1.tick_params('y', labelsize='large')
    yaxis2.tick_params('y', labelsize='large')
    fig.tight_layout(pad=3.0)

fig.savefig('hydrograph_resolution.pdf', bbox_inches='tight', pad_inches=0)
