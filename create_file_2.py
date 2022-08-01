#!/usr/bin/env python
# coding: utf-8
##Importing necessary libraries
import netCDF4 as nc
import numpy as np
import xarray as xr
from datetime import datetime
from datetime import timedelta

#loading files from each simulation and run for the non averaged initialization times
ds_eam_Controli6 = xr.open_dataset('Control-6_init_eam_slice.nc')
ds_elm_Controli6 = xr.open_dataset('Control-6_init_elm_slice.nc')
ds_mos_Controli6 = xr.open_dataset('Control-6_init_mos_slice.nc')
ds_eam_Plus_1Ki6 = xr.open_dataset('Plus1K-6_init_eam_slice.nc')
ds_elm_Plus_1Ki6 = xr.open_dataset('Plus1K-6_init_elm_slice.nc')
ds_mos_Plus_1Ki6 = xr.open_dataset('Plus1K-6_init_mos_slice.nc')
ds_eam_Plus_2Ki6 = xr.open_dataset('Plus2K-6_init_eam_slice.nc')
ds_elm_Plus_2Ki6 = xr.open_dataset('Plus2K-6_init_elm_slice.nc')
ds_mos_Plus_2Ki6 = xr.open_dataset('Plus2K-6_init_mos_slice.nc')
ds_eam_Plus_3Ki6 = xr.open_dataset('Plus3K-6_init_eam_slice.nc')
ds_elm_Plus_3Ki6 = xr.open_dataset('Plus3K-6_init_elm_slice.nc')
ds_mos_Plus_3Ki6 = xr.open_dataset('Plus3K-6_init_mos_slice.nc')
ds_eam_Plus_4Ki6 = xr.open_dataset('Plus4K-6_init_eam_slice.nc')
ds_elm_Plus_4Ki6 = xr.open_dataset('Plus4K-6_init_elm_slice.nc')
ds_mos_Plus_4Ki6 = xr.open_dataset('Plus4K-6_init_mos_slice.nc')
ds_eam_Pre_Indi6 = xr.open_dataset('Pre_Ind-6_init_eam_slice.nc')
ds_elm_Pre_Indi6 = xr.open_dataset('Pre_Ind-6_init_elm_slice.nc')
ds_mos_Pre_Indi6 = xr.open_dataset('Pre_Ind-6_init_mos_slice.nc')

ds_eam_Controli4 = xr.open_dataset('Control-4_init_eam_slice.nc')
ds_elm_Controli4 = xr.open_dataset('Control-4_init_elm_slice.nc')
ds_mos_Controli4 = xr.open_dataset('Control-4_init_mos_slice.nc')
ds_eam_Plus_1Ki4 = xr.open_dataset('Plus1K-4_init_eam_slice.nc')
ds_elm_Plus_1Ki4 = xr.open_dataset('Plus1K-4_init_elm_slice.nc')
ds_mos_Plus_1Ki4 = xr.open_dataset('Plus1K-4_init_mos_slice.nc')
ds_eam_Plus_2Ki4 = xr.open_dataset('Plus2K-4_init_eam_slice.nc')
ds_elm_Plus_2Ki4 = xr.open_dataset('Plus2K-4_init_elm_slice.nc')
ds_mos_Plus_2Ki4 = xr.open_dataset('Plus2K-4_init_mos_slice.nc')
ds_eam_Plus_3Ki4 = xr.open_dataset('Plus3K-4_init_eam_slice.nc')
ds_elm_Plus_3Ki4 = xr.open_dataset('Plus3K-4_init_elm_slice.nc')
ds_mos_Plus_3Ki4 = xr.open_dataset('Plus3K-4_init_mos_slice.nc')
ds_eam_Plus_4Ki4 = xr.open_dataset('Plus4K-4_init_eam_slice.nc')
ds_elm_Plus_4Ki4 = xr.open_dataset('Plus4K-4_init_elm_slice.nc')
ds_mos_Plus_4Ki4 = xr.open_dataset('Plus4K-4_init_mos_slice.nc')
ds_eam_Pre_Indi4 = xr.open_dataset('Pre_Ind-4_init_eam_slice.nc')
ds_elm_Pre_Indi4 = xr.open_dataset('Pre_Ind-4_init_elm_slice.nc')
ds_mos_Pre_Indi4 = xr.open_dataset('Pre_Ind-4_init_mos_slice.nc')

ds_eam_Controli5 = xr.open_dataset('Control-5_init_eam_slice.nc')
ds_elm_Controli5 = xr.open_dataset('Control-5_init_elm_slice.nc')
ds_mos_Controli5 = xr.open_dataset('Control-5_init_mos_slice.nc')
ds_eam_Plus_1Ki5 = xr.open_dataset('Plus1K-5_init_eam_slice.nc')
ds_elm_Plus_1Ki5 = xr.open_dataset('Plus1K-5_init_elm_slice.nc')
ds_mos_Plus_1Ki5 = xr.open_dataset('Plus1K-5_init_mos_slice.nc')
ds_eam_Plus_2Ki5 = xr.open_dataset('Plus2K-5_init_eam_slice.nc')
ds_elm_Plus_2Ki5 = xr.open_dataset('Plus2K-5_init_elm_slice.nc')
ds_mos_Plus_2Ki5 = xr.open_dataset('Plus2K-5_init_mos_slice.nc')
ds_eam_Plus_3Ki5 = xr.open_dataset('Plus3K-5_init_eam_slice.nc')
ds_elm_Plus_3Ki5 = xr.open_dataset('Plus3K-5_init_elm_slice.nc')
ds_mos_Plus_3Ki5 = xr.open_dataset('Plus3K-5_init_mos_slice.nc')
ds_eam_Plus_4Ki5 = xr.open_dataset('Plus4K-5_init_eam_slice.nc')
ds_elm_Plus_4Ki5 = xr.open_dataset('Plus4K-5_init_elm_slice.nc')
ds_mos_Plus_4Ki5 = xr.open_dataset('Plus4K-5_init_mos_slice.nc')
ds_eam_Pre_Indi5 = xr.open_dataset('Pre_Ind-5_init_eam_slice.nc')
ds_elm_Pre_Indi5 = xr.open_dataset('Pre_Ind-5_init_elm_slice.nc')
ds_mos_Pre_Indi5 = xr.open_dataset('Pre_Ind-5_init_mos_slice.nc')


# In[ ]:


#Combining 004, 005, 006 runs for "super ensemble" for non averaged initialization times files
DS_eam_Controli = xr.concat((ds_eam_Controli6, ds_eam_Controli4, ds_eam_Controli5), dim='init')
DS_elm_Controli = xr.concat((ds_elm_Controli6, ds_elm_Controli4, ds_elm_Controli5), dim='init')
DS_mos_Controli = xr.concat((ds_mos_Controli6, ds_mos_Controli4, ds_mos_Controli5), dim='init')

DS_eam_Plus_1Ki = xr.concat((ds_eam_Plus_1Ki6, ds_eam_Plus_1Ki4, ds_eam_Plus_1Ki5), dim='init')
DS_elm_Plus_1Ki = xr.concat((ds_elm_Plus_1Ki6, ds_elm_Plus_1Ki4, ds_elm_Plus_1Ki5), dim='init')
DS_mos_Plus_1Ki = xr.concat((ds_mos_Plus_1Ki6, ds_mos_Plus_1Ki4, ds_mos_Plus_1Ki5), dim='init')

DS_eam_Plus_2Ki = xr.concat((ds_eam_Plus_2Ki6, ds_eam_Plus_2Ki4, ds_eam_Plus_2Ki5), dim='init')
DS_elm_Plus_2Ki = xr.concat((ds_elm_Plus_2Ki6, ds_elm_Plus_2Ki4, ds_elm_Plus_2Ki5), dim='init')
DS_mos_Plus_2Ki = xr.concat((ds_mos_Plus_2Ki6, ds_mos_Plus_2Ki4, ds_mos_Plus_2Ki5), dim='init')

DS_eam_Plus_3Ki = xr.concat((ds_eam_Plus_3Ki6, ds_eam_Plus_3Ki4, ds_eam_Plus_3Ki5), dim='init')
DS_elm_Plus_3Ki = xr.concat((ds_elm_Plus_3Ki6, ds_elm_Plus_3Ki4, ds_elm_Plus_3Ki5), dim='init')
DS_mos_Plus_3Ki = xr.concat((ds_mos_Plus_3Ki6, ds_mos_Plus_3Ki4, ds_mos_Plus_3Ki5), dim='init')

DS_eam_Plus_4Ki = xr.concat((ds_eam_Plus_4Ki6, ds_eam_Plus_4Ki4, ds_eam_Plus_4Ki5), dim='init')
DS_elm_Plus_4Ki = xr.concat((ds_elm_Plus_4Ki6, ds_elm_Plus_4Ki4, ds_elm_Plus_4Ki5), dim='init')
DS_mos_Plus_4Ki = xr.concat((ds_mos_Plus_4Ki6, ds_mos_Plus_4Ki4, ds_mos_Plus_4Ki5), dim='init')

DS_eam_Pre_Indi = xr.concat((ds_eam_Pre_Indi6, ds_eam_Pre_Indi4, ds_eam_Pre_Indi5), dim='init')
DS_elm_Pre_Indi = xr.concat((ds_elm_Pre_Indi6, ds_elm_Pre_Indi4, ds_elm_Pre_Indi5), dim='init')
DS_mos_Pre_Indi = xr.concat((ds_mos_Pre_Indi6, ds_mos_Pre_Indi4, ds_mos_Pre_Indi5), dim='init')


# In[ ]:


#saving ensemble as a netCDF file
DS_eam_Controli.to_netcdf('ds_eam_Controli_ens.nc')
DS_elm_Controli.to_netcdf('ds_elm_Controli_ens.nc')
DS_mos_Controli.to_netcdf('ds_mos_Controli_ens.nc')

DS_eam_Plus_1Ki.to_netcdf('ds_eam_Plus_1Ki_ens.nc')
DS_elm_Plus_1Ki.to_netcdf('ds_elm_Plus_1Ki_ens.nc')
DS_mos_Plus_1Ki.to_netcdf('ds_mos_Plus_1Ki_ens.nc')

DS_eam_Plus_2Ki.to_netcdf('ds_eam_Plus_2Ki_ens.nc')
DS_elm_Plus_2Ki.to_netcdf('ds_elm_Plus_2Ki_ens.nc')
DS_mos_Plus_2Ki.to_netcdf('ds_mos_Plus_2Ki_ens.nc')

DS_eam_Plus_3Ki.to_netcdf('ds_eam_Plus_3Ki_ens.nc')
DS_elm_Plus_3Ki.to_netcdf('ds_elm_Plus_3Ki_ens.nc')
DS_mos_Plus_3Ki.to_netcdf('ds_mos_Plus_3Ki_ens.nc')

DS_eam_Plus_4Ki.to_netcdf('ds_eam_Plus_4Ki_ens.nc')
DS_elm_Plus_4Ki.to_netcdf('ds_elm_Plus_4Ki_ens.nc')
DS_mos_Plus_4Ki.to_netcdf('ds_mos_Plus_4Ki_ens.nc')

DS_eam_Pre_Indi.to_netcdf('ds_eam_Pre_Indi_ens.nc')
DS_elm_Pre_Indi.to_netcdf('ds_elm_Pre_Indi_ens.nc')
DS_mos_Pre_Indi.to_netcdf('ds_mos_Pre_Indi_ens.nc')


# In[ ]:


#loading files from each simulation and run for the averaged initialization time files
ds_eam4 = xr.open_dataset('Control-4_eam_slice.nc')
ds_elm4 = xr.open_dataset('Control-4_elm_slice.nc')
ds_mos4 = xr.open_dataset('Control-4_mos_slice.nc')

ds_eam5 = xr.open_dataset('Control-5_eam_slice.nc')
ds_elm5 = xr.open_dataset('Control-5_elm_slice.nc')
ds_mos5 = xr.open_dataset('Control-5_mos_slice.nc')

ds_eam6 = xr.open_dataset('Control-6_eam_slice.nc')
ds_elm6 = xr.open_dataset('Control-6_elm_slice.nc')
ds_mos6 = xr.open_dataset('Control-6_mos_slice.nc')

ds_eam_Plus_1K4 = xr.open_dataset('Plus1K-4_eam_slice.nc')
ds_elm_Plus_1K4 = xr.open_dataset('Plus1K-4_elm_slice.nc')
ds_mos_Plus_1K4 = xr.open_dataset('Plus1K-4_mos_slice.nc')
ds_eam_Plus_2K4 = xr.open_dataset('Plus2K-4_eam_slice.nc')
ds_elm_Plus_2K4 = xr.open_dataset('Plus2K-4_elm_slice.nc')
ds_mos_Plus_2K4 = xr.open_dataset('Plus2K-4_mos_slice.nc')
ds_eam_Plus_3K4 = xr.open_dataset('Plus3K-4_eam_slice.nc')
ds_elm_Plus_3K4 = xr.open_dataset('Plus3K-4_elm_slice.nc')
ds_mos_Plus_3K4 = xr.open_dataset('Plus3K-4_mos_slice.nc')
ds_eam_Plus_4K4 = xr.open_dataset('Plus4K-4_eam_slice.nc')
ds_elm_Plus_4K4 = xr.open_dataset('Plus4K-4_elm_slice.nc')
ds_mos_Plus_4K4 = xr.open_dataset('Plus4K-4_mos_slice.nc')
ds_eam_Pre_Ind4 = xr.open_dataset('Pre_Ind-4_eam_slice.nc')
ds_elm_Pre_Ind4 = xr.open_dataset('Pre_Ind-4_elm_slice.nc')
ds_mos_Pre_Ind4 = xr.open_dataset('Pre_Ind-4_mos_slice.nc')

ds_eam_Plus_1K5 = xr.open_dataset('Plus1K-5_eam_slice.nc')
ds_elm_Plus_1K5 = xr.open_dataset('Plus1K-5_elm_slice.nc')
ds_mos_Plus_1K5 = xr.open_dataset('Plus1K-5_mos_slice.nc')
ds_eam_Plus_2K5 = xr.open_dataset('Plus2K-5_eam_slice.nc')
ds_elm_Plus_2K5 = xr.open_dataset('Plus2K-5_elm_slice.nc')
ds_mos_Plus_2K5 = xr.open_dataset('Plus2K-5_mos_slice.nc')
ds_eam_Plus_3K5 = xr.open_dataset('Plus3K-5_eam_slice.nc')
ds_elm_Plus_3K5 = xr.open_dataset('Plus3K-5_elm_slice.nc')
ds_mos_Plus_3K5 = xr.open_dataset('Plus3K-5_mos_slice.nc')
ds_eam_Plus_4K5 = xr.open_dataset('Plus4K-5_eam_slice.nc')
ds_elm_Plus_4K5 = xr.open_dataset('Plus4K-5_elm_slice.nc')
ds_mos_Plus_4K5 = xr.open_dataset('Plus4K-5_mos_slice.nc')
ds_eam_Pre_Ind5 = xr.open_dataset('Pre_Ind-5_eam_slice.nc')
ds_elm_Pre_Ind5 = xr.open_dataset('Pre_Ind-5_elm_slice.nc')
ds_mos_Pre_Ind5 = xr.open_dataset('Pre_Ind-5_mos_slice.nc')

ds_eam_Plus_1K6 = xr.open_dataset('Plus1K-6_eam_slice.nc')
ds_elm_Plus_1K6 = xr.open_dataset('Plus1K-6_elm_slice.nc')
ds_mos_Plus_1K6 = xr.open_dataset('Plus1K-6_mos_slice.nc')
ds_eam_Plus_2K6 = xr.open_dataset('Plus2K-6_eam_slice.nc')
ds_elm_Plus_2K6 = xr.open_dataset('Plus2K-6_elm_slice.nc')
ds_mos_Plus_2K6 = xr.open_dataset('Plus2K-6_mos_slice.nc')
ds_eam_Plus_3K6 = xr.open_dataset('Plus3K-6_eam_slice.nc')
ds_elm_Plus_3K6 = xr.open_dataset('Plus3K-6_elm_slice.nc')
ds_mos_Plus_3K6 = xr.open_dataset('Plus3K-6_mos_slice.nc')
ds_eam_Plus_4K6 = xr.open_dataset('Plus4K-6_eam_slice.nc')
ds_elm_Plus_4K6 = xr.open_dataset('Plus4K-6_elm_slice.nc')
ds_mos_Plus_4K6 = xr.open_dataset('Plus4K-6_mos_slice.nc')
ds_eam_Pre_Ind6 = xr.open_dataset('Pre_Ind-6_eam_slice.nc')
ds_elm_Pre_Ind6 = xr.open_dataset('Pre_Ind-6_elm_slice.nc')
ds_mos_Pre_Ind6 = xr.open_dataset('Pre_Ind-6_mos_slice.nc')

#creating the ensemble by concatinating the files
Ds_eam_control = xr.concat((ds_eam6, ds_eam4, ds_eam5), dim='run')
ds_eam_Control = Ds_eam_control.mean( dim=("run"), skipna=True, keep_attrs=True )
Ds_elm_control = xr.concat((ds_elm6, ds_elm4, ds_elm5), dim='run')
ds_elm_Control = Ds_elm_control.mean( dim=("run"), skipna=True, keep_attrs=True )
Ds_mos_control = xr.concat((ds_mos6, ds_mos4, ds_mos5), dim='run')
ds_mos_Control = Ds_mos_control.mean( dim=("run"), skipna=True, keep_attrs=True )

Ds_eam_Plus_1K = xr.concat((ds_eam_Plus_1K6, ds_eam_Plus_1K4, ds_eam_Plus_1K5), dim='run')
ds_eam_Plus_1K = Ds_eam_Plus_1K.mean( dim=("run"), skipna=True, keep_attrs=True )
Ds_elm_Plus_1K = xr.concat((ds_elm_Plus_1K6, ds_elm_Plus_1K4, ds_elm_Plus_1K5), dim='run')
ds_elm_Plus_1K = Ds_elm_Plus_1K.mean( dim=("run"), skipna=True, keep_attrs=True )
Ds_mos_Plus_1K = xr.concat((ds_mos_Plus_1K6, ds_mos_Plus_1K4, ds_mos_Plus_1K5), dim='run')
ds_mos_Plus_1K = Ds_mos_Plus_1K.mean( dim=("run"), skipna=True, keep_attrs=True )

Ds_eam_Plus_2K = xr.concat((ds_eam_Plus_2K6, ds_eam_Plus_2K4, ds_eam_Plus_2K5), dim='run')
ds_eam_Plus_2K = Ds_eam_Plus_2K.mean( dim=("run"), skipna=True, keep_attrs=True )
Ds_elm_Plus_2K = xr.concat((ds_elm_Plus_2K6, ds_elm_Plus_2K4, ds_elm_Plus_2K5), dim='run')
ds_elm_Plus_2K = Ds_elm_Plus_2K.mean( dim=("run"), skipna=True, keep_attrs=True )
Ds_mos_Plus_2K = xr.concat((ds_mos_Plus_2K6, ds_mos_Plus_2K4, ds_mos_Plus_2K5), dim='run')
ds_mos_Plus_2K = Ds_mos_Plus_2K.mean( dim=("run"), skipna=True, keep_attrs=True )

Ds_eam_Plus_3K = xr.concat((ds_eam_Plus_3K6, ds_eam_Plus_3K4, ds_eam_Plus_3K5), dim='run')
ds_eam_Plus_3K = Ds_eam_Plus_3K.mean( dim=("run"), skipna=True, keep_attrs=True )
Ds_elm_Plus_3K = xr.concat((ds_elm_Plus_3K6, ds_elm_Plus_3K4, ds_elm_Plus_3K5), dim='run')
ds_elm_Plus_3K = Ds_elm_Plus_3K.mean( dim=("run"), skipna=True, keep_attrs=True )
Ds_mos_Plus_3K = xr.concat((ds_mos_Plus_3K6, ds_mos_Plus_3K4, ds_mos_Plus_3K5), dim='run')
ds_mos_Plus_3K = Ds_mos_Plus_3K.mean( dim=("run"), skipna=True, keep_attrs=True )

Ds_eam_Plus_4K = xr.concat((ds_eam_Plus_4K6, ds_eam_Plus_4K4, ds_eam_Plus_4K5), dim='run')
ds_eam_Plus_4K = Ds_eam_Plus_4K.mean( dim=("run"), skipna=True, keep_attrs=True )
Ds_elm_Plus_4K = xr.concat((ds_elm_Plus_4K6, ds_elm_Plus_4K4, ds_elm_Plus_4K5), dim='run')
ds_elm_Plus_4K = Ds_elm_Plus_4K.mean( dim=("run"), skipna=True, keep_attrs=True )
Ds_mos_Plus_4K = xr.concat((ds_mos_Plus_4K6, ds_mos_Plus_4K4, ds_mos_Plus_4K5), dim='run')
ds_mos_Plus_4K = Ds_mos_Plus_4K.mean( dim=("run"), skipna=True, keep_attrs=True )

Ds_eam_Pre_Ind = xr.concat((ds_eam_Pre_Ind6, ds_eam_Pre_Ind4, ds_eam_Pre_Ind5), dim='run')
ds_eam_Pre_Ind = Ds_eam_Pre_Ind.mean( dim=("run"), skipna=True, keep_attrs=True )
Ds_elm_Pre_Ind = xr.concat((ds_elm_Pre_Ind6, ds_elm_Pre_Ind4, ds_elm_Pre_Ind5), dim='run')
ds_elm_Pre_Ind = Ds_elm_Pre_Ind.mean( dim=("run"), skipna=True, keep_attrs=True )
Ds_mos_Pre_Ind = xr.concat((ds_mos_Pre_Ind6, ds_mos_Pre_Ind4, ds_mos_Pre_Ind5), dim='run')
ds_mos_Pre_Ind = Ds_mos_Pre_Ind.mean( dim=("run"), skipna=True, keep_attrs=True )


# In[ ]:


#saving the ensemble files as a netCDF
ds_eam_Control.to_netcdf('ds_eam_Control_ens.nc')
ds_elm_Control.to_netcdf('ds_elm_Control_ens.nc')
ds_mos_Control.to_netcdf('ds_mos_Control_ens.nc')

ds_eam_Plus_1K.to_netcdf('ds_eam_Plus_1K_ens.nc')
ds_elm_Plus_1K.to_netcdf('ds_elm_Plus_1K_ens.nc')
ds_mos_Plus_1K.to_netcdf('ds_mos_Plus_1K_ens.nc')

ds_eam_Plus_2K.to_netcdf('ds_eam_Plus_2K_ens.nc')
ds_elm_Plus_2K.to_netcdf('ds_elm_Plus_2K_ens.nc')
ds_mos_Plus_2K.to_netcdf('ds_mos_Plus_2K_ens.nc')

ds_eam_Plus_3K.to_netcdf('ds_eam_Plus_3K_ens.nc')
ds_elm_Plus_3K.to_netcdf('ds_elm_Plus_3K_ens.nc')
ds_mos_Plus_3K.to_netcdf('ds_mos_Plus_3K_ens.nc')

ds_eam_Plus_4K.to_netcdf('ds_eam_Plus_4K_ens.nc')
ds_elm_Plus_4K.to_netcdf('ds_elm_Plus_4K_ens.nc')
ds_mos_Plus_4K.to_netcdf('ds_mos_Plus_4K_ens.nc')

ds_eam_Pre_Ind.to_netcdf('ds_eam_Pre_Ind_ens.nc')
ds_elm_Pre_Ind.to_netcdf('ds_elm_Pre_Ind_ens.nc')
ds_mos_Pre_Ind.to_netcdf('ds_mos_Pre_Ind_ens.nc')


# In[ ]:


#loading h1 files and creating an ensemble
ds_eam4 = xr.open_dataset('Control-4_eam_h1_slice.nc')
ds_eam5 = xr.open_dataset('Control-5_eam_h1_slice.nc')
ds_eam6 = xr.open_dataset('Control-6_eam_h1_slice.nc')

ds_eam_Plus_1K4 = xr.open_dataset('Plus1K-4_eam_h1_slice.nc')
ds_eam_Plus_2K4 = xr.open_dataset('Plus2K-4_eam_h1_slice.nc')
ds_eam_Plus_3K4 = xr.open_dataset('Plus3K-4_eam_h1_slice.nc')
ds_eam_Plus_4K4 = xr.open_dataset('Plus4K-4_eam_h1_slice.nc')
ds_eam_Pre_Ind4 = xr.open_dataset('Pre_Ind-4_eam_h1_slice.nc')

ds_eam_Plus_1K5 = xr.open_dataset('Plus1K-5_eam_h1_slice.nc')
ds_eam_Plus_2K5 = xr.open_dataset('Plus2K-5_eam_h1_slice.nc')
ds_eam_Plus_3K5 = xr.open_dataset('Plus3K-5_eam_h1_slice.nc')
ds_eam_Plus_4K5 = xr.open_dataset('Plus4K-5_eam_h1_slice.nc')
ds_eam_Pre_Ind5 = xr.open_dataset('Pre_Ind-5_eam_h1_slice.nc')

ds_eam_Plus_1K6 = xr.open_dataset('Plus1K-6_eam_h1_slice.nc')
ds_eam_Plus_2K6 = xr.open_dataset('Plus2K-6_eam_h1_slice.nc')
ds_eam_Plus_3K6 = xr.open_dataset('Plus3K-6_eam_h1_slice.nc')
ds_eam_Plus_4K6 = xr.open_dataset('Plus4K-6_eam_h1_slice.nc')
ds_eam_Pre_Ind6 = xr.open_dataset('Pre_Ind-6_eam_h1_slice.nc')

#creating the ensemble by concatinating the files
Ds_eam_control = xr.concat((ds_eam6, ds_eam4, ds_eam5), dim='run')
ds_eam_Control = Ds_eam_control.mean( dim=("run"), skipna=True, keep_attrs=True )

Ds_eam_Plus_1K = xr.concat((ds_eam_Plus_1K6, ds_eam_Plus_1K4, ds_eam_Plus_1K5), dim='run')
ds_eam_Plus_1K = Ds_eam_Plus_1K.mean( dim=("run"), skipna=True, keep_attrs=True )

Ds_eam_Plus_2K = xr.concat((ds_eam_Plus_2K6, ds_eam_Plus_2K4, ds_eam_Plus_2K5), dim='run')
ds_eam_Plus_2K = Ds_eam_Plus_2K.mean( dim=("run"), skipna=True, keep_attrs=True )

Ds_eam_Plus_3K = xr.concat((ds_eam_Plus_3K6, ds_eam_Plus_3K4, ds_eam_Plus_3K5), dim='run')
ds_eam_Plus_3K = Ds_eam_Plus_3K.mean( dim=("run"), skipna=True, keep_attrs=True )

Ds_eam_Plus_4K = xr.concat((ds_eam_Plus_4K6, ds_eam_Plus_4K4, ds_eam_Plus_4K5), dim='run')
ds_eam_Plus_4K = Ds_eam_Plus_4K.mean( dim=("run"), skipna=True, keep_attrs=True )

Ds_eam_Pre_Ind = xr.concat((ds_eam_Pre_Ind6, ds_eam_Pre_Ind4, ds_eam_Pre_Ind5), dim='run')
ds_eam_Pre_Ind = Ds_eam_Pre_Ind.mean( dim=("run"), skipna=True, keep_attrs=True )

#saving the ensemble files as an .nc
ds_eam_Control.to_netcdf('ds_eam_Control_h1_ens.nc')
ds_eam_Plus_1K.to_netcdf('ds_eam_Plus_1K_h1_ens.nc')
ds_eam_Plus_2K.to_netcdf('ds_eam_Plus_2K_h1_ens.nc')
ds_eam_Plus_3K.to_netcdf('ds_eam_Plus_3K_h1_ens.nc')
ds_eam_Plus_4K.to_netcdf('ds_eam_Plus_4K_h1_ens.nc')
ds_eam_Pre_Ind.to_netcdf('ds_eam_Pre_Ind_h1_ens.nc')


# In[ ]:


#loading the soil files and creating an ensemble
ds_eam_Controli4 = xr.open_dataset('Control-4_init_soil.nc')
ds_eam_Plus_1Ki4 = xr.open_dataset('Plus1K-4_init_soil.nc')
ds_eam_Plus_2Ki4 = xr.open_dataset('Plus2K-4_init_soil.nc')
ds_eam_Plus_3Ki4 = xr.open_dataset('Plus3K-4_init_soil.nc')
ds_eam_Plus_4Ki4 = xr.open_dataset('Plus4K-4_init_soil.nc')
ds_eam_Pre_Indi4 = xr.open_dataset('Pre_Ind-4_init_soil.nc')

ds_eam_Controli5 = xr.open_dataset('Control-5_init_soil.nc')
ds_eam_Plus_1Ki5 = xr.open_dataset('Plus1K-5_init_soil.nc')
ds_eam_Plus_2Ki5 = xr.open_dataset('Plus2K-5_init_soil.nc')
ds_eam_Plus_3Ki5 = xr.open_dataset('Plus3K-5_init_soil.nc')
ds_eam_Plus_4Ki5 = xr.open_dataset('Plus4K-5_init_soil.nc')
ds_eam_Pre_Indi5 = xr.open_dataset('Pre_Ind-5_init_soil.nc')

ds_eam_Controli6 = xr.open_dataset('Control-6_init_soil.nc')
ds_eam_Plus_1Ki6 = xr.open_dataset('Plus1K-6_init_soil.nc')
ds_eam_Plus_2Ki6 = xr.open_dataset('Plus2K-6_init_soil.nc')
ds_eam_Plus_3Ki6 = xr.open_dataset('Plus3K-6_init_soil.nc')
ds_eam_Plus_4Ki6 = xr.open_dataset('Plus4K-6_init_soil.nc')
ds_eam_Pre_Indi6 = xr.open_dataset('Pre_Ind-6_init_soil.nc')

#Combining 004, 005, 006 runs for "super ensemble" for non averaged initialization times files
DS_eam_Controli = xr.concat((ds_eam_Controli6, ds_eam_Controli4, ds_eam_Controli5), dim='init')

DS_eam_Plus_1Ki = xr.concat((ds_eam_Plus_1Ki6, ds_eam_Plus_1Ki4, ds_eam_Plus_1Ki5), dim='init')

DS_eam_Plus_2Ki = xr.concat((ds_eam_Plus_2Ki6, ds_eam_Plus_2Ki4, ds_eam_Plus_2Ki5), dim='init')

DS_eam_Plus_3Ki = xr.concat((ds_eam_Plus_3Ki6, ds_eam_Plus_3Ki4, ds_eam_Plus_3Ki5), dim='init')

DS_eam_Plus_4Ki = xr.concat((ds_eam_Plus_4Ki6, ds_eam_Plus_4Ki4, ds_eam_Plus_4Ki5), dim='init')

DS_eam_Pre_Indi = xr.concat((ds_eam_Pre_Indi6, ds_eam_Pre_Indi4, ds_eam_Pre_Indi5), dim='init')

#saving ensemble as a netCDF
DS_eam_Controli.to_netcdf('ds_soil_Controli_ens.nc')
DS_eam_Plus_1Ki.to_netcdf('ds_soil_Plus_1Ki_ens.nc')
DS_eam_Plus_2Ki.to_netcdf('ds_soil_Plus_2Ki_ens.nc')
DS_eam_Plus_3Ki.to_netcdf('ds_soil_Plus_3Ki_ens.nc')
DS_eam_Plus_4Ki.to_netcdf('ds_soil_Plus_4Ki_ens.nc')
DS_eam_Pre_Indi.to_netcdf('ds_soil_Pre_Indi_ens.nc')

