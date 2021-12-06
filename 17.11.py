# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 14:55:22 2021

@author: Chapa
"""

from osgeo import gdal, osr
import matplotlib.pyplot as plt
import rasterio
import rasterio as rs
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import xarray as xr
import rioxarray
import cartopy.crs as ccrs # cartographic coordinate reference system
import cartopy.feature as cfeature # features such as land, borders, coastlines



file = r"C:\Users\BabyCheekz\Downloads\TES-Aura_L2-CH4-Nadir_2015-12_v007_Lite-v02.00(1).nc"
tes_ch4 =  xr.open_dataset(file, engine = "netcdf4")
tes_lonlat = xr.open_dataset(file, engine = "netcdf4")

methane = tes_ch4.Species
lat = tes_lonlat.Latitude
lon = tes_lonlat.Longitude



plt.figure(dpi=500) #open a new figure window and set the resolution
ax1 = plt.axes(projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND) #fill in the land areas
ax1.coastlines(linewidth=0.5) #use the default low-resolution coastline


gl = ax1.gridlines(draw_labels=True, linewidth=0.3, color='black') # default is to label all axes.
gl.top_labels=False #turn off two of them.
gl.right_labels=False
gl.bottom_labels=False
gl.left_labels=False


methane.plot(x=lon.any(), y=lat.any(), cmap = "gist_ncar", vmin=methane.min(), vmax=methane.max(),
transform=ccrs.PlateCarree())



######################################################################
#creating new data array
file1 = r"C:\Users\BabyCheekz\Google Drive\Data Science MSc\Final Project\Satellite Images\2TROPESS_CrIS-SNPP_L2_Standard_CH4_20210521_MUSES_R1p12_FS_F0p1.nc"

tropairs = xr.open_dataset(file1, engine = "h5netcdf")

#x group holding values for methane emissions
methane = tropairs.x.target

lat = methane.latitude
lon = methane.longitude
 
da = xr.DataArray(data=methane, dims=["lat", "lon"], coords=[lat,lon])


