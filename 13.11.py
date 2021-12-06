# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 17:41:10 2021

@author: Chapa
"""

from osgeo import gdal, osr
import matplotlib.pyplot as plt
import rasterio
import rasterio as rs
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import xarray as xr
import rioxarray
import cartopy.crs as ccrs # cartographic coordinate reference system
import cartopy.feature as cfeature # features such as land, borders, coastlines


file = r"C:\Users\BabyCheekz\Google Drive\Data Science MSc\Final Project\Satellite Images\2TROPESS_CrIS-SNPP_L2_Standard_CH4_20210521_MUSES_R1p12_FS_F0p1.nc"
ds = xr.open_dataset(file, engine = "h5netcdf")

ds = ds.x


dssel = ds.where((-125 < ds.longitude) & (ds.longitude < -114)
                 & (49 < ds.latitude) & (ds.latitude < 55), drop = True)

plt.figure()
dssel.values()




########################################################

f = Dataset(file)
tropairs = xr.open_dataset(file, engine = "h5netcdf")

coords = tropairs.set_coords(["x"])
swap_dims = tropairs.swap_dims
stack = tropairs.stack("latitude", "longitude", "x")

print(f.variables.keys())

#methane = f.variables["x"]
methane = tropairs.x
print(methane)

for d in f.dimensions.items():
    print(d)
    print()
    
methane.dimensions
methane.shape

# lat = f.variables["latitude"]
# lon = f.variables["longitude"]
lons = tropairs.longitude
lats = tropairs.latitude

methane.plot(x=tropairs.longitude , y=tropairs.latitude,
cmap='gist_ncar', vmin=methane.min(), vmax=methane.max(),
transform=ccrs.PlateCarree())

plt.plot(lat, lon, linewidth=0.1)


########################################################################
# lat = f.variables["latitude"]
# lon = f.variables["longitude"]



da = xr.DataArray(data = methane,
                  dims = ["lats", "lons"],
                  coords = dict(
                      lon =(["x", "y"], lons),
                      lat = (["x", "y"], lats),
                      x = tropairs.x
                  ),
                  attrs=dict(
                      description = "Methane emissions",
                      units = "1"))

lat = tropairs.longitude
lon = tropairs.latitude












####################################################################
methane1 = coords.x[0:]

plt.plot(methane1, x="latitude", y="longitude",
cmap='gist_ncar', 
transform=ccrs.PlateCarree())

plt.plot(lat, lon, 'bo', linewidth=0.01)

plt.imshow(methane)

plt.plot(methane,
transform=ccrs.PlateCarree())

plt.scatter(lat, lon, s=0.1, cmap="gist_ncar")