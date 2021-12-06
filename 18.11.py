# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 13:38:33 2021

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



######################################################################
file1 = r"C:\Users\BabyCheekz\Google Drive\Data Science MSc\Final Project\Satellite Images\2TROPESS_CrIS-SNPP_L2_Standard_CH4_20210521_MUSES_R1p12_FS_F0p1.nc"

tropairs = xr.open_dataset(file1, engine = "h5netcdf")

methane = tropairs.x
level = tropairs.level[8]
target = tropairs.target


plt.contourf(methane, level, target)
plt.show()


tropairs.plot(x=target, y=level,
cmap='gist_ncar', vmin=tropairs.min(), vmax=tropairs.max(),
transform=ccrs.PlateCarree(), add_colorbar=False)



lat = tropairs.x.longitude
lon = tropairs.x.latitude

lat.shape #(24057)
lon.shape #(24057)

#putting them together to make 2D array of coordinate points
result = np.vstack((lat, lon)).T
print(result)
