# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 16:33:14 2021

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


############################################ PLOTTING SENTINEL-5P DATA - nc file ############################################

file1 = r"C:\Users\BabyCheekz\Google Drive\Data Science MSc\Final Project\Satellite Images\s5p combine\s5p 1.nc"
file2 = r"C:\Users\BabyCheekz\Google Drive\Data Science MSc\Final Project\Satellite Images\s5p combine\s5p2.nc"
file3 = r"C:\Users\BabyCheekz\Google Drive\Data Science MSc\Final Project\Satellite Images\s5p combine\s5p 3.nc"

###------------------------------------------------------------------------

#looking at the dataset and into the PRODUCT group in the hierarchical groups as this is where the data on methane
#is stored
s5p_1 = xr.open_dataset(file1, engine = "netcdf4", group = "PRODUCT")
s5p_2 = xr.open_dataset(file2, engine = "netcdf4", group = "PRODUCT")
s5p_3 = xr.open_dataset(file3, engine = "netcdf4", group = "PRODUCT")
#s5p_3 = s5p_3.drop('time')

#selecting the feature in the PRODUCT group which we want to plot
#now trying to plot 3 data products to get a full view of methane emissions over north america for a given day
methane1 = s5p_1.methane_mixing_ratio_bias_corrected[0,0:-1,:]
methane2 = s5p_2.methane_mixing_ratio_bias_corrected[0,0:-1,:]
methane3 = s5p_3.methane_mixing_ratio_bias_corrected[0,0:-1,:]


plt.figure(dpi=500) #open a new figure window and set the resolution
ax1 = plt.axes(projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND) #fill in the land areas
ax1.coastlines(linewidth=0.5) #use the default low-resolution coastline


gl = ax1.gridlines(draw_labels=True, linewidth=0.3, color='black') # default is to label all axes.
gl.top_labels=False #turn off two of them.
gl.right_labels=False
gl.bottom_labels=False
gl.left_labels=False

methane1.plot(x='longitude', y='latitude',
cmap='gist_ncar', vmin=methane1.min(), vmax=methane1.max(),
transform=ccrs.PlateCarree(), add_colorbar=False)

#source: https://polar.ncep.noaa.gov/ngmmf_python/Python_tutorial_grib.pdf
#add_colorbar=False to stop the colorbar from the 3 stacked plots from appearing


############################################ PLOTTING TROPESS-AIRS DATA ############################################
#trying with tropess-airs satellite image
file4 = r"C:\Users\BabyCheekz\Google Drive\Data Science MSc\Final Project\Satellite Images\2TROPESS_CrIS-SNPP_L2_Standard_CH4_20210521_MUSES_R1p12_FS_F0p1.nc"

tropairs = xr.open_dataset(file4, engine = "h5netcdf")

methane4 = tropairs.x[0:-1]

plt.figure(dpi=500) #open a new figure window and set the resolution
ax2 = plt.axes(projection=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND) #fill in the land areas
ax2.coastlines(linewidth=0.5) #use the default low-resolution coastline

gl = ax2.gridlines(draw_labels=True, linewidth=0.3, color='black') # default is to label all axes.
gl.top_labels=False #turn off two of them.
gl.right_labels=False
gl.bottom_labels=False
gl.left_labels=False


lons = tropairs.longitude
lats = tropairs.latitude

#lon,lat= np.meshgrid(lons,lats)

methane4.plot(x="target", y="level",
cmap='gist_ncar', vmin=methane4.min(), vmax=methane4.max(),
transform=ccrs.PlateCarree())
#image is weird, I don't think I'm plotting the right things

#plt.scatter(lat, lon, s=0.1)


######################################################################################################
file5 = r"C:\Users\BabyCheekz\Downloads\TROPESS_CrIS-JPSS1_L2_Standard_CH4_20211104_MUSES_R1p12_FS_F0p3.nc"
tropairs2 = xr.open_dataset(file5, engine = "h5netcdf")

methane5 = tropairs2.x[0:-1]

plt.figure(dpi=500) #open a new figure window and set the resolution
ax2 = plt.axes(projection=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND) #fill in the land areas
ax2.coastlines(linewidth=0.5) #use the default low-resolution coastline

gl = ax2.gridlines(draw_labels=True, linewidth=0.3, color='black') # default is to label all axes.
gl.top_labels=False #turn off two of them.
gl.right_labels=False
gl.bottom_labels=False
gl.left_labels=False


lons2 = tropairs2.longitude
lats2 = tropairs2.latitude

latlon = np.dstack((lons2, lats2))

#lon,lat= np.meshgrid(lons,lats)

methane5.plot(x=latlon, y="target",
cmap='gist_ncar', vmin=methane5.min(), vmax=methane5.max(),
transform=ccrs.PlateCarree())




################################################################################################

methane1.plot(x='longitude', y='latitude',
cmap='gist_ncar', vmin=methane1.min(), vmax=methane1.max(),
transform=ccrs.PlateCarree(), add_colorbar=False)

methane2.plot(x='longitude', y='latitude',
cmap='gist_ncar', vmin=methane2.min(), vmax=methane2.max(),
transform=ccrs.PlateCarree(), add_colorbar=False)

methane3.plot(x='longitude', y='latitude',
cmap='gist_ncar', vmin=methane3.min(), vmax=methane3.max(),
transform=ccrs.PlateCarree())
plt.show()




################################################################################################
#error saying to make coords into numpy array first
trop_xy = np.array(tropairs.coords)
s5p_xy = np.array(s5p_1.coords)




trop_lats = []

for i in lat, lon:
    coord = (lat[i], lon[i])
    trop_lats.append(coord)
    
lat = np.array(lat)
lon = np.array(lon)



tropairs.transpose('latitude', 'longitude', 'x', ...)


###############################################################
#trying something I saw online
da = xr.open_rasterio(file4)


methane4.plot(tropairs.x)
    
lats = file4.variables["latitude"][:]
lons = file4.variables["longitude"][:]
methane5 = file4.variables["x"][:]



fig, axs = plt.subplots(figsize=(15, 10), nrows=2,ncols=1,gridspec_kw={'height_ratios': [20,1.5]},constrained_layout=True)
pcm=axs[0].pcolormesh(lons,lats,methane5,cmap='viridis')
cbar=fig.colorbar(pcm,cax=axs[1], extend='both', orientation='horizontal')
cbar.set_label('PM 2.5 [$\mu$g m$^{-3}]$')


DS_new = xr.merge([methane3, methane4], compat="override")


DS_new.plot.scatter(x='longitude', y='latitude',
cmap='gist_ncar', vmin=DS_new.min(), vmax=DS_new.max(),
transform=ccrs.PlateCarree())



