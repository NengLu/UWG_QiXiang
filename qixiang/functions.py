#!/usr/bin/env python
# coding: utf-8
#
# Neng Lu
# nengl@student.unimelb.edu.au
# ANU & Unimelb
# Canberra, Australia
# 
# Version: 1.0
# First version 14 May, 2020
# Last modified 22 May, 2020


import numpy as np
import math
from osgeo import gdal
from osgeo import osr

def testimport():
    print("It works!")

#-----------------------------------------------------------#
# data=(ny,nx) version
def array2geotiff_yx(fname, data, latRange, lonRange, dtype):   
    """
    save GeoTiff file from the array of dem data
    input:
    fname: save file name
    data: elevation data, an array in size of (n_lat,n_lon) 
    latRange: range of latitude, an array as [minlat,maxlat]
    lonRange: range of longitude, an array as [minlon,maxlon]
    dtype: dtype in gdal, as gdal.GDT_Byte or gdal.GDT_Float32
    """   
    nx = data.shape[1]
    ny = data.shape[0]
    xmin,xmax,ymin,ymax = [lonRange[0],lonRange[1],latRange[0],latRange[1]]
    dx = (xmax - xmin) / float(nx)
    dy = (ymax - ymin) / float(ny)
    geotransform = (xmin, dx, 0, ymax, 0, -dy)
    dst = gdal.GetDriverByName('GTiff').Create(fname, nx, ny, 1, dtype)
    dst.SetGeoTransform(geotransform) 
    dst.GetRasterBand(1).WriteArray(data)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    dst.SetProjection(srs.ExportToWkt())  
    dst.FlushCache() 

def get_extent(extent_s,res_deg):
    """
    Get the real coordinate extent of the area in the grid of earth2014
    -----------
    Input:
    extent_s: the rough coordinate extent, a tuple like (-180,180,-90,90)
    res_deg: the grid interval of earth2014 data
    -----------
    Output:
    extent_r: the read coordinate extent 
    """
    lats = np.arange((-90+res_deg/2),(90-res_deg/4),res_deg)
    lons = np.arange((-180+res_deg/2),(180-res_deg/4),res_deg)
    nlat = len(lats) 
    nlon = len(lons)
    minlon1,maxlon1,minlat1,maxlat1 = (lons.min(),lons.max(),lats.min(),lats.max())
    
    minlon2,maxlon2,minlat2,maxlat2 = extent_s
    
    minX_index = np.around((minlon2-minlon1)/res_deg).astype(int)
    maxX_index = np.around((maxlon2-minlon1)/res_deg).astype(int)
    if (maxX_index-minX_index)%2 == 0:
        maxX_index = maxX_index - 1
    minY_index = np.around((minlat2-minlat1)/res_deg).astype(int) 
    maxY_index = np.around((maxlat2-minlat1)/res_deg).astype(int)
    if (maxY_index-minY_index)%2 == 0:
        maxY_index = maxY_index - 1
    
    extent_t = (lons[minX_index],lons[maxX_index],lats[minY_index],lats[maxY_index])
    return extent_t
    
def get_data(data_s,extent_s,res_deg):
    """
    Get the data of the target area
    -----------
    Input:
    data_s: the data of the source area, an array of [ny,nx]
    extent_s: the rough coordinate extent, a tuple like (-180,180,-90,90)
    res_deg: the grid interval of earth2014 data, a value, unit: degree
    -----------
    Output:
    extent_r: the read coordinate extent of the target area
    data_r: the data of the target area
    """
    lats = np.arange((-90+res_deg/2),(90-res_deg/4),res_deg)
    lons = np.arange((-180+res_deg/2),(180-res_deg/4),res_deg)
    nlat = len(lats) 
    nlon = len(lons)
    minlon1,maxlon1,minlat1,maxlat1 = (lons.min(),lons.max(),lats.min(),lats.max())
    
    minlon2,maxlon2,minlat2,maxlat2 = extent_s
    
    minX_index = np.around((minlon2-minlon1)/res_deg).astype(int)
    maxX_index = np.around((maxlon2-minlon1)/res_deg).astype(int)
    if (maxX_index-minX_index)%2 == 0:
        maxX_index = maxX_index - 1
    minY_index = np.around((minlat2-minlat1)/res_deg).astype(int) 
    maxY_index = np.around((maxlat2-minlat1)/res_deg).astype(int)
    if (maxY_index-minY_index)%2 == 0:
        maxY_index = maxY_index - 1
    
    extent_t = (lons[minX_index],lons[maxX_index],lats[minY_index],lats[maxY_index])
    
    data_m = np.flipud(data_s.copy())
    data_t = data_m[minY_index:(maxY_index+1),minX_index:(maxX_index+1)]
    data_t = np.flipud(data_t)

    return extent_t, data_t
    

#-----------------------------------------------------------#
def cal_dis_LngLat(lon1,lat1,lon2,lat2):
    latitude1 = (math.pi/180)*lat1
    latitude2 = (math.pi/180)*lat2
    longitude1 = (math.pi/180)*lon1
    longitude2= (math.pi/180)*lon2
    #{arccos[sinb*siny+cosb*cosy*cos(a-x)]}*R
    R = 6378.137
    d = math.acos(math.sin(latitude1)*math.sin(latitude2)+ math.cos(latitude1)*math.cos(latitude2)*math.cos(longitude2-longitude1))*R
    return d

def cal_azi_LngLat(lon1,lat1,lon2,lat2):
    lat1_rad = lat1 * math.pi / 180
    lon1_rad = lon1 * math.pi / 180
    lat2_rad = lat2 * math.pi / 180
    lon2_rad = lon2 * math.pi / 180

    y = math.sin(lon2_rad - lon1_rad) * math.cos(lat2_rad)
    x = math.cos(lat1_rad) * math.sin(lat2_rad) - \
        math.sin(lat1_rad) * math.cos(lat2_rad) * math.cos(lon2_rad - lon1_rad)

    azi = math.atan2(y, x) * 180 / math.pi
    azi = float((azi + 360.0) % 360.0)
    return azi

def cal_azi(x1,y1,x2,y2):
    y = y2-y1
    x = x2-x1
    azi = math.atan2(y, x) * 180 / math.pi
    azi = float((-azi + 90.0) % 360.0)
    return azi

def cal_dis(x1,y1,x2,y2):
    y = y2-y1
    x = x2-x1
    d = math.sqrt(y**2+x**2) 
    return d
    
def cal_azi_river_LngLat(river_xy):
    river_x = river_xy[:,0]
    river_y = river_xy[:,1]
    N = len(river_x)
    azi = np.zeros(N)
    for i in range(0,N-1):
        azi[i] = cal_azi_LngLat(river_x[i],river_y[i],river_x[i+1],river_y[i+1])
    return azi

def cal_azi_river(river_xy):
    river_x = river_xy[:,0]
    river_y = river_xy[:,1]
    N = len(river_x)
    azi = np.zeros(N)
    for i in range(0,N-1):
        azi[i] = cal_azi(river_x[i],river_y[i],river_x[i+1],river_y[i+1])
    return azi
    
    
    
    
    
    
    
