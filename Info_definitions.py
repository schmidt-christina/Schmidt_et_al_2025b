from collections import OrderedDict
import xarray as xr
import numpy as np
import cosima_cookbook as cc
import cf_xarray
from metpy.interpolate import cross_section
import pyproj
import matplotlib.path as mpath

path_database = '/g/data/e14/cs6673/Ross_salinity/Python_scripts_published/'
path_output = '/g/data/e14/cs6673/Ross_salinity/data_analysis/'
path_plots = '/g/data/e14/cs6673/Ross_salinity/plots/'

exptdict = OrderedDict([
    ('ctrl',  # control
     {'expt': '01deg_jra55v13_ryf9091_rerun_for_easterlies',
      'expt_str': 'CONTROL',
      'col': 'k',
      'u_file': '/g/data/ik11/inputs/JRA-55/RYF/v1-3/RYF.u_10.1990_1991.nc',
      'v_file': '/g/data/ik11/inputs/JRA-55/RYF/v1-3/RYF.v_10.1990_1991.nc'}),
    ('wind_50_down_zonal', # zonal winds between Ross and Amundsen Sea decreased by 50%
     {'expt': 'access-om2-01_ryf_wind_50_down_zonal',
      'expt_str': 'WIND$-$',
      'pert_lon': [-167, -115],
      'mask_wind': [
          [360-167, 360-167, 360-140, 360-115, 360-115, 360-140, 360-167],
          [-80.5, -72, -69 , -69, -76, -77, -80.5]],
      'col': 'r',
      'u_file': ('/home/142/cs6673/work/Ross_salinity/forcing_perturbed/' +
                 'RYF_winds_decreased_50.u_10.1990_1991.nc'),
      'v_file': '/g/data/ik11/inputs/JRA-55/RYF/v1-3/RYF.v_10.1990_1991.nc'}),
    ('mw_50_down', # meltwater input in Amundsen Sea decreased by 50%
     {'expt': 'access-om2-01_ryf_meltwater_50_down',
      'expt_str': 'MW$-$',
      'pert_lon': [-150, -88],
      'col': 'b'}),
    ('mw_50_down_over_1_yr', # meltwater input in Amundsen Sea linearly 
     # decreased over 1 year by 50%
     {'expt': 'access-om2-01_ryf_meltwater_50_down_over_1_yr',
      'expt_str': '',
      'pert_lon': [-150, -88],
      'col': 'limegreen'}),
    ('era5',
     {'expt_str': 'ERA5',
      'u_file': '/g/data/rt52/era5/single-levels/monthly-averaged/10u/*/10u_era5_moda_sfc*.nc',
      'v_file': '/g/data/rt52/era5/single-levels/monthly-averaged/10v/*/10v_era5_moda_sfc*.nc'})
])

def yearly_mean(var):
    month_length = var.time.dt.days_in_month
    weights_month = (month_length.groupby('time.year') /
                     month_length.groupby('time.year').sum())
    var = (var * weights_month).groupby('time.year').sum()
    var = var.rename({'year': 'time'})
    return var

def shelf_mask_isobath(var, output_mask=False):
    '''
    Masks ACCESS-OM2-01 variables by the region polewards of the 1000m isobath as computed using 
    a script contributed by Adele Morrison.
    Only to be used with ACCESS-OM2-0.1 output!
    '''
    contour_file = np.load('/g/data/ik11/grids/Antarctic_slope_contour_1000m.npz')
    
    shelf_mask = contour_file['contour_masked_above']
    yt_ocean = contour_file['yt_ocean']
    xt_ocean = contour_file['xt_ocean']
    
    # in this file the points along the isobath are given a positive value, the points outside (northwards) 
    # of the isobath are given a value of -100 and all the points on the continental shelf have a value of 0 
    # so we mask for the 0 values 
    shelf_mask[np.where(shelf_mask!=0)] = np.nan
    shelf_mask = shelf_mask+1
    shelf_map = np.nan_to_num(shelf_mask)
    shelf_mask = xr.DataArray(shelf_mask, coords = [('yt_ocean', yt_ocean), ('xt_ocean', xt_ocean)])
    shelf_map = xr.DataArray(shelf_map, coords = [('yt_ocean', yt_ocean), ('xt_ocean', xt_ocean)])
    
    # then we want to multiply the variable with the mask so we need to account for the shape of the mask. 
    # The mask uses a northern cutoff of 59S.
    masked_var = var.sel(yt_ocean = slice(-90, -59.03)) * shelf_mask
    
    if output_mask == True:
        return masked_var, shelf_map
    else:
        return masked_var

def select_bottom_values(var, lat_north=-59):
    var = var.sel(yt_ocean=slice(-90, lat_north))
    var = var.where(var != 0)
    
    ht = cc.querying.getvar(
        exptdict['ctrl']['expt'], 'ht', session=cc.database.create_session(), n=1) 
    ht = ht.sel(yt_ocean=slice(-90, lat_north))
    land_mask = (ht*0).fillna(1)

    # select bottom values
    depth_array = var*0 + var.st_ocean
    max_depth = depth_array.max(dim='st_ocean', skipna=True)

    var_bottom = var.where(depth_array.st_ocean >= max_depth)
    var_bottom = var_bottom.sum(dim='st_ocean')
    var_bottom = var_bottom.where(land_mask == 0)
    return var_bottom