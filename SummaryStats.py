#Module Imports
import glob
import rasterstats
from rasterstats import zonal_stats
import pandas as pd
import gdal
import numpy as np
import os
import rasterio
import geopandas
from shapely.geometry import box

#Shared Functions

def next_path(path_pattern):
    """
    Finds the next free path in an sequentially named list of files

    e.g. path_pattern = 'file-%s.txt':

    file-1.txt
    file-2.txt
    file-3.txt

    Runs in log(n) time where n is the number of existing files in sequence
    """
    i = 1

    # First do an exponential search
    while os.path.exists(path_pattern % i):
        i = i * 2

    # Result lies somewhere in the interval (i/2..i]
    # We call this interval (a..b] and narrow it down until a + 1 = b
    a, b = (i // 2, i)
    while a + 1 < b:
        c = (a + b) // 2 # interval midpoint
        a, b = (c, b) if os.path.exists(path_pattern % c) else (a, c)

    return path_pattern % b

def raster_stats(source_files, buffer_file, stats):
    """
    Loops through files and calculates specified raster stats for points which overlap with spatial extend of the file

    """
    # read buffer points
    gdf = geopandas.read_file(buffer_file)
    # prepare output
    gdf_out = gdf[['geometry']].copy()
    if isinstance(stats, str):
        stats = stats.split()
    for name in stats:
        gdf_out[name] = np.nan  # initialize statistic columns
    for fn in source_files:
        # open raster file
        with rasterio.open(fn, 'r') as src:
            affine = src.transform
            nodata = src.nodata
            bounds = box(*src.bounds)
            idxs = gdf.sindex.query(bounds)  # get points which overlap with the file spatial extent
            if len(idxs) > 0:  # only if overlap read data
                array = src.read(1)
                stats_calc = zonal_stats(
                    gdf_out.loc[idxs].geometry, array, stats=stats, affine=affine, nodata=nodata
                )
                gdf_out.loc[idxs, stats] = pd.DataFrame(index=idxs, data=stats_calc)
    return gdf_out

#GSW Euclidean Processing
def reclass(raster, out_path):
    """
    Reclassifies 1 to 0 and vice-versa in input raster

    """
    file = gdal.Open(raster)
    band = file.GetRasterBand(1)
    reRaster = band.ReadAsArray()
    where_0 = np.where(reRaster == 0)
    where_1 = np.where(reRaster == 1)
    reRaster[where_0] = 1
    reRaster[where_1] = 0
    driver = gdal.GetDriverByName('GTiff')
    path_outDs = next_path(out_path)
    file2 = driver.Create(path_outDs, file.RasterXSize, file.RasterYSize, 1, gdal.GPI_RGB)
    file2.GetRasterBand(1).WriteArray(reRaster)

    # spatial ref system
    proj = file.GetProjection()
    georef = file.GetGeoTransform()
    meta = file.GetMetadata()
    colors = file.GetRasterBand(1).GetRasterColorTable()

    file2.SetProjection(proj)
    file2.SetGeoTransform(georef)
    file2.SetMetadata(meta)
    file2.GetRasterBand(1).SetRasterColorTable(colors)

# Creating a list of the gsw file names for euclidean distance calculation
gsw_files = glob.glob("E:/dirk test run/GSW2/*.tif")

euc_out_path ="E:/dirk test run/GSW_Reclass2/GSWReclass%g.tif"
# applying reclassify to the gsw files
for gsw_file in gsw_files:
    reclass(gsw_file, euc_out_path)

# Creating a list of the reclassed files
gsw_reclass_files = glob.glob("E:/dirk test run/GSW_Reclass2/*.tif")
gsw_euc_folder = "E:/dirk test run/GSW_Euc2"  # NOTE: see note below
print(gsw_reclass_files)

# Calculating euclidean distance for the reclassed files
for in_path in gsw_reclass_files:
    src_ds = gdal.Open(in_path)
    srcband = src_ds.GetRasterBand(1)

    basename = os.path.basename(in_path)
    dst_filename = os.path.join(gsw_euc_folder, basename.replace('.tif', '_euc.tif'))
    print(dst_filename)
    drv = gdal.GetDriverByName('GTiff')
    dst_ds = drv.Create(dst_filename,
                        src_ds.RasterXSize, src_ds.RasterYSize, 1,
                        gdal.GetDataTypeByName('Int16'))

    dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
    dst_ds.SetProjection(src_ds.GetProjectionRef())

    dstband = dst_ds.GetRasterBand(1)

    gdal.ComputeProximity(srcband, dstband, ["VALUES='1'", "DISTUNITS=PIXEL"])
    srcband = None
    dstband = None
    src_ds = None
    dst_ds = None

#Raster Stats

#DEM
#Defence
dem_files = glob.glob("E:/dirk test run/DEM/*.tif")
defences_buffer = "E:/dirk test run/Point buffer/ptbuff2.shp"
dem_raster_stats = "mean std"

dem_stats_def = raster_stats(dem_files, defences_buffer, dem_raster_stats)
dem_stats_def = dem_stats_def.rename(columns={'mean': 'DEM_Mean', 'std': 'DEM_std'})
print(dem_stats_def)

# #Not Defence
# not_defences_buffer = "E:/dirk test run/Point buffer/ptbuff2.shp"
#
# dem_stats_not = raster_stats(dem_files, not_defences_buffer, dem_raster_stats)
# dem_stats_not = dem_stats_not.rename(columns={'mean': 'DEM_Mean', 'std': 'DEM_std'})
# print(dem_stats_not)

#HAND
#Defence
hand_files = glob.glob("E:/dirk test run/HAND/*.tif")
hand_raster_stats = "mean range std"

hand_stats_def = raster_stats(hand_files, defences_buffer, hand_raster_stats)
hand_stats_def = hand_stats_def.rename(columns={'mean': 'HAND_Mean', 'range': 'HAND_range', 'std': 'HAND_std'}).drop(['geometry'], axis=1)
print(hand_stats_def)

# # Not Defence
# hand_stats_not = raster_stats(hand_files, not_defences_buffer, hand_raster_stats)
# hand_stats_not = hand_stats_not.rename(columns={'mean': 'HAND_Mean', 'range': 'HAND_range', 'std': 'HAND_std'}).drop(['geometry'], axis=1)
# print(hand_stats_not)

#GSW
#Defence
gsw_files = glob.glob("E:/dirk test run/GSW/*.tif")
gsw_raster_stats = "mean"

gsw_stats_def = raster_stats(gsw_files, defences_buffer, gsw_raster_stats)
gsw_stats_def = gsw_stats_def.rename(columns={'mean': 'GSW_Mean'}).drop(['geometry'], axis=1)
print(gsw_stats_def)

# # Not Defence
# gsw_stats_not = raster_stats(gsw_files, not_defences_buffer, gsw_raster_stats)
# gsw_stats_not = gsw_stats_not.rename(columns={'mean': 'GSW_Mean'}).drop(['geometry'], axis=1)
# print(gsw_stats_not)

#GSW Euclidean
#Defence
gsw_euc_files = glob.glob("E:/dirk test run/GSW_Euc2/*.tif")
gsw_euc_raster_stats = "mean"

gsw_euc_stats_def = raster_stats(gsw_euc_files, defences_buffer, gsw_euc_raster_stats)
gsw_euc_stats_def = gsw_euc_stats_def.rename(columns={'mean': 'GSW_Euc_Mean'}).drop(['geometry'], axis=1)
print(gsw_euc_stats_def)

# # Not Defence
# gsw_euc_stats_not = raster_stats(gsw_euc_files, not_defences_buffer, gsw_euc_raster_stats)
# gsw_euc_stats_not = gsw_euc_stats_not.rename(columns={'mean': 'GSW_Euc_Mean'}).drop(['geometry'], axis=1)
# print(gsw_euc_stats_not)

#GFP Euclidean
#Defence
gfp_euc_files = glob.glob("E:/dirk test run/GFP_Euc/*.tif")
gfp_euc_raster_stats = "mean"

gfp_euc_stats_def = raster_stats(gfp_euc_files, defences_buffer, gfp_euc_raster_stats)
gfp_euc_stats_def = gfp_euc_stats_def.rename(columns={'mean': 'GFP_Euc_Mean'}).drop(['geometry'], axis=1)
print(gfp_euc_stats_def)

# # Not Defence
# gfp_euc_stats_not = raster_stats(gfp_euc_files, not_defences_buffer, gfp_euc_raster_stats)
# gfp_euc_stats_not = gfp_euc_stats_not.rename(columns={'mean': 'GFP_Euc_Mean'}).drop(['geometry'], axis=1)
# print(gfp_euc_stats_not)

#Raster Stats Combined
combined_stats_def = pd.concat([dem_stats_def, hand_stats_def, gsw_stats_def, gsw_euc_stats_def, gfp_euc_stats_def], axis=1)
print(combined_stats_def)

combined_stats_def.to_csv(r'E:/dirk test run/Test1.csv')