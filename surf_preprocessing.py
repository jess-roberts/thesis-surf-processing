import glob
import os
import rasterio
import gdal
import numpy as np
import osr 
import matplotlib.pyplot as plt
from rasterio.merge import merge
from rasterio.fill import fillnodata
from sklearn.preprocessing import Normalizer, StandardScaler

def findRasters(tifDir,AOI_name):
    # criteria to find the geotiffs
    search_criteria = str(AOI_name)+'*'
    tifs = os.path.join(tifDir, search_criteria)

    # taking all the files that match and globbing together
    all_tifs = glob.glob(tifs)

    # opening each of these file names and putting them into a list for the merge
    tifs_4_mosaic = []
    for fp in all_tifs:
        src = rasterio.open(fp)
        tifs_4_mosaic.append(src)
    
    print(len(tifs_4_mosaic),'GeoTiffs found')
        
    return all_tifs

class tiffHandle(object):
    def __init__(self,in_filename,out_filename):
        self.input = in_filename
        self.output = out_filename
        self.classed = False
        self.norm = False
        self.readReprojRaster()

    def readReprojRaster(self):
        '''
        Read a geotiff in to RAM
        '''
        print("-- Reading in raster --")
        # open a dataset object
        self.ds = gdal.Warp(str(self.output), self.input, xRes=30, yRes=-30, dstSRS='EPSG:3857', outputType=gdal.GDT_Float32, srcNodata=np.nan, dstNodata=np.nan)
    
        # read data from geotiff object
        self.nX=self.ds.RasterXSize             # number of pixels in x direction
        self.nY=self.ds.RasterYSize             # number of pixels in y direction
        # geolocation tiepoint
        transform_ds = self.ds.GetGeoTransform()# extract geolocation information
        self.xOrigin=transform_ds[0]       # coordinate of x corner
        self.yOrigin=transform_ds[3]       # coordinate of y corner
        self.pixelWidth=transform_ds[1]    # resolution in x direction
        self.pixelHeight=transform_ds[5]   # resolution in y direction
        # read data. Returns as a 2D numpy array
        self.data=self.ds.GetRasterBand(1).ReadAsArray(0,0,self.nX,self.nY)

    def classRaster(self):
        """
            Reclassing 99th percentile data
        """
        print('-- Reclassing raster --')
        self.p99 = np.nanpercentile(self.data, 99)
        self.data[(self.data >= self.p99)] = 255 # impervious
        self.data[(self.data <= self.p99)] = 0 # pervious

        self.classed = True 

    def writeTiff(self,epsg=3857):
        """
            Create output raster
        """
        # set geolocation information (note geotiffs count down from top edge in Y)
        geotransform = (self.xOrigin, self.pixelWidth, 0, self.yOrigin, 0, self.pixelHeight)

        # load data in to geotiff object
        dst_ds = gdal.GetDriverByName('GTiff').Create(self.output, self.nX, self.nY, 1, gdal.GDT_Float32)
        dst_ds.SetGeoTransform(geotransform)    # specify coords
        srs = osr.SpatialReference()            # establish encoding
        srs.ImportFromEPSG(epsg)                # set crs
        dst_ds.SetProjection(srs.ExportToWkt()) # export coords to file
        dst_ds.GetRasterBand(1).WriteArray(self.data)  # write image to the raster
        dst_ds.GetRasterBand(1).SetNoDataValue(np.nan)  # set no data value
        dst_ds.FlushCache()                     # write to disk
        dst_ds = None
        
        print("Image written to",self.output)

######### Processing ISEI layer ###########
filepw = '/Volumes/Jaffa 1TB/Uni/Thesis/ISEI/'
country = 'MALAWI'

# Finding and standardising each country scene file
tif_list = findRasters(str(filepw)+str(country)+'/',str(country))
tifs_4_mosaic = []
for tif in tif_list:
    name = str(tif[:-4])+'_CLASSED.tif'
    isei_tif_obj = tiffHandle(in_filename=tif,out_filename=name)
    isei_tif_obj.classRaster()
    isei_tif_obj.writeTiff()
    src = rasterio.open(name)
    tifs_4_mosaic.append(src)

# Merging all the files found
print('-- Merging all processed files --')
surf_mosaic, out_trans = merge(tifs_4_mosaic) # merge returns single array

# Updating the metadata
out_meta = src.meta.copy()
out_meta.update({"driver": "Gtiff",
                        "height": surf_mosaic.shape[1],
                        "width": surf_mosaic.shape[2],
                        "transform": out_trans})

# Writing out the merged file
print('-- Writing merged file --')
surf_output = str(filepw)+str(country)+'/'+str(country)+'_ISEI_MERGED.tif'

with rasterio.open(surf_output, "w", **out_meta) as dest:
        dest.write(surf_mosaic)

print('Impervious surface image written to',surf_output)

