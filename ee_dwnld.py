import ee
import pandas as pd 
import numpy as np
import geopandas as gpd
import folium
import csv
import argparse
import time

def readCommands():
  '''
  Read commandline arguments
  '''
  p = argparse.ArgumentParser(description=("Creating ISEI tiles of AOI"))
  p.add_argument("--AOI-geopackage", dest ="AOI_gpkg", type=str, default='../../other_data/country_extents/malawi_extent.gpkg', help=("Input AOI geopackage pathway"))
  p.add_argument("--AOI-name", dest ="AOI_name", type=str, default='MALAWI', help=("AOI name to be used in file output"))
  p.add_argument("--EE-collection", dest ="EE_collection", type=str, default='LANDSAT/LC08/C01/T1_TOA', help=("Earth Engine image collection desired"))
  p.add_argument("--date-start", dest ="date_start", type=str, default='2015-06-01', help=("Start date to obtain imagery"))
  p.add_argument("--date-end", dest ="date_end", type=str, default='2017-06-30', help=("End date to obtain imagery"))
  p.add_argument("--RS-shapefile", dest ="rs_shp", type=str, default='../../other_data/WRS2_descending_0/wrs2_descending.shp', help=("Input reference system shapefile pathway (eg WRS for Landsat)"))
  p.add_argument("--out-resolution", dest ="resolution", type=int, default=30, help=("Output resolution (m)"))
  p.add_argument("--cloudMask", dest="cloudMask", type=bool, default=True, help=("Apply cloud masking"))
  p.add_argument("--origRGB", dest="origRGB", type=bool, default=True, help=("Output original RGB image"))
  p.add_argument("--outMap", dest="outMap", type=bool, default=True, help=("Output webpage of scene, AOI and overlap map"))
  p.add_argument("--outSURF", dest="outSURF", type=str, default='ISEI', help=("Output SURF image"))
  cmdargs = p.parse_args()
  return cmdargs

class SURF_Creation(object):
    """
        SURF data handler
    """
    def __init__(self,AOI_gpkg,EE_collection,date_start,date_end,rs_shp,resolution,cloudMask,origRGB,outSURF,AOI_name,outMap):
        """
            Connecting to EE server
        """
        try:
            ee.Initialize()
            print('Earth Engine initialized successfully')
            try: # Begin processing
                self.ImageryIntersectAOI(EE_collection,AOI_gpkg,rs_shp,date_start,date_end,resolution,cloudMask,origRGB,outSURF,AOI_name,outMap) # Find area for processing
            except:
                print('Image processing failed')
        except:
            print('Earth Engine initialisation failed')

    def ImageryIntersectAOI(self,EE_collection,AOI_gpkg,rs_shp,date_start,date_end,resolution,cloudMask,origRGB,outSURF,AOI_name,outMap):
        """
            Finding my AOI/Landsat intersection to process 
        """
        # Define image collection to search with date
        self.collection = ee.ImageCollection(EE_collection).filterDate(date_start,date_end)

        # Read in bounds of shapefiles: AOI bounds and WRS scenes to find overlap
        self.bounds = gpd.read_file(AOI_gpkg)
        self.wrs = gpd.GeoDataFrame.from_file(rs_shp)
        self.wrs_intersection = self.wrs[self.wrs.intersects(self.bounds.geometry[0])]

        # Create these return values as variables
        self.paths, self.rows = self.wrs_intersection['PATH'].values, self.wrs_intersection['ROW'].values
        print(len(self.paths),'Intersecting Scenes Found')

        self.SceneImageryProcess(resolution,cloudMask,origRGB,outSURF,AOI_name,outMap) # Begin processing

    def SceneImageryProcess(self,resolution,cloudMask,origRGB,outSURF,AOI_name,outMap):
        """
            Processing per scene of intersection
        """
        x = 0 # Start count
        for i, scene in self.wrs_intersection.iterrows():

                # Select scene for processing
                print('-- Processing Scene',x+1,'of path:',scene.PATH, 'row:',scene.ROW,' --')
                        
                # Extract geometry JSON and ee.Geometry of WRS scene
                row_gpd = gpd.GeoDataFrame()
                row_gpd['geometry'] = None
                row_gpd.loc[0, 'geometry'] = scene.geometry
                scene_geo = row_gpd.__geo_interface__  
                self.scene_extent = ee.Geometry(scene_geo['features'][0]['geometry']) 
                        
                # Extract geometry JSON and ee.Geometry of intersection between 
                # the WRS scene and input shapefile as the AOI to process
                poly_intersection = gpd.overlay(self.bounds,row_gpd)
                overlay_geom = poly_intersection.__geo_interface__
                self.AOI = ee.Geometry(overlay_geom['features'][0]['geometry'])
                        
                if(outMap):
                    ############################################################
                    # Create webpage of scenes and AOIs being processed  #
                    ############################################################
                    # Setting up Folium map parameters
                    xy = np.asarray(self.bounds.centroid[0].xy).squeeze()
                    center = list(xy[::-1])
                    zoom = 7
                    self.m = folium.Map(location=center, zoom_start=zoom, control_scale=True)

                    # Add the AOI geometry object
                    self.m.add_child(folium.GeoJson(self.bounds.__geo_interface__, name='Area of Study', 
                                    style_function=lambda x: {'color': 'red', 'alpha': 0}))

                    # Add the scene geometry object
                    name = 'path: %03d, row: %03d' % (scene.PATH, scene.ROW)                        
                    self.g = folium.GeoJson(scene.geometry.__geo_interface__, name=name)
                    self.g.add_child(folium.Popup(name))
                    
                    # Add the AOI-scene intersection geometry object
                    self.m.add_child(folium.GeoJson(poly_intersection.__geo_interface__, name='Area of Overlap', 
                                                    style_function=lambda x: {'color': 'green', 'alpha': 0}))  
                    self.g.add_to(self.m)

                    # Add all layers created to the folium map
                    folium.LayerControl().add_to(self.m)
                            
                    # Export webpage of scene and intersection AOI completed 
                    try:
                        self.m.save('./wrs_AOI.html') # overwrite per scene
                        print('Success: webpage of WRS Scene',str(x+1),'and AOI created.')
                    except:
                        print('Failed: unable to create WRS Scene and AOI website.')

                # Create image collection at that scene
                path_use = self.collection.filter(ee.Filter.eq('WRS_PATH',int(scene.PATH)))    
                self.pathrow = path_use.filter(ee.Filter.eq('WRS_ROW',int(scene.ROW)))

                ############################################################
                # Begin image selection and SURF processing  #
                ############################################################

                ### STEP ONE: Find least cloudy image for that scene 
                print('-- Sorting image collection by cloud cover --')
                self.sort_by_cloud = self.pathrow.sort('CLOUD_COVER_LAND')
                self.least_cloudy = self.sort_by_cloud.first()
                self.id = self.least_cloudy.getInfo()['id']
                self.id = self.id[24:]
                cloudiness = self.least_cloudy.getInfo()['properties']['CLOUD_COVER_LAND']
                print('Selected image',self.id,'with cloud cover of',cloudiness,'for processing')

                if(origRGB):
                    # Create RGB of surface reflectance image
                    sr_image = ee.Image('LANDSAT/LC08/C01/T1_SR/'+str(self.id))
                    print('-- Preparing SR RGB image for export --')
                    orig_image_RGB = sr_image.visualize(bands=['B4','B3','B2'],min=0,max=3000)
                    output_orig_RGB = orig_image_RGB.clip(self.AOI)
                    output1 = ee.batch.Export.image.toDrive(image=output_orig_RGB,description=str(AOI_name)+'_SR_RGB_R'+str(scene.ROW)+'_P'+str(scene.PATH)+'_SC'+str(x+1),folder='SR_RGB_data',scale=resolution)
                    output1.start()
                    print("SR RGB Export: process sent to your drive")
                
               
                # Create SURF image
                if(cloudMask):
                    scored = ee.Algorithms.Landsat.simpleCloudScore(self.least_cloudy) # cloud score band
                    mask = scored.select(['cloud']).lte(15) # 15% threshold

                    # Apply the cloud mask to the image 
                    print('-- Applying cloud mask --')
                    self.masked = self.least_cloudy.updateMask(mask)
                else:
                    self.masked = self.least_cloudy 

                ### STEP TWO: Extract bands needed for SURF calculations
                self.b2 = self.masked.select('B2') # BLUE
                self.b3 = self.masked.select('B3') # GREEN
                self.b4 = self.masked.select('B4') # RED
                self.nir = self.masked.select('B5') # NIR
                self.swir1 = self.masked.select('B6') # SWIR1
                self.swir2 = self.masked.select('B7') #SWIR2

                if outSURF == 'ISEI':
                    ### STEP THREE: Calculate Impervious Surface Extent Index (ISEI)
                    print('-- Calculating ISEI --')
                    # Concatenate all bands so far into one image
                    self.img = ee.Image.cat([self.b4, self.b3, self.b2, self.nir, self.swir1, self.swir2])

                    # ISEI = (pRed*pBlue) / [(pNIR/pSWIR1)+(pNIR/pSWIR2)] from Twumasi, et al (2020)
                    ISEI = self.img.expression('(RED*BLUE)/((NIR/SWIR1)+(NIR/SWIR2))',{
                        'RED': self.img.select('B4'),
                        'BLUE': self.img.select('B2'),
                        'NIR': self.img.select('B5'),
                        'SWIR1': self.img.select('B6'),
                        'SWIR2': self.img.select('B7')
                    })

                    ### STEP FOUR: Export data to Google Drive 
                    print('-- Preparing ISEI data for export --')
                    self.ISEI_final = ISEI.select(['B4'],['ISEI']) # rename for later use
                    ISEI_float = self.ISEI_final.float()
                    outputSURF = ISEI_float.clip(self.AOI)
                
                elif outSURF == 'ENDISI':
                    ### STEP THREE: Calculate Enhanced Normalised Difference Impervious Surface Index (ENDISI)
                    print('-- Calculating ENDISI --')
                    # Concatenate all bands so far into one image
                    self.img = ee.Image.cat([self.b4, self.b3, self.b2, self.nir, self.swir1, self.swir2])

                    # MNDWI = (pGreen - pSWIR1) / (pGreen + pSWIR1)
                    mndwi = self.img.expression('(GREEN - SWIR1)/(GREEN + SWIR1)',{
                                        'GREEN': self.img.select('B3'), 
                                        'SWIR1': self.img.select('B6')
                    })
                    mndwi_renamed = mndwi.select(['B3'],['MNDWI'])  # rename the band for later use
                                    
                    ### Calculating Alpha ###
                    print('-- Calculating Alpha --')
                    # pBlue.mean
                    mean_B2 = self.img.reduceRegion(reducer=ee.Reducer.mean(), geometry=self.scene_extent, scale=30, maxPixels=1e9).get('B2')   
                    
                    # pSWIR1/pSWIR2
                    swir_ratio = self.swir1_toa.divide(self.swir2_toa)
                    swir_ratio_renamed = swir_ratio.select(['B6'],['SWIR_RATIO']) # renaming resultant band
                    
                    # (pSWIR1/pSWIR2).mean
                    mean_swir_ratio = swir_ratio_renamed.reduceRegion(reducer=ee.Reducer.mean(), geometry=self.scene_extent, scale=resolution, maxPixels=1e9).get('SWIR_RATIO') 
                    
                    # MNDWI^2
                    mndwi_sq = mndwi_renamed.multiply(mndwi_renamed)
                    
                    # (MNDWI^2).mean
                    mean_mndwi_sq = mndwi_sq.reduceRegion(reducer=ee.Reducer.mean(), geometry=self.scene_extent, scale=resolution, maxPixels=1e9).get('MNDWI') 
                    
                    # a = (2 * pBlue.mean) / ((pSWIR1/pSWIR2).mean + (MNDWI^2).mean)
                    alpha = (2*mean_B2.getInfo())/(mean_swir_ratio.getInfo()+mean_mndwi_sq.getInfo())
                    
                    # ENDISI = [pBlue - a*(pSWIR1/pSWIR2 + MNDWI^2)] / [pBlue + a*(pSWIR1/pSWIR2 + MNDWI^2)] from Chen, et al (2019)
                    ENDISI = self.img.expression('(B2 - ('+str(alpha)+'*(SWIR_RATIO + MNDWI_SQ)))/(B2 + ('+str(alpha)+'*(SWIR_RATIO + MNDWI_SQ)))',{
                                'B2': self.img.select('B2'),
                                'SWIR_RATIO': swir_ratio_renamed.select('SWIR_RATIO'),
                                'MNDWI_SQ': mndwi_sq.select('MNDWI')
                    })
                    
                    ### STEP FOUR: Export data to Google Drive 
                    print('-- Preparing ENDISI data for export --')
                    self.ENDISI_final = ENDISI.select(['B4'],['ENDISI']) # rename for later use
                    ENDISI_float = self.ENDISI_final.float()
                    outputSURF = ENDISI_float.clip(self.AOI)
                else:
                    print('There\'s a problem with your IS index choice')
                    break
                
                output2 = ee.batch.Export.image.toDrive(image=outputSURF,description=str(AOI_name)+'_'+str(outSURF)+'_R'+str(scene.ROW)+'_P'+str(scene.PATH)+'_SC'+str(x+1),folder='SURF_data',scale=resolution)
                output2.start()
                print("SURF Export: process sent to your drive")
            
                x += 1

        print("All scenes completed!")

if __name__=="__main__":
    start_time = time.time()
    cmd = readCommands()
    
    make_SURF = SURF_Creation(AOI_gpkg=cmd.AOI_gpkg,
    EE_collection=cmd.EE_collection,
    date_start=cmd.date_start,
    date_end=cmd.date_end,
    rs_shp=cmd.rs_shp,
    resolution=cmd.resolution,
    origRGB=cmd.origRGB,
    cloudMask=cmd.cloudMask,
    outSURF=cmd.outSURF,
    outMap=cmd.outMap,
    AOI_name=cmd.AOI_name)

    print("--- %s seconds ---" % (time.time() - start_time))