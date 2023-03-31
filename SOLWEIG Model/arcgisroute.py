from skimage.graph import route_through_array
import numpy as np
from osgeo import gdal
from pyproj import Transformer

def raster2array(rasterfn):
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1)
    array = band.ReadAsArray()
    return array

# TAKES IN LAT LONG AND RETURNS PIXEL COORDS
def coord2pixelOffset(rasterfn, lat, long):
    # convert from wgs84 to 26912 to pixel
    raster = gdal.Open(rasterfn)
    #https://epsg.io/transform#s_srs=4326&t_srs=26912&ops=1728&x=-111.9296140&y=33.4207069
    transformer = Transformer.from_crs("EPSG:4326","EPSG:26912", always_xy=True)
    lat, long = transformer.transform(lat, long)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    xOffset = int((lat - originX) / pixelWidth)
    yOffset = int((long - originY) / pixelHeight)
    return xOffset, yOffset

def createPath(CostSurfacefn, costSurfaceArray, startCoord, stopCoord):
    # coordinates to array index
    startCoordX = startCoord[0]
    startCoordY = startCoord[1]
    startIndexX, startIndexY = coord2pixelOffset(CostSurfacefn, startCoordX, startCoordY)

    stopCoordX = stopCoord[0]
    stopCoordY = stopCoord[1]
    stopIndexX, stopIndexY = coord2pixelOffset(CostSurfacefn, stopCoordX, stopCoordY)

    # create path
    indices, weight = route_through_array(costSurfaceArray, (startIndexY, startIndexX), (stopIndexY, stopIndexX),
                                          geometric=True, fully_connected=True)
    print(indices) # coordinates in pixel space
    print(weight) # total cost of route
    return indices, weight

def getRouteCoords(CostSurfacefn, startCoord, stopCoord, outputPathfn=None):
    costSurfaceArray = raster2array(CostSurfacefn)  # creates array from cost surface raster
    pixel_path, total_weight = createPath(CostSurfacefn, costSurfaceArray, startCoord, stopCoord)  # creates path array

    # convert pixel_path to WGS84
    raster = gdal.Open(CostSurfacefn)
    gt = raster.GetGeoTransform()

    x_min = gt[0]
    x_size = gt[1]
    y_min = gt[3]
    y_size = gt[5]
    coords = [(mx * x_size + x_min, my * y_size + y_min) for (mx,my) in pixel_path] # return cords in NAD83
    return coords

# costsurfacefn is the solweig mrt tif file generated
# start/stop is a cooridnate in wgs84 (long, lat)
# outputPath is a file name for the output tif file
def getRoute(startCoord, stopCoord, dateTimeString):
    CostSurfacefn = 'output/{0}_mrt.tif'.format(dateTimeString) # expected format is "2023-03-30_12:00"
    startCoord = startCoord  # start point in wgs84 (long, lat)
    stopCoord = stopCoord  # psych north in wgs84 (long, lat)
    # outputPathfn = 'Maps/PathFromBrickyardToNobleTestNAD.tif'
    return getRouteCoords(CostSurfacefn, startCoord, stopCoord)