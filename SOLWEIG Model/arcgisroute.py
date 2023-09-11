from skimage.graph import route_through_array
import rasterio
from pyproj import Transformer
import networkx as nx

def raster2array(rasterfn):
    raster = rasterio.open(rasterfn)
    band = raster.read(1)
    #array = band.ReadAsArray()
    return band

# TAKES IN LAT LONG AND RETURNS PIXEL COORDS
def coord2pixelOffset(rasterfn, lat, long):

    # convert from wgs84 to 26912 to pixel
    raster = rasterio.open(rasterfn)
    # https://epsg.io/transform#s_srs=4326&t_srs=26912&ops=1728&x=-111.9296140&y=33.4207069
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:26912", always_xy=True)
    lat, long = transformer.transform(lat, long)
    geotransform = raster.transform
    lat, long = geotransform * (lat, long)

    return lat, long

def createPath(CostSurfacefn, costSurfaceArray, startCoord, stopCoord):
    # coordinates to array index
    startCoordX = startCoord[0]
    startCoordY = startCoord[1]
    startIndexX, startIndexY = coord2pixelOffset(CostSurfacefn, startCoordX, startCoordY)
    print(startIndexX, startIndexY)
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
    raster = rasterio.open(CostSurfacefn)
    gt = raster.transform

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

import geopandas as gpd
from rasterio.features import geometry_mask
import networkx as nx
import rasterio
import matplotlib.pyplot as plt
from pyproj import Proj, transform
import time

def convertWGS84ToPixelSpace(long, lat, raster):
    # convert from wgs84 to 26912 to pixel
    raster_crs = raster.crs.to_proj4()
    lat_lon_crs = Proj(init='epsg:4326')  # Assuming WGS84 coordinates (EPSG:4326)

    # Transform latitude and longitude to the CRS of the raster
    x, y = transform(lat_lon_crs, raster_crs, long, lat)

    # Get the pixel coordinates from the transformed x and y
    pixel_x, pixel_y = raster.index(x, y)

    return pixel_x, pixel_y

def downsample_raster(raster, factor):
    """
    Down-sample a raster by a given factor.

    Args:
        raster (numpy.ndarray): The original raster data.
        factor (int): The down-sampling factor.

    Returns:
        numpy.ndarray: The down-sampled raster data.
    """
    downsampled_data = raster[::factor, ::factor]
    return downsampled_data

def makeRoutingGraph(mean_radiant_temperature, obstacle_mask, dateTimeString, downsample_factor=10):
    downsampled_mrt = downsample_raster(mean_radiant_temperature, downsample_factor)
    downsampled_obstacle_mask = downsample_raster(obstacle_mask, downsample_factor)

    G = nx.grid_2d_graph(downsampled_mrt.shape[0], downsampled_mrt.shape[1])
    print('making graph')
    # Add weights to the graph based on the cost matrix
    for node in G.nodes:
        i, j = node
        G.nodes[node]['cost'] = downsampled_mrt[i, j] + downsampled_obstacle_mask[i, j] * float('inf')
    nx.write_gpickle(G, 'output/{0}_graph.gpickle'.format(dateTimeString))
    return G

import os
import numpy as np
def getRouteAStar(startCoord, stopCoord, dateTimeString):
    print('starting route')
    mrt_file_path = 'output/{0}_mrt.tif'.format(dateTimeString)  # expected format is "2023-03-30_12:00"
    startCoord = startCoord  # start point in wgs84 (long, lat)
    stopCoord = stopCoord  # psych north in wgs84 (long, lat)
    down_sample_rate = 3

    print('reading files')
    # Load mean radiant temperature raster
    with rasterio.open(mrt_file_path) as src:
        mean_radiant_temperature = src.read(1)
        print("Raster CRS:", src.crs)
        startCoordX = startCoord[0]
        startCoordY = startCoord[1]
        startIndexX, startIndexY = convertWGS84ToPixelSpace(startCoordX, startCoordY, src)
        startIndexX = int(startIndexX / down_sample_rate)
        startIndexY = int(startIndexY / down_sample_rate)
        stopCoordX = stopCoord[0]
        stopCoordY = stopCoord[1]
        stopIndexX, stopIndexY = convertWGS84ToPixelSpace(stopCoordX, stopCoordY, src)
        stopIndexX = int(stopIndexX / down_sample_rate)
        stopIndexY = int(stopIndexY / down_sample_rate)

    downsampled_mrt = downsample_raster(mean_radiant_temperature, down_sample_rate)

    # Load obstacle mask shape file as a raster

    # Load the obstacle shapefile
    obstacle_shapefile = 'Maps/Tempe_MaskedBuildingsRoads.shp'
    obstacle_gdf = gpd.read_file(obstacle_shapefile)

    # Create a mask from the shapefile that matches the mean radiant temperature raster
    with rasterio.open(mrt_file_path) as src:
        obstacle_mask = geometry_mask(obstacle_gdf.geometry, out_shape=src.shape, transform=src.transform, invert=True)

    # Now, obstacle_mask contains the obstacle information where True represents obstacles

    # Use obstacle_mask in your pathfinding algorithm
    print('files read')
    print(mean_radiant_temperature.shape)

    print(startIndexX, startIndexY)
    print(stopIndexX, stopIndexY)
    # if the graph file exists, load it, otherwise make it
    graph_path = 'output/{0}_graph.gpickle'.format(dateTimeString)
    if os.path.exists(graph_path):
        G = nx.read_gpickle(graph_path)
    else:
        # Create a graph representing the raster cells and their connections
        G = makeRoutingGraph(mean_radiant_temperature, obstacle_mask, dateTimeString, down_sample_rate)
    # convert cells to pixel coordinates
    print('graph made')
    print('finding path')
    path_time = time.time()
    # Use a pathfinding algorithm to find the most comfortable route
    path = nx.astar_path(G, source=(startIndexX, startIndexY), target=(stopIndexX, stopIndexY), weight='cost')
    print(f'path found in {time.time() - path_time}')
    # Create a copy of the mean radiant temperature raster
    output_path = 'path_result.tif'
    print('writing file')
    with rasterio.open(mrt_file_path) as src:
        profile = src.profile
        # You can modify the profile as needed, e.g., set the data type to integer for path visualization:
        profile['dtype'] = 'int32'
        profile['count'] = 2  # Increase the band count by 1
        with rasterio.open(output_path, 'w', **profile) as dst:
            # Copy the mean radiant temperature data to the output raster
            dst.write(mean_radiant_temperature, 1)

            # Create a path mask with a unique value for the path cells (e.g., 9999)
            path_mask = downsampled_mrt.copy()
            path_mask.fill(0)
            for node in path:
                i, j = node
                path_mask[i, j] = 100

            dst.write(path_mask, 2)

    print('file written')
    # Visualize the result (mean radiant temperature with the path)
    with rasterio.open(output_path) as src:
        fig, ax = plt.subplots(figsize=(8, 8))
        mrt = src.read(1)
        # Determine the minimum and maximum values in the mean radiant temperature data
        vmin = max(np.min(mrt),0)
        vmax = np.max(mrt)
        print(vmin, vmax)
        plt.imshow(mrt, cmap='coolwarm', vmin=vmin, vmax=vmax)

        # Overlay the path on top of the mean radiant temperature
        path_mask = src.read(2)
        marker_coords = np.argwhere(path_mask != 0)
        if len(marker_coords) > 0:
            x_coords = marker_coords[:, 1]  # x-coordinates of markers
            y_coords = marker_coords[:, 0]  # y-coordinates of markers
            plt.scatter(x_coords, y_coords, c='purple', s=5, marker='o', label='Path Markers')

        plt.title('Mean Radiant Temperature with Path')
        plt.colorbar(label='Mean Radiant Temperature', ticks=[vmin, vmax])
        plt.show()


startCoord = (-111.93952587328305, 33.423795079832) # brickyard in wgs84 (long, lat)
stopCoord = (-111.92961401637315, 33.42070688780706) #psych north in wgs84 (long, lat)
getRouteAStar(startCoord, stopCoord, '2023-4-8-2100')