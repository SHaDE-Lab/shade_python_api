import math
import os
from random import random

import rasterio
from pyproj import Transformer
import time
import osmnx as ox
from shapely.geometry import LineString, Point
import simplekml
import geojson
import networkx as nx
import geopandas as gp
import pandas as pd

default_mrt_file_path = 'output/default_mrt.tif'
default_graph_path = 'output/default_graph_networked.graphml'
def convert_to_pixel(lat, long, raster):
    # Get the pixel coordinates from the lat long
    return raster.index(lat, long)

def make_walking_network_graph(mean_radiant_temp, date_time_string):
    # bounding box of tempe campus 
    # TODO CHANGE TO REALLY BE THE BOUNDING BOX
    north = 33.42796506379531
    west = -111.94015084151006

    south = 33.41592636331673
    east = -111.92634308211977
    granularity = 1

    # ox.settings.use_cache = True
    #G = ox.graph_from_bbox(north, south, east, west, network_type='walk')
    G = convert_shp2graph3('Maps/tempeShape/tempeCampusSelected.shp', make_G_bidi=True, name='tempe')
    G = ox.project_graph(G)
    # ox.plot_graph(G)
    mrt_data = mean_radiant_temp.read(1)

    for u, v, data in G.edges(data=True):
        # num samples is the edge distance
        if 'geometry' in data:
            # Get the edge's geometry as a LINESTRING
            edge_geometry = data['geometry']
            # create a LineString object
            edge_line = LineString(edge_geometry)
            num_samples = int(edge_line.length * granularity)
            num_samples = num_samples if num_samples > 0 else 1
            # Sample points along the edge's LineString
            total_mrt = 0
            for i in range(num_samples + 1):
                alpha = i / num_samples
                # Interpolate the point along the LineString
                point = edge_line.interpolate(alpha, normalized=True)
                coords = (point.x, point.y)
                u = convert_to_pixel(coords[0], coords[1], mean_radiant_temp)
                if u[0] < 0 or u[0] >= mean_radiant_temp.height or u[1] < 0 or u[1] >= mean_radiant_temp.width:
                    continue
                total_mrt += max(mrt_data[u], 0)
        else:
            # If no 'geometry' key, interpolate between u and v
            total_mrt = 0
            u_coords = (G.nodes[u]['x'], G.nodes[u]['y'])
            v_coords = (G.nodes[v]['x'], G.nodes[v]['y'])
            edge_length = ox.distance.euclidean(u_coords[1], u_coords[0], v_coords[1], v_coords[0])
            num_samples = int(edge_length * granularity)
            num_samples = num_samples if num_samples > 0 else 1
            for i in range(num_samples + 1):
                alpha = i / num_samples
                coords = (u_coords[0] + alpha * (v_coords[0] - u_coords[0]), v_coords[1] + alpha * (v_coords[1] - u_coords[1]))
                u = convert_to_pixel(coords[0], coords[1], mean_radiant_temp)
                # make sure in bounds of raster
                if u[0] < 0 or u[0] >= mean_radiant_temp.height or u[1] < 0 or u[1] >= mean_radiant_temp.width:
                    continue
                total_mrt += max(mrt_data[u], 0)

        mean_mrt_value = total_mrt / (num_samples + 1)

        # gross hack to make sure theres no 0 weight edges
        # check if its NAN
        if mean_mrt_value == 0 or total_mrt == 0 or math.isnan(mean_mrt_value) or math.isnan(total_mrt):
            mean_mrt_value = 1000000000
            total_mrt = 10000000000
        # Add the custom MRT attribute and cost attribute to the edge
        data['mrt'] = mean_mrt_value
        data['cost'] = total_mrt
    ox.save_graphml(G, 'output/{0}_graph_{1}.graphml'.format(date_time_string, 'networked'))
    return G


def get_route(start_coord, stop_coord, date_time_string):
    # if the graph file exists, load it, otherwise make it
    graph_path = 'output/{0}_graph_{1}.graphml'.format(date_time_string, 'networked')
    mrt_file_path = 'output/{0}_mrt.tif'.format(date_time_string)  # expected format is "2023-03-30-1200"
    attribute_types = {'cost': float, 'mrt': float, 'oneway': str}
    node_types = {'osmid': float}

    if not os.path.exists(mrt_file_path):
        mrt_file_path = default_mrt_file_path
        # G = ox.load_graphml(default_graph_path, edge_dtypes=attribute_types)
        G = make_walking_network_graph(rasterio.open(mrt_file_path), date_time_string)
    elif os.path.exists(graph_path):
        G = ox.load_graphml(graph_path, edge_dtypes=attribute_types, node_dtypes=node_types)
        ox.plot_graph(G)
    else:
        # Create a graph representing the raster cells and their connections
        with rasterio.open(mrt_file_path) as src:
            G = make_walking_network_graph(src, date_time_string)
    path_time = time.time()
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:26912", always_xy=True)
    # Find the nearest network nodes to the origin and destination
    long, lat = transformer.transform(start_coord[0], start_coord[1])
    orig_node = ox.distance.nearest_nodes(G, long, lat)
    long, lat = transformer.transform(stop_coord[0], stop_coord[1])
    dest_node = ox.distance.nearest_nodes(G, long, lat)
    # Calculate the route using Dijkstra's algorithm (shortest path)
    route = ox.routing.shortest_path(G, orig_node, dest_node, weight='cost')
    # route is list of node IDs constituting the shortest path
    print('path found in {0}'.format(time.time() - path_time))
    # Convert the route to a GeoDataFrame
    print('route: {0}'.format(route))
    # Plot the graph with the route highlighted
    ox.plot_graph_route(G, route)
    # convert the route to a kml file

    kml = convert_to_kml(G, route)
    # calculate stats from the route
    statistics = calculate_statistics(G, route)
    print(statistics)

    geojson = convert_to_geoJSON(G, route)

    return kml, statistics, geojson

def convert_to_geoJSON(G, route):
    # for each node in the route, get the lat/long
    route_coords = []
    for node in route:
        node_coords = ((G.nodes[node]['lon'], G.nodes[node]['lat']))
        route_coords.append(node_coords)
    # convert the route to a geojson
    route_geojson = geojson.LineString(route_coords)
    return route_geojson

def convert_to_kml(G, route):
    kml = simplekml.Kml()
    for i in range(len(route) - 1):
        source_node = route[i]
        target_node = route[i + 1]
        # Use .edges() to get edge data for the current pair of nodes
        edge_data = G.get_edge_data(source_node, target_node)[0]
        if edge_data is not None:
            # Access the edge attributes
            # Create a placemark for the edge with name and coordinates
            name = edge_data.get('name', 'unnamed')
            # Get the coordinates of the source and target nodes
            source_coords = (G.nodes[source_node]['x'], G.nodes[source_node]['y'])
            target_coords = (G.nodes[target_node]['x'], G.nodes[target_node]['y'])

            # Create a LineString in the KML
            coords = [source_coords, target_coords]
            kml.newlinestring(name=name, coords=coords)
    kml.save('output/route.kml')
    # Save the KML to a string
    kml_string = kml.kml()

    return kml_string

def calculate_statistics(G, route):
    # calculate the length of the route
    length = 0.0
    # calculate the total mrt of the route
    mrt = 0.0
    for i in range(len(route) - 1):
        source_node = route[i]
        target_node = route[i + 1]
        # Use .edges() to get edge data for the current pair of nodes
        edge = G.get_edge_data(source_node, target_node)[0]
        length += float(edge['length'])
        mrt += float(edge['mrt'])
    average_mrt = mrt / (len(route) - 1)
    # return stats dictionary
    return {'length': length, 'mrt': mrt, 'average_mrt': average_mrt}


def convert_shp2graph(p, make_G_bidi=True, name='unamed'):
    """
    Converts shapefile to routable networkx graph.

    Parameters
    ----------
    p : str, File path - allowed formats geojson and ESRI Shapefile and other formats Fiona can read and write
    make_G_bidi : bool, if True, assumes linestrings are bidirectional
    name : str, Optional name of graph

    Returns
    -------
    G : graph
    """

    # Load shapefile into GeoDataFrame
    gdf = gp.read_file(p)

    # shapefile needs to include minimal: geometry linestring and the length computed (e.g. in QGIS)
    if 'LineLength' not in gdf.columns:
        raise Exception('Shapefile is invalid: length not in attributes:\n{}'.format(gdf.columns))

    if not gdf.geometry.map(lambda x: isinstance(x, LineString)).all():
        s_invalid_geo = gdf.geometry[~gdf.geometry.map(lambda x: isinstance(x, LineString))]
        raise Exception('Shapefile is invalid: geometry not all linestring \n{}'.format(s_invalid_geo))

    # Initialize empty list to store interpolated points
    interpolated_points = []

    # Interpolate points along LineStrings
    for line in gdf.geometry:
        # Extract points from LineString
        points = list(line.coords)

        # Interpolate points between consecutive pairs
        for i in range(len(points) - 1):
            # Interpolate 10 points between each consecutive pair of points
            interpolated_points.extend(LineString([points[i], points[i + 1]]).interpolate(x) for x in range(1, 11))

    # Create a GeoDataFrame from the interpolated points
    gdf_interpolated = gp.GeoDataFrame(geometry=interpolated_points, crs=gdf.crs)

    # Combine original GeoDataFrame with interpolated points GeoDataFrame
    gdf_combined = gp.GeoDataFrame(geometry=gp.GeoSeries(gdf.geometry.tolist() + gdf_interpolated.geometry.tolist()),
                                    crs=gdf.crs)

    # Convert to MultiDiGraph
    G = nx.MultiDiGraph(name=name, crs=gdf_combined.crs)

    # Add nodes to graph
    for idx, point in enumerate(gdf_combined.geometry):
        x, y = point.xy
        G.add_node(idx, x=x[0], y=y[0])

    # Add edges to graph
    for idx, line in gdf.iterrows():
        for i in range(len(line.geometry.coords) - 1):
            start_point = Point(line.geometry.coords[i])
            end_point = Point(line.geometry.coords[i + 1])

            # Find start and end nodes in the combined GeoDataFrame
            start_node = gdf_combined[gdf_combined.intersects(start_point)].index[0]
            end_node = gdf_combined[gdf_combined.intersects(end_point)].index[0]

            G.add_edge(start_node, end_node, **line.drop('geometry').to_dict())

    if make_G_bidi:
        # Adding bidirectional edges
        for idx, line in gdf.iterrows():
            for i in range(len(line.geometry.coords) - 1):
                start_point = Point(line.geometry.coords[i])
                end_point = Point(line.geometry.coords[i + 1])

                # Find start and end nodes in the combined GeoDataFrame
                start_node = gdf_combined[gdf_combined.intersects(start_point)].index[0]
                end_node = gdf_combined[gdf_combined.intersects(end_point)].index[0]
                G.add_edge(end_node, start_node, **line.drop('geometry').to_dict())

    return G

def convert_shp2graph3(p, make_G_bidi = True, name='unamed'):
    # Load shapefile into GeoDataFrame
    gdf = gp.read_file(p)

    # shapefile needs to include minimal: geometry linestring and the length computed (e.g. in QGIS)
    if 'LineLength' not in gdf.columns:
        raise Exception('Shapefile is invalid: length not in attributes:\n{}'.format(gdf.columns))

    if not gdf.geometry.map(lambda x: isinstance(x, LineString)).all():
        s_invalid_geo = gdf.geometry[~gdf.geometry.map(lambda x: isinstance(x, LineString))]
        raise Exception('Shapefile is invalid: geometry not all linestring \n{}'.format(s_invalid_geo))

    # Create a directed graph from the edges GeoDataFrame
    G_dir = ox.graph_from_gdfs(gdf_edges=gdf)

    # Extract nodes and edges GeoDataFrames from the directed graph
    gdf_nodes, gdf_edges = ox.graph_to_gdfs(G_dir)

    # Create graph
    G = ox.utils_graph.graph_from_gdfs(gdf_nodes, gdf_edges)

    if make_G_bidi:
        # Make graph bidirectional
        G = G.to_undirected()
    return G

def convert_shp2graph2(p, make_G_bidi = True, name='unamed'):
    """
    Converts shapefile to routable networkx graph.
    
    Parameters
    ----------
    p : str, File path - allowed formats geojson and ESRI Shapefile and other formats Fiona can read and write
    make_G_bidi : bool, if True, assumes linestrings are bidirectional
    name : str, Optional name of graph
    
    Returns
    -------
    G : graph
    """

    # Load shapefile into GeoDataFrame
    gdf = gp.read_file(p)
    
    # shapefile needs to include minimal: geometry linestring and the length computed (e.g. in QGIS)
    if 'LineLength' not in gdf.columns:
        raise Exception('Shapefile is invalid: length not in attributes:\n{}'.format(gdf.columns))

    if  not gdf.geometry.map(lambda x: type(x) ==  LineString).all():
        s_invalid_geo = gdf.geometry[gdf.geometry.map(lambda x: type(x) ==  LineString)]
        raise Exception('Shapefile is invalid: geometry not all linestring \n{}'.format(s_invalid_geo))
   
    # Compute the start- and end-position based on linestring 
    gdf['Start_pos'] = gdf.geometry.apply(lambda x: x.coords[0])
    gdf['End_pos'] = gdf.geometry.apply(lambda x: x.coords[-1])

    # Create Series of unique nodes and their associated position
    s_points = gdf.Start_pos._append(gdf.End_pos).reset_index(drop=True)
    s_points = s_points.drop_duplicates()   

    # Add index of start and end node of linestring to geopandas DataFrame
    df_points = pd.DataFrame(s_points, columns=['Start_pos'])
    df_points['FNODE_'] = df_points.index
    gdf = pd.merge(gdf, df_points, on='Start_pos', how='inner')

    df_points = pd.DataFrame(s_points, columns=['End_pos'])
    df_points['TNODE_'] = df_points.index
    gdf = pd.merge(gdf, df_points, on='End_pos', how='inner')
    
    # Bring nodes and their position in form needed for osmnx (give arbitrary osmid (index) despite not osm file)
    df_points.columns = ['pos', 'osmid']
    # assign random base 10 int as osmid
    df_points['osmid'] = df_points['osmid'].apply(lambda x: int(float(random() * 10000)))

    df_points[['x', 'y']] = df_points['pos'].apply(pd.Series)
    df_node_xy = df_points.drop('pos', axis=1)

    # Create Graph Object
    G = nx.MultiDiGraph(name=name, crs=gdf.crs)
    
    # Add nodes to graph
    for node, data in df_node_xy.T.to_dict().items():
        G.add_node(node, **data)
        
    # Add edges to graph
    for i, row  in gdf.iterrows():
        dict_row  = row.to_dict()
        if 'geometry' in dict_row: del dict_row['geometry']
        G.add_edge(u_for_edge=dict_row['FNODE_'], v_for_edge=dict_row['TNODE_'], **dict_row)
        
    if make_G_bidi:
        gdf.rename(columns={'Start_pos': 'End_pos',
                   'End_pos': 'Start_pos',
                   'FNODE_': 'TNODE_', 
                   'TNODE_': 'FNODE_', }, inplace=True)
        
        # Add edges to graph
        for i, row  in gdf.iterrows():
            dict_row  = row.to_dict()
            if 'geometry' in dict_row: del dict_row['geometry']
            G.add_edge(u_for_edge=dict_row['FNODE_'], v_for_edge=dict_row['TNODE_'], **dict_row)
        
#         G = G.to_undirected() # Some function in osmnx do not work anymore
        
    # Log information
#     log('Graph has been successfully generated /n {}'.format(nx.info(G)))
#     log('Show graph data structure EDGE'.format(G.get_edge_data(*list(G.edges())[0])))
#     log('Show graph data structure NO
    return G

# TODO - add function to determine if an edge is walk only
def isWalkOnly(edge):
    return False

if __name__ == '__main__':
    brickyard = (-111.93952587328305, 33.423795079832)  # brickyard in wgs84 (long, lat)
    psych_north = (-111.92961401637315, 33.42070688780706)  # psych north in wgs84 (long, lat)
    date_time_string = 'default'
    # read tempeGraph.xml
    G = ox.graph_from_xml('Maps/tempeGraph.xml')
    with rasterio.open(default_mrt_file_path) as src:
        G = make_walking_network_graph(src, date_time_string)
    get_route(brickyard, psych_north, date_time_string)

