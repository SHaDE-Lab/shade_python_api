import math
import os
import rasterio
from pyproj import Transformer
import time
import osmnx as ox
from shapely.geometry import LineString
import simplekml
import geojson
import osmnx as ox
import geopandas as gpd

default_mrt_file_path = 'output/default_mrt.tif'
default_graph_path = 'output/default_graph_networked.graphml'
geopackage_path = 'Maps/downloadedShapeFile.gpkg'
node_shapefile_path = 'Maps/downloadedShapeFile/nodes.shp'
edge_shapefile_path = 'Maps/downloadedShapeFile/edges.shp'
granularity = 1

def readFromShapeFile(node_filepath, edge_filepath):
    # Read the GeoPackage file into GeoDataFrames
    gdf_nodes = gpd.read_file(node_filepath, layer='nodes')
    gdf_edges = gpd.read_file(edge_filepath, layer='edges')

    # Set the index of gdf_nodes to 'osmid'
    gdf_nodes.set_index('osmid', inplace=True)

    # Set the multi-index of gdf_edges to ['u', 'v', 'key']
    gdf_edges.set_index(['u', 'v', 'key'], inplace=True)
    # Convert the GeoDataFrames to a NetworkX graph
    G = ox.graph_from_gdfs(gdf_nodes=gdf_nodes, gdf_edges=gdf_edges)
    return G

def readFromGeopackage(file_path):
    # Read the GeoPackage file into GeoDataFrames
    gdf_nodes = gpd.read_file(file_path, layer='nodes')
    gdf_edges = gpd.read_file(file_path, layer='edges')

    # Set the index of gdf_nodes to 'osmid'
    gdf_nodes.set_index('osmid', inplace=True)

    # Set the multi-index of gdf_edges to ['u', 'v', 'key']
    gdf_edges.set_index(['u', 'v', 'key'], inplace=True)

    # Convert the GeoDataFrames to a NetworkX graph
    G = ox.graph_from_gdfs(gdf_nodes=gdf_nodes, gdf_edges=gdf_edges)
    return G


def convert_to_pixel(lat, long, raster):
    # Get the pixel coordinates from the lat long
    return raster.index(lat, long)

def make_walking_network_graph(mean_radiant_temp, date_time_string):
    # bounding box of tempe campus
    # G = readFromGeopackage(geopackage_path)
    # G = readFromShapeFile(node_shapefile_path, edge_shapefile_path)
    North = 33.42943
    South = 33.409729
    East = -111.917768
    West = -111.941743

    # ox.settings.use_cache = True
    # G = ox.graph_from_bbox(North, South, East, West, network_type='walk')
    G = readFromGeopackage(geopackage_path)
    G = ox.project_graph(G)
    ox.plot_graph(G)
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
            mean_mrt_value = 10000
            total_mrt = 10000
        # Add the custom MRT attribute and cost attribute to the edge
        data['mrt'] = mean_mrt_value
        data['cost'] = total_mrt
    ox.save_graphml(G, 'output/{0}_graph_{1}.graphml'.format(date_time_string, 'networked'))
    return G


def get_route(start_coord, stop_coord, date_time_string):
    # if the graph file exists, load it, otherwise make it
    graph_path = 'output/{0}_graph_{1}.graphml'.format(date_time_string, 'networked')
    mrt_file_path = 'output/{0}_mrt.tif'.format(date_time_string)  # expected format is "2023-03-30-1200"
    attribute_types = {'cost': float, 'mrt': float}

    # if graph does not exist, use default graph
    if not os.path.exists(graph_path):
        print('Graph does not exist for {0}, using default graph'.format(graph_path))
        graph_path = default_graph_path
    G = ox.load_graphml(graph_path, edge_dtypes=attribute_types)

    path_time = time.time()
    transformer = Transformer.from_crs("EPSG:4326", G.graph['crs'], always_xy=True)
    # Find the nearest network nodes to the origin and destination
    long, lat = transformer.transform(start_coord[0], start_coord[1])
    orig_node = ox.distance.nearest_nodes(G, long, lat)
    long, lat = transformer.transform(stop_coord[0], stop_coord[1])
    dest_node = ox.distance.nearest_nodes(G, long, lat)
    # Calculate the route using Dijkstra's algorithm (shortest path)
    print(orig_node, dest_node)

    route = ox.routing.k_shortest_paths(G, orig_node, dest_node, k=1, weight='cost')
    # convert route from generator object
    route = list(route)[0]
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


# TODO - add function to determine if an edge is walk only
def isWalkOnly(edge):
    return False

if __name__ == '__main__':
    brickyard = (-111.93952587328305, 33.423795079832)  # brickyard in wgs84 (long, lat)
    psych_north = (-111.935056, 33.4148208)  # psych north in wgs84 (long, lat)
    make_walking_network_graph(rasterio.open(default_mrt_file_path), 'default')
    get_route(brickyard, psych_north, 'default')

