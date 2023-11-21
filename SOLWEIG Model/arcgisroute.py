import os

import rasterio
from pyproj import Transformer
import time
import osmnx as ox
from shapely.geometry import LineString
import simplekml
import geojson

default_mrt_file_path = 'output/2023-4-8-2100_mrt.tif'
default_graph_path = 'output/2023-4-8-2100_graph_networked.graphml'
def convert_to_pixel(lat, long, raster):
    # Get the pixel coordinates from the lat long
    return raster.index(lat, long)

def make_walking_network_graph(mean_radiant_temp, date_time_string):
    # bounding box of tempe campus
    north = 33.42796506379531
    west = -111.94015084151006

    south = 33.41592636331673
    east = -111.92634308211977
    granularity = 1

    ox.settings.use_cache = True
    G = ox.graph_from_bbox(north, south, east, west, network_type='walk')
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
            total_mrt = 0.0
            for i in range(num_samples + 1):
                alpha = i / num_samples
                # Interpolate the point along the LineString
                point = edge_line.interpolate(alpha, normalized=True)
                coords = (point.x, point.y)
                u = convert_to_pixel(coords[0], coords[1], mean_radiant_temp)
                total_mrt += max(mrt_data[u], 0)
        else:
            # If no 'geometry' key, interpolate between u and v
            total_mrt = 0.0
            u_coords = (G.nodes[u]['x'], G.nodes[u]['y'])
            v_coords = (G.nodes[v]['x'], G.nodes[v]['y'])
            edge_length = ox.distance.euclidean(u_coords[1], u_coords[0], v_coords[1], v_coords[0])
            num_samples = int(edge_length * granularity)
            num_samples = num_samples if num_samples > 0 else 1
            for i in range(num_samples + 1):
                alpha = i / num_samples
                coords = (u_coords[0] + alpha * (v_coords[0] - u_coords[0]), v_coords[1] + alpha * (v_coords[1] - u_coords[1]))
                u = convert_to_pixel(coords[0], coords[1], mean_radiant_temp)
                total_mrt += max(mrt_data[u], 0)

        mean_mrt_value = total_mrt / (num_samples + 1)

        # Add the custom MRT attribute and cost attribute to the edge
        data['mrt'] = mean_mrt_value
        data['cost'] = total_mrt
    ox.save_graphml(G, 'output/{0}_graph_{1}.graphml'.format(date_time_string, 'networked'))
    return G


def get_route(start_coord, stop_coord, date_time_string):
    print('starting route')

    print('reading files')

    # if the graph file exists, load it, otherwise make it
    graph_path = 'output/{0}_graph_{1}.graphml'.format(date_time_string, 'networked')
    mrt_file_path = 'output/{0}_mrt.tif'.format(date_time_string)  # expected format is "2023-03-30_12:00"
    attribute_types = {'cost': float}

    if not os.path.exists(mrt_file_path):
        print('no mrt file found')
        mrt_file_path = default_mrt_file_path
        G = ox.load_graphml(default_graph_path, edge_dtypes=attribute_types)
    elif os.path.exists(graph_path):
        G = ox.load_graphml(graph_path, edge_dtypes=attribute_types)
        print('graph loaded')
    else:
        # Create a graph representing the raster cells and their connections
        with rasterio.open(mrt_file_path) as src:
            G = make_walking_network_graph(src, date_time_string)
            print('made graph')
    print('finding path')
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
    print(route)
    print('path found in {0}'.format(time.time() - path_time))
    # Convert the route to a GeoDataFrame# Plot the graph with the route highlighted
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

brickyard = (-111.93952587328305, 33.423795079832)  # brickyard in wgs84 (long, lat)
psych_north = (-111.92961401637315, 33.42070688780706)  # psych north in wgs84 (long, lat)
date_time_string = '2023-4-8-2100'

if __name__ == '__main__':
    get_route(brickyard, psych_north, date_time_string)
