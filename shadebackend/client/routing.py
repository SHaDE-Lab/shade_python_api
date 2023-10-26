import os

from pyproj import Transformer
import networkx as nx
import rasterio
import time
import osmnx as ox

def convert_to_pixel(lat, long, raster):
    # Get the pixel coordinates from the lat long
    return raster.index(lat, long)


def make_walking_network_graph(mean_radiant_temp, date_time_string):
    # bounding box of tempe campus
    north = 33.42796506379531
    west = -111.94015084151006

    south = 33.41592636331673
    east = -111.92634308211977

    ox.settings.use_cache = True
    G = ox.graph_from_bbox(north, south, east, west, network_type='walk')
    G = ox.project_graph(G)
    ox.plot_graph(G)
    mrt_data = mean_radiant_temp.read(1)
    for u, v, data in G.edges(data=True):
        # Get the coordinates of the edge's nodes
        u_coords = G.nodes[u]['x'], G.nodes[u]['y']
        v_coords = G.nodes[v]['x'], G.nodes[v]['y']
        # Convert coordinates to indices based on your grid or array
        u = convert_to_pixel(u_coords[0], u_coords[1], mean_radiant_temp)
        v = convert_to_pixel(v_coords[0], v_coords[1], mean_radiant_temp)
        # Calculate the mean MRT value along the edge
        mean_mrt_value = (max(mrt_data[u], 0) + max(mrt_data[v], 0)) / 2.0
        edge_cost = mean_mrt_value
        # Add the custom MRT attribute and cost attribute to the edge
        data['mrt'] = mean_mrt_value
        data['cost'] = edge_cost
    ox.save_graphml(G, 'output/{0}_graph_{1}.graphml'.format(date_time_string, 'networked'))
    return G


def calculate_route(start_coord, stop_coord, date_time_string):
    print('starting route')

    print('reading files')

    # if the graph file exists, load it, otherwise make it
    graph_path = 'output/{0}_graph_{1}.graphml'.format(date_time_string, 'networked')
    if os.path.exists(graph_path):
        attribute_types = {'cost': float}
        G = ox.load_graphml(graph_path, edge_dtypes=attribute_types)
        print('graph loaded')
    else:
        mrt_file_path = 'output/{0}_mrt.tif'.format(date_time_string)  # expected format is "2023-03-30-1200"
        print(mrt_file_path)
        print(os.path.exists(mrt_file_path))
        print(os.path.exists('/output/2023-4-8-2100_mrt.tif'))
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
    # Calculate the route using NetworkX's shortest_path function
    route = nx.shortest_path(G, orig_node, dest_node, weight='cost')
    print('path found in {0}'.format(time.time() - path_time))
    # Convert the route to a GeoDataFrame# Plot the graph with the route highlighted
    ox.plot_graph_route(G, route)
    print(route)
    return route


# brickyard = (-111.93952587328305, 33.423795079832)  # brickyard in wgs84 (long, lat)
# psych_north = (-111.92961401637315, 33.42070688780706)  # psych north in wgs84 (long, lat)
# calculate_route(brickyard, psych_north, '2023-4-8-2100')
