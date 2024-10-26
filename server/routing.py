import os

import rasterio
from pyproj import Transformer
import time
import osmnx as ox
from shapely.geometry import LineString
import geojson
import geopandas as gpd

import pandas as pd
from datetime import datetime, timezone

default_mrt_file_path = 'output/2024-04-26-2100_mrt.tif'
default_graph_path = 'output/2024-04-26-2100_graph_networked.graphml' 
default_geopackage_path = 'Maps/walkOnlyGraph.gpkg'
def get_data_from_raster(lat, long, raster):
    # gets the data from the raster at the given lat/long coordinates (should be in same crs as raster)
    sample = raster.sample([(lat, long)])
    return next(sample, 0)

def readFromGeopackage(file_path):
    # Read the GeoPackage file into GeoDataFrames
    file_path = os.path.join(os.getcwd(), file_path)
    gdf_nodes = gpd.read_file(file_path, layer='nodes')
    gdf_edges = gpd.read_file(file_path, layer='edges')

    # Set the index of gdf_nodes to 'osmid'
    gdf_nodes.set_index('osmid', inplace=True)

    # Set the multi-index of gdf_edges to ['u', 'v', 'key']
    gdf_edges.set_index(['u', 'v', 'key_'], inplace=True)

    # Convert the GeoDataFrames to a NetworkX graph
    G = ox.graph_from_gdfs(gdf_nodes=gdf_nodes, gdf_edges=gdf_edges)
    # weird work around for the edge data being a weird character on one edge in the graph
    edge_data_key  = 'lanes'
    for u, v, key, data in G.edges(keys=True, data=True):
        if edge_data_key in data and data[edge_data_key] == '':
            data[edge_data_key] = 2
            break

    return G

def make_walking_network_graph(mean_radiant_temp, date_time_string):
    granularity = 1
   
    G = readFromGeopackage(default_geopackage_path)
    G = ox.project_graph(G, to_crs=mean_radiant_temp.crs)
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
                mrt_at_point = get_data_from_raster(coords[0], coords[1], mean_radiant_temp)
                total_mrt += max(mrt_at_point, 0)
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
                mrt_at_point = get_data_from_raster(coords[0], coords[1], mean_radiant_temp)
                total_mrt += max(mrt_at_point, 0)

        mean_mrt_value = total_mrt / (num_samples + 1)

        # Add the custom MRT attribute and cost attribute to the edge
        try:
            if isinstance(mean_mrt_value, float):
                data['mrt'] = mean_mrt_value
                data['cost'] = total_mrt
            else:
                data['mrt'] = mean_mrt_value[0]
                data['cost'] = total_mrt[0]
        except:
            print(mean_mrt_value)

    ox.save_graphml(G, 'historical_mrt_data/{0}_graph_{1}.graphml'.format(date_time_string, 'networked'))
    return G


def get_route(start_coord, stop_coord, date_time_string, weight='cost'):
    # if the graph file exists, load it, otherwise make it
    graph_path = 'historical_mrt_data/{0}_graph_{1}.graphml'.format(date_time_string, 'networked')
    attribute_types = {'cost': float, 'oneway': str}

    if os.path.exists(graph_path):
        G = ox.load_graphml(graph_path, edge_dtypes=attribute_types)
    else:
        G = ox.load_graphml(default_graph_path, edge_dtypes=attribute_types)
        print(f'using default graph for {date_time_string}', flush=True)
    G = ox.utils_graph.get_undirected(G) 
    ox.plot_graph(G)
    path_time = time.time()
    transformer = Transformer.from_crs("EPSG:4326", G.graph['crs'], always_xy=True)
    # Find the nearest network nodes to the origin and destination
    long, lat = transformer.transform(start_coord[0], start_coord[1])
    orig_node = ox.distance.nearest_nodes(G, long, lat)
    long, lat = transformer.transform(stop_coord[0], stop_coord[1])
    dest_node = ox.distance.nearest_nodes(G, long, lat)
    # Calculate the route using Dijkstra's algorithm (shortest path)
    route = ox.routing.shortest_path(G, orig_node, dest_node, weight=weight)
    # route is list of node IDs constituting the shortest path
    print(f'path {route} found in {time.time() - path_time}')  

    # calculate stats from the route
    statistics = calculate_statistics(G, route)

    geojson, mrt = convert_to_geoJSON(G, route)

    return statistics, geojson, mrt

def convert_to_geoJSON(G, route):
    # for each node in the route, get the lat/long
    route_coords = []
    mrt_values = []
    for node in route:
        node_coords = ((G.nodes[node]['lon'], G.nodes[node]['lat']))
        route_coords.append(node_coords)
        if route.index(node) < len(route) - 1:
            # avg mrg value between the two nodes in the route
            edge_data = G.get_edge_data(node, route[route.index(node) + 1])[0]
            mrt_values.append(float(edge_data['mrt']))
    # convert the route to a geojson
    route_geojson = geojson.LineString(route_coords)
    return route_geojson, mrt_values

def calculate_statistics(G, route):
    if route and len(route) < 2:
        print("Route is too short to calculate statistics.")

    # calculate the length of the route
    length = 0.0
    # calculate the total mrt of the route
    mrt = 0.0
    weighted_sum = 0.0
    for i in range(len(route) - 1):
        source_node = route[i]
        target_node = route[i + 1]
        # Use .edges() to get edge data for the current pair of nodes
        edge = G.get_edge_data(source_node, target_node)[0]
        length += float(edge['length'])
        mrt += float(edge['mrt'])
        weighted_sum += float(edge['mrt']) * float(edge['length'])
    average_mrt = weighted_sum / length
    # return stats dictionary
    return {'length': length, 'mrt': mrt, 'average_mrt': average_mrt}

if __name__ == '__main__':
    G = readFromGeopackage(default_geopackage_path)
    ox.plot_graph(G)
    # brickyard = (-111.93952587328305, 33.423795079832)  # brickyard in wgs84 (long, lat)
    # psych_north = (-111.92961401637315, 33.42070688780706)  # psych north in wgs84 (long, lat)
    # date_time_string = '2024-06-15-1100'
    # res1 = get_route(brickyard, psych_north, date_time_string, 'cost')
    # res2 = get_route(brickyard, psych_north, date_time_string, 'length')

    selected_dates = [
        [2023, 10, 15],
        [2023, 11, 12],
        [2023, 12, 15],
        [2024, 1, 16],
        [2024, 2, 15],
        [2024, 3, 13],
        [2024, 4, 15],
        [2024, 5, 15],
        [2024, 6, 15],
        [2024, 7, 17],
        [2024, 8, 15],
        [2024, 9, 16],
    ]

    answer = []

    locations = pd.read_csv("final_landmarks.csv")
    for i in range(1, len(locations)):
        src = (locations.iloc[i]["lon"], locations.iloc[i]["lat"])
        for j in range(i + 1, len(locations)):
            dest = (locations.iloc[j]["lon"], locations.iloc[j]["lat"])

            for dates in selected_dates:
                target_date = datetime(dates[0], dates[1], dates[2], tzinfo=timezone.utc)
                target_date = target_date.replace(minute=0, second=0, microsecond=0)

                for hour in range(24):
                    target_date_ts = pd.Timestamp(target_date).replace(hour=hour)
                    timekey = target_date_ts.strftime('%Y-%m-%d-%H00')

                    res1 = get_route(src, dest, timekey, 'cost')
                    res2 = get_route(src, dest, timekey, 'length')

                    answer.append({
                        "source_name": locations.iloc[i]["name"],
                        "source_lon": locations.iloc[i]["lon"],
                        "source_lat": locations.iloc[i]["lat"],

                        "dest_name": locations.iloc[j]["name"],
                        "dest_lon": locations.iloc[j]["lon"],
                        "dest_lat": locations.iloc[j]["lat"],

                        "optimal_length": round(res1[0]["length"], 5),
                        "optimal_mrt": round(res1[0]["mrt"], 5),
                        "optimal_avg_mrt": round(res1[0]["average_mrt"], 5),

                        "shortest_path_length": round(res2[0]["length"], 5),
                        "shortest_path_mrt": round(res2[0]["mrt"], 5),
                        "shortest_path_avg_mrt": round(res2[0]["average_mrt"], 5)
                    })
    
    pd.DataFrame(answer).to_csv("random_cool_routes_optimal_vs_shortest.csv", sep=",", index=False)