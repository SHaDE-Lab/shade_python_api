import osmnx as ox
import geopandas as gpd
# TODO CHANGE TO REALLY BE THE BOUNDING BOX
North = 33.42943
South = 33.409729
East = -111.917768
West = -111.941743

# ox.settings.use_cache = True
G = ox.graph_from_bbox(North, South, East, West, network_type='walk')
G = ox.utils_graph.remove_isolated_nodes(G)
ox.plot_graph(G)

# saves graph as geopackage and shapefile
ox.save_graph_geopackage(G, encoding='utf-8', filepath='downloadedShapeFile.gpkg')
ox.save_graph_shapefile(G, filepath='downloadedShapeFile', encoding='utf-8')

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