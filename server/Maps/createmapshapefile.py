import osmnx as ox
import geopandas as gpd
# TODO CHANGE TO REALLY BE THE BOUNDING BOX
north = 33.42796506379531
west = -111.94015084151006

south = 33.41592636331673
east = -111.92634308211977
granularity = 1

# ox.settings.use_cache = True
G = ox.graph_from_bbox(north, south, east, west)

ox.plot_graph(G)

ox.save_graph_geopackage(G, encoding='utf-8', filepath='downloadedShapeFile.gpkg')

# # read in the downloadedgeopackage
# load the shapefile
nodes = gpd.read_file('downloadedShapeFile/nodes.shp')
edges = gpd.read_file('downloadedShapeFile/edges.shp')

# Set the index of gdf_nodes to 'osmid'
nodes.set_index('osmid', inplace=True)

# Set the multi-index of gdf_edges to ['u', 'v', 'key']
edges.set_index(['u', 'v', 'key'], inplace=True)

# Specify the file path to the GeoPackage file
file_path = 'downloadedShapeFile.gpkg'

# Read the GeoPackage file into GeoDataFrames
gdf_nodes = gpd.read_file(file_path, layer='nodes')
gdf_edges = gpd.read_file(file_path, layer='edges')


# Convert the GeoDataFrames to a NetworkX graph
G = ox.graph_from_gdfs(gdf_nodes=nodes, gdf_edges=edges)


ox.plot_graph(G)