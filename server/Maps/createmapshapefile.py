import osmnx as ox
# TODO CHANGE TO REALLY BE THE BOUNDING BOX
north = 33.42796506379531
west = -111.94015084151006

south = 33.41592636331673
east = -111.92634308211977
granularity = 1

# ox.settings.use_cache = True
G = ox.graph_from_bbox(north, south, east, west, network_type='walk')

ox.plot_graph(G)
ox.save_graph_shapefile(G, filepath='downloadedShapeFile', encoding='utf-8')