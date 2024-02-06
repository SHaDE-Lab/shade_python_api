import networkx as nx
import geopandas as gp
import osmnx as ox
import pandas as pd
from shapely.geometry import LineString

def convert_shp2graph(p, make_G_bidi = True, name='unamed'):
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
    if 'length' not in gdf.columns:
        raise Exception('Shapefile is invalid: length not in attributes:\n{}'.format(gdf.columns))

    if  not gdf.geometry.map(lambda x: type(x) ==  LineString).all():
        s_invalid_geo = gdf.geometry[gdf.geometry.map(lambda x: type(x) ==  LineString)]
        raise Exception('Shapefile is invalid: geometry not all linestring \n{}'.format(s_invalid_geo))
   
    # Compute the start- and end-position based on linestring 
    gdf['Start_pos'] = gdf.geometry.apply(lambda x: x.coords[0])
    gdf['End_pos'] = gdf.geometry.apply(lambda x: x.coords[-1])
    
    # Create Series of unique nodes and their associated position
    s_points = gdf.Start_pos.append(gdf.End_pos).reset_index(drop=True)
    s_points = s_points.drop_duplicates()   
#     log('GeoDataFrame has {} elements (linestrings) and {} unique nodes'.format(len(gdf),len(s_points)))
    
    # Add index of start and end node of linestring to geopandas DataFrame
    df_points = pd.DataFrame(s_points, columns=['Start_pos'])
    df_points['FNODE_'] = df_points.index
    gdf = pd.merge(gdf, df_points, on='Start_pos', how='inner')

    df_points = pd.DataFrame(s_points, columns=['End_pos'])
    df_points['TNODE_'] = df_points.index
    gdf = pd.merge(gdf, df_points, on='End_pos', how='inner')
    
    # Bring nodes and their position in form needed for osmnx (give arbitrary osmid (index) despite not osm file)
    df_points.columns = ['pos', 'osmid'] 
    df_points[['x', 'y']] = df_points['pos'].apply(pd.Series)
    df_node_xy = df_points.drop('pos', 1)
    
    # Create Graph Object
    G = nx.MultiDiGraph(name=name, crs=gdf.crs)
    
    # Add nodes to graph
    for node, data in df_node_xy.T.to_dict().items():
        G.add_node(node, **data)
        
    # Add edges to graph
    for i, row  in gdf.iterrows():
        dict_row  = row.to_dict()
        if 'geometry' in dict_row: del dict_row['geometry']
        G.add_edge(u=dict_row['FNODE_'], v=dict_row['TNODE_'], **dict_row)
        
    if make_G_bidi:
        gdf.rename(columns={'Start_pos': 'End_pos',
                   'End_pos': 'Start_pos',
                   'FNODE_': 'TNODE_', 
                   'TNODE_': 'FNODE_', }, inplace=True)
        
        # Add edges to graph
        for i, row  in gdf.iterrows():
            dict_row  = row.to_dict()
            if 'geometry' in dict_row: del dict_row['geometry']
            G.add_edge(u=dict_row['FNODE_'], v=dict_row['TNODE_'], **dict_row)
        
#         G = G.to_undirected() # Some function in osmnx do not work anymore
        
    # Log information
#     log('Graph has been successfully generated /n {}'.format(nx.info(G)))
#     log('Show graph data structure EDGE'.format(G.get_edge_data(*list(G.edges())[0])))
#     log('Show graph data structure NO
    return G 