import rasterio as rio #to handle rasters
import os
import fnmatch
import solweig_mrt as sol #the main solweig script which imports all the otehr needed scripts and functions
#from functools import partial # uncomment for  parallel computing
#import multiprocess as mp # uncomment for  parallel computing
import warnings
warnings.filterwarnings('ignore')


svf_walls_paths= '../sample_data/tempe_camp_svfs/' #path with all input rasters (DSM,DEM, CDSM....)
out_path = '../sample_results/'

#put all the paths in a list
files_list = [svf_walls_paths+files for files in os.listdir(svf_walls_paths)]

