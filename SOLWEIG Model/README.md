# SOLWEIG-Implementation
scripts for running solweig in jupyter notebooks

Import the solweig_mrt.py into a jupyter notebook with "from solweig_mrt import*"

# Routing 
ArcGisRoute.py contains all the routing logic

# Solweig
To run solweig just run the solweig_run.py script with `python solweig_run.py`

# Server
We have switched to Flask. make sure its installed and then you can run the server.py file.
When adding api routes add it to this file try to abstract most of the logic to its own file

Make sure flask is installed `pip install flask`
You can run the server by cding to the directory and running `python server.py`