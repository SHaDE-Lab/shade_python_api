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

# Hourly Jobs 

Running solweig and making the graphs are controlled by the `run_solweig_jobs.py` script.
There are 3 params that it can take in:
1. `hourly` this will run the solweig job and make the graphs for the current hour
2. `daily` this will run the solweig job and make the graphs for the current day
3. `build_buffer` this will run the job and make the graphs for the next 72 hours (3 days)

When having a cold start you should run `python run_solweig_jobs.py build_buffer` to make sure that the buffer is built up.
Then you can run `python run_solweig_jobs.py hourly` to run the job for the current hour.
Then you can run `python run_solweig_jobs.py daily` to run the job for the current day.

For automation purposes you could have a job that runs the `hourly` job every hour and the `daily` job every day at midnight.