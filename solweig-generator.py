import os
import time
from datetime import datetime, timedelta, timezone
import pandas as pd
from run_solweig_jobs import run_solweig_hourly, run_solweig_buildup

def startup_function(output_folder):
    print("File generator started!")
    os.makedirs(output_folder, exist_ok=True)
    run_solweig_buildup()

def generate_file(output_folder):
    hour_from_now = datetime.now(timezone.utc) + timedelta(hours=1)
    hour_from_now = pd.Timestamp(hour_from_now).replace(minute=0, second=0, microsecond=0)
    print(f"Generating file for {hour_from_now}")
    run_solweig_hourly(hour_from_now)

def hourly_function(output_folder):
    while True:
        generate_file(output_folder)
        # Sleep for one hour (3600 seconds)
        time.sleep(3600)

if __name__ == "__main__":
    # Specify the output folder (in the base directory)
    print("Starting up solweig generator...")
    output_folder = "output"

    # Call the startup function
    startup_function(output_folder)

    # Run the hourly function
    hourly_function(output_folder)
