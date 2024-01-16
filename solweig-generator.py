import time
from datetime import datetime, timedelta, timezone
import pandas as pd
from run_solweig_jobs import run_solweig_hourly, run_solweig_buildup, run_solweig_daily
import schedule
def startup_function(output_folder):
    print("File generator started!", flush=True)
    run_solweig_buildup()
    print("Finished startup function", flush=True)

def daily_function(output_folder):
    print("Running daily function", flush=True)
    run_solweig_daily(datetime.now(timezone.utc) + timedelta(days=3))

def hourly_function(output_folder):
    hour_from_now = datetime.now(timezone.utc) + timedelta(hours=1)
    hour_from_now = pd.Timestamp(hour_from_now).replace(minute=0, second=0, microsecond=0)
    print(f"Generating file for {hour_from_now}")
    run_solweig_hourly(hour_from_now)

if __name__ == "__main__":
    # Specify the output folder (in the base directory)
    print("Starting up solweig generator...", flush=True)
    output_folder = "output"

    # Call the startup function
    startup_function(output_folder)

    # Schedule the job to run every hour
    schedule.every().hour.do(hourly_function(output_folder))
    schedule.every().day.at("00:00").do((daily_function(output_folder)))
    while True:
        # Run pending jobs
        print("Checking for jobs", flush=True)
        schedule.run_pending()
        time.sleep(60*60)
