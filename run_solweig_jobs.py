#!/usr/bin/env python3

from solweig_run import run_solweig
from routing import make_walking_network_graph
from datetime import datetime, timezone, timedelta
import sys
import pandas as pd
import rasterio
import os

def main():
    # if the arg is daily the run a daily job
    if len(sys.argv) > 1:
        print("running now")
        if sys.argv[1] == 'daily':
            run_solweig_daily()
        if sys.argv[1] == 'build_buffer':
            # run the solweig job for the next 72 hours
            run_solweig_buildup()
        elif sys.argv[1] == "hourly":
            # add an hour to today
            today = datetime.now(timezone.utc) + timedelta(hours=1)
            print(today)
            # gets the timestamp but zeros out the minutes, seconds, and microseconds
            today_ts = pd.Timestamp(today).replace(minute=0, second=0, microsecond=0)
            run_solweig_hourly(today_ts)
    else:
        # Get the Date
        today = datetime.now(timezone.utc)
        print(today)
        # gets the timestamp but zeros out the minutes, seconds, and microseconds
        today_ts = pd.Timestamp(today).replace(minute=0, second=0, microsecond=0)
        run_solweig_hourly(today_ts)
    file_cleanup()

def run_solweig_buildup():
    for day in range(3):
        today = datetime.now(timezone.utc)
        today = today.replace(day=today.day + day)
        for hour in range(24):
            # gets the timestamp but zeros out the minutes, seconds, and microseconds
            today_ts = pd.Timestamp(today).replace(hour=hour, minute=0, second=0, microsecond=0)
            run_solweig_hourly(today_ts)

def run_solweig_daily():
    # run solweig for all of the hours in the day
    for hour in range(24):
        # Get the Date
        today = datetime.now(timezone.utc)
        # gets the timestamp but zeros out the minutes, seconds, and microseconds
        today_ts = pd.Timestamp(today).replace(hour=hour, minute=0, second=0, microsecond=0)

        run_solweig_hourly(today_ts)


def run_solweig_hourly(today_ts):
    print("running solweig for {0}".format(today_ts))
    timekey = today_ts.strftime('%Y-%m-%d-%H00')
    # run solweig
    run_solweig(today_ts)

    # get the newly generated mean radiant temperature raster
    # the timestamp is usually in UTC so we need to convert it to remove the extra info
    mean_radiant_temp = rasterio.open('output/{0}_mrt.tif'.format(timekey))
    print("finished making mrt for {0}".format(timekey))

    # generate the graph
    make_walking_network_graph(mean_radiant_temp, timekey)
    # Add the following lines at the end to create a signal file
    with open('/app/.reboot_done', 'w') as signal_file:
        signal_file.write('Reboot job completed\n')


def file_cleanup():
    # clean up any existing files in the output directory
    for file in os.listdir('output'):
        # delete a file if the time stamp is before the current time
        if file.endswith('.tif') and file < datetime.now().strftime('%Y-%m-%d-%H00'):
            print('deleting {0}'.format(file))
            # TODO DUMP TO LONG TERM STORAGE
            os.remove(os.path.join('output', file))


if __name__ == '__main__':
    main()