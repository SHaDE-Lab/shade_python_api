from solweig_run import run_solweig
from arcgisroute import make_walking_network_graph
from datetime import date
import sys
import pandas as pd

def main():
    # if the arg is daily the run a daily job
    if sys.argv[1] == 'daily':
        run_solweig_daily()
    elif sys.argv[1] == 'build_buffer':
        # run the solweig job for the next 72 hours
        for day in range(3):
            today = date.today()
            today = today.replace(day=today.day + day)
            for hour in range(24):
                # gets the timestamp but zeros out the minutes, seconds, and microseconds
                today_ts = pd.Timestamp(today).replace(hour=hour, minute=0, second=0, microsecond=0)
                run_solweig_hourly(today_ts)
    else:
        # Get the Date
        today = date.today()
        # gets the timestamp but zeros out the minutes, seconds, and microseconds
        today_ts = pd.Timestamp(today).replace(minute=0, second=0, microsecond=0)
        run_solweig_hourly(today_ts)

def run_solweig_daily():
    # run solweig for all of the hours in the day
    for hour in range(24):
        # Get the Date
        today = date.today()
        # gets the timestamp but zeros out the minutes, seconds, and microseconds
        today_ts = pd.Timestamp(today).replace(hour=hour, minute=0, second=0, microsecond=0)

        run_solweig_hourly(today_ts)


def run_solweig_hourly(today_ts):
    # run solweig
    run_solweig(today_ts)

    datetime = datetime.datetime.now()
    # 0 out the minutes and seconds
    datetime = datetime.replace(minute=0, second=0, microsecond=0)
    # get the newly generated mean radiant temperature raster
    mean_radiant_temp = rasterio.open('output/{0}_mean_radiant_temperature.tif'.format(datetime))

    # generate the graph
    make_walking_network_graph(mean_radiant_temp, datetime)

    # clean up any existing files in the output directory
    for file in os.listdir('output'):
        # delete a file if the time stamp is before the current time
        if file.endswith('.tif') and file < datetime.datetime.now().strftime('%Y-%m-%d_%H:%M'):
            os.remove(os.path.join('output', file))


if __name__ == '__main__':
    main()