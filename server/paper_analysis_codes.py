import rasterio
import pandas as pd
from routing import make_walking_network_graph
from datetime import datetime, timezone, timedelta

if __name__ == '__main__':
    selected_dates = [
        [2023, 10, 15],
        [2023, 11, 12],
        [2023, 12, 15],
        [2024, 1, 16],
        [2024, 2, 15],
        [2024, 3, 13],
        [2024, 4, 15],
        [2024, 5, 15],
        [2024, 6, 15],
        [2024, 7, 17],
        [2024, 8, 15],
        [2024, 9, 16],
    ]

    for dates in selected_dates:
        target_date = datetime(dates[0], dates[1], dates[2], tzinfo=timezone.utc)
        target_date = target_date.replace(minute=0, second=0, microsecond=0)


        for hour in range(24):
            target_date_ts = pd.Timestamp(target_date).replace(hour=hour)
            timekey = target_date_ts.strftime('%Y-%m-%d-%H00')
            mean_radiant_temp = rasterio.open('historical_mrt_data/{0}_mrt.tif'.format(timekey))
            make_walking_network_graph(mean_radiant_temp, timekey)
