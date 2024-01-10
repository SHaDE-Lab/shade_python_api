from datetime import datetime, timezone, timedelta
from flask import Flask, jsonify, request, Response
from routing import get_route
from flask_cors import CORS
import json
import pandas as pd
from flask_apscheduler import APScheduler
from threading import Thread

from run_solweig_jobs import run_solweig_hourly, run_solweig_buildup

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes
scheduler = APScheduler()

@app.route('/api/route', methods=['GET'])
def route_handler():
    try:
        # Convert the JSON data from the URL to a Python dictionary
        # Get the JSON data from the query string
        json_data = request.args.get('json_data')

        # Parse the JSON data into a Python dictionary
        data_dict = json.loads(json_data)

        # Extract required values from the dictionary
        start_lat = float(data_dict['startPoint'][0])
        start_long = float(data_dict['startPoint'][1])
        end_lat = float(data_dict['endPoint'][0])
        end_long = float(data_dict['endPoint'][1])
        start = (start_long, start_lat)
        end = (end_long, end_lat)
        date_time = data_dict['dateTime']
        print(f"Received route request from {start} to {end} at {date_time}")
        kml, stats, geojson = get_route(start, end, date_time)
        return jsonify({'kml': kml, 'stats': stats, 'geojson': json.dumps(geojson)})
    except Exception as e:
        # Handle errors appropriately
        print(e)
        return jsonify({'error': str(e)})


# Define a route to serve the KML file
@app.route('/get_kml', methods=['GET'])
def get_kml():
    print("Received kml request")
    try:
        # Convert the JSON data from the URL to a Python dictionary
        # Get the JSON data from the query string
        json_data = request.args.get('json_data')

        # Parse the JSON data into a Python dictionary
        data_dict = json.loads(json_data)

        # Extract required values from the dictionary
        start_lat = float(data_dict['startPoint'][0])
        start_long = float(data_dict['startPoint'][1])
        end_lat = float(data_dict['endPoint'][0])
        end_long = float(data_dict['endPoint'][1])
        start = (start_long, start_lat)
        end = (end_long, end_lat)
        date_time = data_dict['dateTime']
        print(f"Received kml request from {start} to {end} at {date_time}")
        kml_data, _ = get_route(start, end, date_time)
        print(len(kml_data))
        return Response(kml_data, mimetype='application/vnd.google-earth.kml+xml')
    except Exception as e:
        # Handle errors appropriately
        print(e)
        return jsonify({'error': str(e)})


@app.route('/')
def hello_world():
    return {'message': 'Hello, World!'}

def run_after_startup():
    print("running solweig buildup")
    run_solweig_buildup()
    print("finished running solweig buildup")

if __name__ == '__main__':
    print('starting server')
    scheduler.api_enabled = True
    scheduler.init_app(app)

    # Start the scheduler in a separate thread
    scheduler_thread = Thread(target=scheduler.start)
    scheduler_thread.start()

    # Start a new thread for run_after_startup
    startup_thread = Thread(target=run_after_startup)
    startup_thread.start()
    app.run(host="0.0.0.0",port=5000,debug=True)


@scheduler.task('interval', id='my_job', minutes=60)
def my_job():
    hour_from_now = datetime.now(timezone.utc) + timedelta(hours=1)
    hour_from_now = pd.Timestamp(hour_from_now).replace(minute=0, second=0, microsecond=0)
    run_solweig_hourly(hour_from_now)