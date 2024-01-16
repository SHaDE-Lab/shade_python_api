from flask import Flask, jsonify, request, Response
from routing import get_route
from flask_cors import CORS
import json

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

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


if __name__ == '__main__':
    print('starting server')
    app.run(host="0.0.0.0",port=5000,debug=True)
