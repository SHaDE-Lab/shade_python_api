from flask import Flask, jsonify
from arcgisroute import get_route

app = Flask(__name__)

@app.route('/api/getRoute/<start>/<end>/<date_time>')
def get_route_api(start, end, date_time):
    kml, stats = get_route(start, end, date_time)
    return jsonify({'kml': kml, 'stats': stats})

@app.route('/')
def hello_world():
    return {'message': 'Hello, World!'}

if __name__ == '__main__':
    print('starting server')
    app.run(debug=True)