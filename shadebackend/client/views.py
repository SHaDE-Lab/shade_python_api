from django.http import HttpResponse
from django.shortcuts import render
from rest_framework.mixins import (
    CreateModelMixin, ListModelMixin, RetrieveModelMixin, UpdateModelMixin
)
from rest_framework.viewsets import GenericViewSet
# import arcgisroute
import json


# Create your views here.
def get_route(request, coords):
    load = json.loads(coords)
    start = tuple(load['startPoint'])
    end = tuple(load['endPoint'])
    date_time = load['dateTime']

    # path = arcgisroute.getRoute(start, end, date_time)

    # line below is temporary and purely for testing purposes, delete it when model is hooked up
    path = [[1.1, 2.2], [3.3, 4.4], [5.5, 6.6], [7.7, 8.8]]

    response = json.dumps(path)
    return HttpResponse(response)


# class MRTViewSet(GenericViewSet,
#                 CreateModelMixin,
#                 RetrieveModelMixin,
#                 UpdateModelMixin,
#                 ListModelMixin):
#    serializer_class = MRTSerializer
#    queryset = MRTSerializer.objects.all()
    
