from rest_framework.serializers import ModelSerializer
from .models import RouteRequest


class RouteRequestSerializer(ModelSerializer):
    class Meta:
        model = RouteRequest
        fields = ('startPoint', 'endPoint', 'dateTime')

