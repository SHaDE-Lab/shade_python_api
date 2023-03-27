from rest_framework.serializers import ModelSerializer
from .models import MRT

class MRTSerializer(ModelSerializer):
    class Meta:
        model = MRT
        fields = (
            'date_time'
        )