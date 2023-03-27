from django.shortcuts import render
from rest_framework.mixins import (
    CreateModelMixin, ListModelMixin, RetrieveModelMixin, UpdateModelMixin
)
from rest_framework.viewsets import GenericViewSet
from .models import MRT
from .serializers import MRTSerializer

# Create your views here.
class MRTViewSet(GenericViewSet,
                 CreateModelMixin,
                 RetrieveModelMixin,
                 UpdateModelMixin,
                 ListModelMixin):
    serializer_class = MRTSerializer
    queryset = MRTSerializer.objects.all()
    
