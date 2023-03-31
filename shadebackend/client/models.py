from django.db import models


# Create your models here.
# class MRT(models.Model):
#    date_time = models.DateTimeField(max_length=255)

#    def __str__(self):
#        return self.date_time


class RouteRequest(models.Model):
    # startPoint: {long, lat} in WGS84
    startPoint: models.TextField(max_length=4095)
    # endPoint: {long, lat} in WGS84
    endPoint: models.TextField(max_length=4095)
    dateTime = models.DateTimeField(max_length=255)

    class Meta:
        app_label = 'client'


class RouteResponse(models.Model):
    # path: [{long, lat} in NAD83]
    path: models.TextField(max_length=4095)
