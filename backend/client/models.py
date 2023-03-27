from django.db import models

# Create your models here.
class MRT(models.Model):
    date_time = models.DateTimeField(max_length=255)

    def __str__(self):
        return self.date_time
