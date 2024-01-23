# Use an official Python runtime as a parent image
FROM python:3.10

# Install GDAL dependencies
RUN apt-get update \
    && apt-get install -y libgdal-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory to /app
WORKDIR /app

# Copy only the requirements file to the working directory
COPY requirements.txt /app/

# Install any needed packages specified in requirements.txt
RUN pip install -r requirements.txt

# Copy the current directory contents into the container at /app
COPY . /app

# Make port 5000 available to the world outside this container

EXPOSE 5000

# Run app.py when the container launches
CMD ["python", "server.py"]
