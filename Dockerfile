# Set base image (host OS)
FROM python:3.10

# By default, listen on port 5000
EXPOSE 5000/tcp

# Set the working directory in the container
WORKDIR /app

# Copy the dependencies file to the working directory
COPY requirements.txt .

# Install any dependencies
RUN apt-get update

RUN apt-get install -y libgdal-dev

RUN pip install GDAL==3.6.4

RUN pip install -r requirements.txt

# Copy the content of the local src directory to the working directory
COPY server.py .

# Specify the command to run on container start
CMD [ "python", "./app.py" ]