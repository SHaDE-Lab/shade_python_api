# Use an official Python runtime as a parent image
FROM andrejreznik/python-gdal:py3.10.0-gdal3.2.3

# Set the working directory to /app
WORKDIR /app

# Copy only the requirements file to the working directory
COPY requirements.txt /app/

RUN pip install --upgrade pip

# Install any needed packages specified in requirements.txt
RUN pip install -r requirements.txt

# Copy the current directory contents into the container at /app
COPY . /app

# Make port 5000 available to the world outside this container

EXPOSE 5000

# Run app.py when the container launches
CMD ["python", "server.py"]
