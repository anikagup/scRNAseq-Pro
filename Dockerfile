# Use an official Python runtime as a parent image
FROM python:3.10

# Set the working directory
WORKDIR /app

# Copy requirements and install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Install Shiny for Python
RUN pip install shiny matplotlib scanpy

# Copy all files
COPY . .

# Create a directory for serving images
RUN mkdir -p /app/www

# Expose the app port
EXPOSE 8000

# Run the Shiny app
CMD ["shiny", "run", "--port", "8000", "app.py"]
