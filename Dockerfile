# Use an official Python runtime as a parent image
FROM python:3.10

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Set the working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libomp-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install dependencies
COPY requirements.txt .

RUN pip install --upgrade pip setuptools wheel \
    && pip install --no-cache-dir -r requirements.txt \
    && pip install shiny matplotlib scanpy

#explicitly set container's internal path
ENV PATH="/usr/local/bin:${PATH}"

# Copy all files
COPY . .

# Create a directory for serving images
RUN mkdir -p /app/www

# Create directories needed by the app
RUN mkdir -p /app/www /app/uploads /app/processed_data /app/figures

# Expose the app port
EXPOSE 8000

# Run the Shiny app
CMD ["shiny", "run", "UI/anikaapp.py", "--host", "0.0.0.0", "--port", "8000"]

