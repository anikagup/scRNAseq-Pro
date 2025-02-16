# Use Python 3.10 as the base image
FROM python:3.10

# Set the working directory
WORKDIR /app

# Copy dependencies first
COPY requirements.txt .

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy source code
COPY src/ ./src
COPY main.py .

# Set the default command
CMD ["python", "src/main.py"]
