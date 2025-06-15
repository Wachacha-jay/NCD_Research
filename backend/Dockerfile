
# Simple Dockerfile for the Research Agent Backend (Flask)
FROM python:3.11

WORKDIR /app

# Install dependencies
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Copy code
COPY . .

# Expose port
EXPOSE 5000

CMD ["python", "app.py"]
