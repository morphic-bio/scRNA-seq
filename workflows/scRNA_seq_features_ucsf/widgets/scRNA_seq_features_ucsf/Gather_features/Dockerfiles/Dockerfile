FROM python:3.11.13-slim-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends build-essential && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir "numpy<2.0" scanpy pyroe

COPY ./gather_features.py /usr/local/bin/gather_features.py
RUN chmod +x /usr/local/bin/gather_features.py
