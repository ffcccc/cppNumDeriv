# Dockerfile used by .devcontainer/devcontainer.json

FROM ubuntu:24.04

# …other build dependencies…

# eigen is header‑only; the apt package is 3.3.x (3.3.9 in focal-backports)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libeigen3-dev \
        mc wget \
    && rm -rf /var/lib/apt/lists/*

# alternatively, fetch a specific release if you really need 3.3.9
# RUN wget https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz \
#     && tar xzf eigen-3.3.9.tar.gz -C /usr/local/include --strip-components=1 \
#     && rm eigen-3.3.9.tar.gz

# …rest of your Dockerfile…