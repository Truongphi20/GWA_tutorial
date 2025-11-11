FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Build argument for Docker group ID
ARG DOCKER_GID=137

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    unzip \
    git \
    graphviz \
    openjdk-17-jdk \
    docker.io \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user for VS Code
RUN useradd -ms /bin/bash devuser

# Modify existing docker group GID and add devuser
RUN groupmod -g $DOCKER_GID docker \
    && usermod -aG docker devuser

# Install Nextflow CLI
RUN curl -s https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/ \
    && chmod +x /usr/local/bin/nextflow

USER devuser
WORKDIR /workspace
