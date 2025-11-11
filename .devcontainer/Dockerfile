FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    curl \
    unzip \
    git \
    graphviz \
    openjdk-17-jdk \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow
    
# Setup docker
RUN apt-get update && apt-get install -y docker.io

# Create non-root user for VS Code
RUN useradd -ms /bin/bash devuser
RUN groupadd -g 999 docker || true \
    && usermod -aG docker devuser

USER devuser
WORKDIR /workspace

