FROM ubuntu:20.04

# Install necessary dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    software-properties-common \
    python3 \
    python3-pip \
    r-base \
    r-cran-tidyverse \
    r-cran-biostrings \
    r-cran-rsamtools \
    minimap2 \
    samtools \
    porechop \
    nanopack \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN apt-get update && apt-get install -y \
    openjdk-11-jdk \
    && rm -rf /var/lib/apt/lists/*
RUN curl -s https://get.nextflow.io | bash
ENV PATH="/root/nextflow:$PATH"

# Copy the project files
COPY . /app
WORKDIR /app

# Set the entrypoint to run the Nextflow pipeline
ENTRYPOINT ["nextflow", "run", "main.nf"]