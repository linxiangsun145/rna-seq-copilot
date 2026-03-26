FROM python:3.11-slim

# Install R and required packages
RUN apt-get update && apt-get install -y \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages (cached layer)
RUN Rscript -e "install.packages(c('BiocManager', 'optparse', 'jsonlite', 'ggplot2', 'ggrepel', 'pheatmap', 'RColorBrewer'), repos='https://cloud.r-project.org')"
RUN Rscript -e "BiocManager::install(c('DESeq2'), ask=FALSE)"

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy backend source
COPY . /app/

# Ensure jobs dir exists at runtime
RUN mkdir -p /app/jobs

EXPOSE 8000

# Run from /app (backend dir) so non-prefixed imports resolve correctly
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
