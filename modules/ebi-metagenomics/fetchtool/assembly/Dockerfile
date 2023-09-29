FROM python:3.9-slim

LABEL maintainer="Microbiome Informatics"
LABEL version="0.9.0"
LABEL description="EBI Fetch Tool Docker Image."

# We need curl to download aspera and ps for nextflow monitoring
ENV DEBIAN_FRONTEND=noninteractive

RUN apt update && apt install -y curl procps && rm -rf /var/lib/apt/lists/*

COPY . .

RUN pip install --no-cache-dir .

# Aspera is an IBM library for data sharing
RUN ./install-aspera.sh

RUN export PATH=$PATH:/aspera-cli/cli/bin

CMD [ fetch-read-tool ]
