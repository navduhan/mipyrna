ARG MIPYRNA_PLATFORM=linux/amd64
FROM --platform=${MIPYRNA_PLATFORM} condaforge/miniforge3:latest

LABEL maintainer="Naveen Duhan <naveen.duhan@usu.edu>" \
      description="miPyRNA: end-to-end small RNA-seq analysis toolkit" \
      version="1.0.0" \
      license="MIT"

RUN apt-get update && \
    apt-get install -y --no-install-recommends bzip2 ca-certificates curl procps wget && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /mipyrna
COPY mipyrna_environment.yaml /mipyrna/mipyrna_environment.yaml

RUN mamba env create -y -f /mipyrna/mipyrna_environment.yaml && \
    mamba clean -afy

COPY . /mipyrna
RUN /opt/conda/envs/mipyrna/bin/pip install --no-deps . && \
    echo "conda activate mipyrna" >> /root/.bashrc

ENV PATH="/opt/conda/envs/mipyrna/bin:${PATH}"
WORKDIR /data

CMD ["mipyrna", "--help"]
