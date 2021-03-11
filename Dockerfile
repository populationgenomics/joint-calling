# Inheriting from a service image that has conda and hail installed
FROM australia-southeast1-docker.pkg.dev/vlad-dev/test-ci/hailbatch:latest
MAINTAINER Centre for Population Genomics "https://github.com/populationgenomics"

RUN mkdir -p /work
WORKDIR /work
COPY setup.py /work/
COPY scripts /work/scripts
COPY cpg_qc /work/cpg_qc
COPY test /work/test
COPY conda /work/conda
COPY README.md /work/README.md

RUN conda install conda-build conda-verify anaconda-client
RUN conda build conda/cpg-qc
RUN conda create --use-local -n env cpg-qc
ENV PATH /miniconda/envs/env/bin:$PATH

RUN wget https://broad.io/install-gcs-connector

# Clean up
RUN rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/tmp/* && \
    cd /usr/local && \
    apt-get clean && \
    rm -rf /.cpanm
