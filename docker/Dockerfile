FROM mambaorg/micromamba:debian13-slim

COPY python/requirements.txt python/requirements-dev.txt ./

RUN micromamba install -y -n base -c conda-forge \
      python=3.13 numpy biopython scipy \
      pip && \
      micromamba clean --all --yes
