# Minimal, reproducible image built from the unified envs/tracktx.yaml
FROM mambaorg/micromamba:1.5.10-jammy

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG TARGETPLATFORM
ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_NO_BANNER=1 \
    PIP_NO_CACHE_DIR=1

SHELL ["/bin/bash", "-lc"]

# Install Conda env
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/tracktx.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && micromamba clean -a -y

# (Optional) small OS tools that don’t exist in Conda packages
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates procps tini \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Runtime niceties
ENV PATH=/opt/conda/bin:$PATH \
    LC_ALL=C.UTF-8 LANG=C.UTF-8
WORKDIR /data
LABEL org.opencontainers.image.source="https://github.com/serhataktay/tracktx" \
      org.opencontainers.image.title="TrackTx" \
      org.opencontainers.image.description="Transcription tracking pipeline image (Nextflow)" \
      org.opencontainers.image.licenses="Apache-2.0"

ENTRYPOINT ["/usr/bin/tini","--"]
CMD ["bash"]