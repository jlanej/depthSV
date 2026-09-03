# depthSV — analysis pipeline image
#
#   docker build -t depthsv:dev .
#   docker run --rm -v "$PWD:/work" depthsv:dev tests/smoke_test.sh /work/scratch
#
# rocker/r-ver pins both the R version and a dated CRAN snapshot, so a rebuild
# resolves the same package versions rather than whatever is current. That is
# the property that makes a container worth having here. The Debian packages
# below come from the base image's pinned distribution, so their versions
# move only when the base image does.
#
# CI publishes this image for linux/amd64 on every release tag
# (ghcr.io/<owner>/depthsv:<version>) and prints the digest; workflow inputs
# should pin that digest.

FROM rocker/r-ver:4.4.1

ARG DEPTHSV_VERSION=dev
ARG DEPTHSV_COMMIT=unknown
LABEL org.opencontainers.image.title="depthSV" \
      org.opencontainers.image.description="Read-depth association testing at biobank scale" \
      org.opencontainers.image.source="https://github.com/jlanej/depthSV" \
      org.opencontainers.image.version="${DEPTHSV_VERSION}" \
      org.opencontainers.image.revision="${DEPTHSV_COMMIT}" \
      org.opencontainers.image.licenses="MIT"

# bgzip and tabix come from the tabix package; procps provides the `ps` that
# GNU parallel uses to size its job pool; curl fetches the GCS token when a
# workflow streams the matrix from a bucket (rocker removes it).
RUN apt-get update && apt-get install -y --no-install-recommends \
        tabix \
        parallel \
        procps \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# survival ships with R itself; only these two are extra.
RUN install2.r --error --skipinstalled optparse data.table

# GNU parallel prints a citation notice to stderr on first run in each new
# HOME. A world-readable PARALLEL_HOME silences it for whichever user the
# workflow engine runs the task as, not only root.
RUN mkdir -p /opt/parallel && touch /opt/parallel/will-cite && chmod -R a+rwX /opt/parallel
ENV PARALLEL_HOME=/opt/parallel

COPY VERSION   /opt/depthsv/VERSION
COPY lib/      /opt/depthsv/lib/
COPY R/        /opt/depthsv/R/
COPY scripts/  /opt/depthsv/scripts/
COPY tests/    /opt/depthsv/tests/
COPY conf/     /opt/depthsv/conf/
COPY workflows/ /opt/depthsv/workflows/

RUN chmod +x /opt/depthsv/scripts/*.sh /opt/depthsv/tests/*.sh \
    && printf '%s\n' "${DEPTHSV_COMMIT}" > /opt/depthsv/COMMIT

ENV PATH=/opt/depthsv/scripts:${PATH}
WORKDIR /work

# No ENTRYPOINT: the stages are separate commands and a workflow engine calls
# them individually. `docker run ... scripts/correct.sh --region chr1 ...`
CMD ["bash"]
