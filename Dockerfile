# depthSV — analysis pipeline image
#
#   docker build -t depthsv:dev .
#   docker run --rm -v "$PWD:/work" depthsv:dev tests/smoke_test.sh /work/scratch
#
# rocker/r-ver pins both the R version and a dated CRAN snapshot, so a rebuild
# resolves the same package versions rather than whatever is current. That is
# the property that makes a container worth having here.

FROM rocker/r-ver:4.4.1

# bgzip and tabix come from the tabix package; procps provides the `ps` that
# GNU parallel uses to size its job pool.
RUN apt-get update && apt-get install -y --no-install-recommends \
        tabix \
        parallel \
        procps \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# survival ships with R itself; only these two are extra.
RUN install2.r --error --skipinstalled optparse data.table

# GNU parallel prints a citation notice to stderr on first run in each new HOME.
# Silencing it keeps stage logs to what the pipeline actually reports.
RUN mkdir -p /root/.parallel && touch /root/.parallel/will-cite
ENV PARALLEL_HOME=/root/.parallel

COPY lib/    /opt/depthsv/lib/
COPY R/      /opt/depthsv/R/
COPY scripts/ /opt/depthsv/scripts/
COPY tests/  /opt/depthsv/tests/
COPY conf/   /opt/depthsv/conf/

RUN chmod +x /opt/depthsv/scripts/*.sh /opt/depthsv/tests/*.sh

ENV PATH=/opt/depthsv/scripts:${PATH}
WORKDIR /work

# No ENTRYPOINT: the stages are separate commands and a workflow engine calls
# them individually. `docker run ... scripts/correct.sh --region chr1 ...`
CMD ["bash"]
