## Declares build arguments
ARG NB_USER=rstudio
ARG NB_UID=1001

FROM rocker/binder:4.4.1

ENV DEBIAN_FRONTEND=noninteractive
USER root

## Install linux packages
RUN apt-get update --fix-missing > /dev/null \
        && apt-get install --yes \
          libncurses5-dev \
        && apt-get clean > /dev/null \
        && rm -rf /var/lib/apt/lists/*

## Install SAMTOOLS
RUN curl -kL https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 \
   | tar -C /tmp/ -jxf - && cd /tmp/samtools-1.19.2/ \
   && ./configure && make -j4 && make install \
   && rm -rf /tmp/samtools-1.19.2/

## Install BCFTOOLS
RUN curl -kL https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 \
   | tar -C /tmp/ -jxf - && cd /tmp/bcftools-1.19/ \
   && ./configure && make -j4 && make install \
   && rm -rf /tmp/bcftools-1.19/

## Install MINIMAP2
RUN curl -kL https://github.com/lh3/minimap2/releases/download/v2.27/minimap2-2.27.tar.bz2 \
   | tar -C /tmp/ -jxf - \
   && cd /tmp/minimap2-2.27/ \
   && (if [ "$TARGETPLATFORM" = "linux/arm64" ]; then \
      make arm_neon=1 aarch64=1; \
  else \
      make; \
  fi) \
  && cp ./minimap2 ./misc/paftools.js /usr/local/bin/ \
  && rm -rf /tmp/minimap2-2.27/


USER ${NB_USER}

## Install R packages
RUN Rscript -e '\
  install.packages(c("tidyverse","rmarkdown","kableExtra","bs4Dash","BiocManager"),Ncpus=4L);\
  BiocManager::install(c("Biostrings","GenomicAlignments","Rsamtools","GenomicRanges"),Ncpus=4L)'


## Copy files to home folder
COPY --chown=${NB_USER} ./ ${HOME}/
