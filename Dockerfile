FROM bioconductor/bioconductor_docker:RELEASE_3_22

LABEL name="jorainer/xcms4gnps2" \
      url="https://github.com/jorainer/xcms4gnps2" \
      maintainer="johannes.rainer@gmail.com" \
      description="Docker container to run xcms-based preprocessing for GNPS2. This version bases on the Bioconductor release 3.22 docker image." \
      license="Artistic-2.0"

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

## Global installation of required packages
RUN Rscript -e "BiocManager::install(c('xcms', 'MsExperiment', 'mzR') , ask = FALSE, dependencies = c('Depends', 'Imports'), build_vignettes = FALSE)"

USER rstudio

## Download the data and store it to the local folder
RUN wget -r ftp://massive-ftp.ucsd.edu/v04/MSV000090156/peak/mzml/POS_MSMS/Lab_2/* -P /home/rstudio/data

## root user needed for rstudio server properly working
USER root

## Install the current package with vignettes
RUN Rscript -e "devtools::install('.', dependencies = c('Depends', 'Imports'), type = 'source', build_vignettes = TRUE, repos = BiocManager::repositories())"

## Clean up
RUN find vignettes/ -name "*.html" -type f -delete && find vignettes/ -name "*_files" -type d -exec rm -r {} + && \
    rm -rf /tmp/*
