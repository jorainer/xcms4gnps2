FROM bioconductor/bioconductor_docker:RELEASE_3_23

LABEL name="jorainer/xcms4gnps2" \
      url="https://github.com/jorainer/xcms4gnps2" \
      maintainer="johannes.rainer@gmail.com" \
      description="Docker container to run xcms-based preprocessing for GNPS2. This version bases on the Bioconductor release 3.22 docker image." \
      license="Artistic-2.0"

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

## Global installation of required packages
RUN Rscript -e "BiocManager::install(c('xcms', 'MsExperiment', 'mzR', 'remotes', 'pak') , ask = FALSE, dependencies = c('Depends', 'Imports'), build_vignettes = FALSE)"

## Install MsBackendMassIVE from GitHub; change for RELEASE_3_24
RUN Rscript -e "BiocManager::install('RforMassSpectrometry/MsBackendMassIVE', ref = 'gabri')"

## Cache the MS data set
RUN Rscript -e "library(Spectra); library(MsBackendMassIVE); Spectra('MSV000090156', filePattern = 'Lab_2/Interlab-LC-MS_Lab2.*mzML$', source = MsBackendMassIVE())"

## Move the cache to the rstudio user and cross-link
RUN mkdir -p /home/rstudio/.cache/R/BiocFileCache && \
    mv /root/.cache/R/BiocFileCache /home/rstudio/.cache/R/ && \
    chown -R rstudio:rstudio /home/rstudio/.cache && \
    mkdir -p /root/.cache/R/ && \
    ln -s /home/rstudio/.cache/R/BiocFileCache /root/.cache/R/BiocFileCache

## Install package dependencies
RUN Rscript -e "pak::local_install_dev_deps(ask = FALSE)"

## Install the current package with vignettes
RUN Rscript -e "pak::local_install()"

## Clean up
RUN find vignettes/ -name "*.html" -type f -delete && find vignettes/ -name "*_files" -type d -exec rm -r {} + && \
    rm -rf /tmp/* && \
    chown -R rstudio:rstudio /home/rstudio/.cache
