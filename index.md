# *xcms*-based LC-MS/MS data preprocessing for FBMN with GNPS2

[![License: CC BY-NC
4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)\
[![Docker Image Version (latest by
date)](https://img.shields.io/docker/v/jorainer/xcms4gnps2?label=docker%20image)](https://hub.docker.com/repository/docker/jorainer/xcms4gnps2)\
[![DOI](https://zenodo.org/badge/1099034068.svg)](https://doi.org/10.5281/zenodo.17800013)

This repository contains an example workflow for preprocessing and
preparation of an LC-MS/MS data set for feature-based molecular
networking (FBMN) with GNPS2.

## Analysis workflow

### [xcms-based preprocessing for FBMS with GNPS](https://jorainer.github.io/xcms4gnps2/articles/MSV000090156-preprocessing.html)

This workflow explains preprocessing of a public LC-MS/MS data set with
*xcms* and export of the results for feature-based molecular networking
with GNPS.

🎥 a video recording of a presentation of this workshop is available at
<https://youtu.be/yc6fsegFg-k>.

------------------------------------------------------------------------

## 📌 Reproducibility & Usage

The workflow is available as pre-rendered webpage
[xcms4gnps](https://jorainer.github.io/xcms4gnps2). In addition, a
[docker](https://docker.com) image is available allowing to run the
workflow interactively:

- If you don’t already have, install [docker](https://www.docker.com/).
  Find installation information
  [here](https://docs.docker.com/desktop/).
- Get the [docker image](https://hub.docker.com/r/jorainer/xcms4gnps2)
  of this tutorial e.g. from the command line with:

&nbsp;

    docker pull jorainer/xcms4gnps2:latest

- Start the docker container, either through the Docker Desktop, or on
  the command line with

&nbsp;

    docker run -e PASSWORD=bioc -p 8787:8787 jorainer/xcms4gnps2:latest

- Enter [`http://localhost:8787`](http://localhost:8787) in a web
  browser and log in with username `rstudio` and password `bioc`.
- In the RStudio server version: open any of the Quarto files in the
  *vignettes* folder and evaluate the R code blocks in that document.

ℹ️ macOS users might need to emulate an AMD64 CPU to run the docker
image:

    docker pull --platform linux/amd64 jorainer/xcms4gnps2:latest
    docker run --platform linux/amd64 -e PASSWORD=bioc -p 8787:8787 jorainer/xcms4gnps2:latest

As an alternative to the docker-based workshop, it is also be possible
to run and evaluate the workshop *natively* in a local R installation (R
version \>= 4.6):

``` r

#' install the workshop and all required dependencies
install.packages(c("BiocManager", "remotes"))
BiocManager::install("jorainer/xcms4gnps2", dependencies = TRUE)
BiocManager::install("RforMassSpectrometry/MsBackendMassIVE")
```

The
[*vignettes/MSV000090156-preprocessing.qmd*](https://jorainer.github.io/xcms4gnps2/vignettes/MSV000090156-preprocessing.qmd)
file could then be opened with any editor (after downloading locally)
and the individual code lines evaluated in an R session.

## 🤝 Contribution

Interested in contributing? Please check out the [**RforMassSpectrometry
Contributions
Guide**](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#contributions).

### 📜 Code of Conduct

We follow the [**RforMassSpectrometry Code of
Conduct**](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#code-of-conduct)
to maintain an inclusive and respectful community.

------------------------------------------------------------------------

## 🙌 Acknowledgements

![EU
Logo](https://github.com/rformassspectrometry/Metabonaut/raw/main/vignettes/images/EULogo.jpg)

EU Logo

Part of this work is funded by the **European Union** under the
**HORIZON-MSCA-2021** project **101073062: HUMAN – Harmonising and
Unifying Blood Metabolic Analysis Networks**.

🔗 Learn more: [HUMAN Project Website](https://human-dn.eu/)
