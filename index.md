# *xcms*-based LC-MS/MS data preprocessing for FBMN with GNPS2

[![License: CC BY-NC
4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)  
[![Docker Image Version (latest by
date)](https://img.shields.io/docker/v/jorainer/xcms4gnps2?label=docker%20image)](https://hub.docker.com/repository/docker/jorainer/xcms4gnps2)

This repository contains an example workflow for preprocessing and
preparation of an LC-MS/MS data set for feature-based molecular
networking (FBMN) with GNPS2.

## Analysis workflow

### [xcms-based preprocessing for FBMS with GNPS](https://jorainer.github.io/xcms4gnps2/articles/MSV000090156-preprocessing.html)

This workflow explains preprocessing of a public LC-MS/MS data set with
*xcms* and export of the results for feature-based molecular networking
with GNPS.

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
