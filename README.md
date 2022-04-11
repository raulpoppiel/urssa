# `urssa` Unsupervised Routine Soil Spectral Analysis  <img align="right" width="160px" src="https://user-images.githubusercontent.com/28542789/162770776-30b1032c-248e-4d2b-b221-af4ac92c1856.png"/>

<p align='left'>
  <a href="#"><img src="https://img.shields.io/badge/repo%20status-100%25%20ready-green"></a>
    <a href="#"><img src="https://img.shields.io/badge/Last%20update-05--Mar--2022-blueviolet"></a>
  <a href="#"><img src="https://img.shields.io/github/downloads/raulpoppiel/urssa/ursa-v1.0.0/total?color=green&label=downloads&style=flat"></a>
  
  <a href="#"><img align='left' src="https://badges.pufler.dev/visits/raulpoppiel/urssa"></a> pedometricians have come here

\
<em><p align="right"> Bridging the gap between soil spectroscopy and traditional laboratory </p></em>

## About
The `urssa` code provides basic unsupervised functionalities for the use of spectra in laboratory routines. 
This repository provides a csv file conainting 350-2500 nm soil spectra and codes (divided into modules) following a reproducible example.


## Core functionalities (routines)
- Unsupervised spectral clustering of samples
- Internal quality control (outlier detection and removal) of analytical results
- Correlation analysis between soil data and spectra
- Samples selection for traditional analysis
- Samples selection for soil prediction (cost reduction)

## Data available for download
The `Soil_data_spectra.csv` file contains:
* Soil Attributes data: 
    - Clay (g kg-1)
    - Sand (g kg-1) 
    - Organic Matter - OM (g kg-1) 
    - Cation Exchange Capacity - CEC (mmolc kg-1) 
    - Base Saturation - V (%)
* Soil reflectance spectra: 
    - From 350 to 2500 nm at 1 nm resolution

## R Codes
`urssa_01.R`:
- MODULE 1: Importing and pre-processing data
- MODULE 2: Unsupervised spectral clustering of samples
- MODULE 3: Variable importance
- MODULE 4: Identification of soil attribute outliers
- MODULE 5: Correlation analysis
- MODULE 6: Plotting boxplot and spectra by cluster and Laboratory

`urssa_02.R`:
- MODULE 1: Quantification of outliers

`urssa_03.R`:
- MODULE 1: Splitting data into training and test subsets
- MODULE 2: Assessment of subsets
- MODULE 3: Soil attributes modeling with CUBIST
- MODULE 4: Assessment of prediction results

## Methodological flowchart
<div align="center">
<img src="https://user-images.githubusercontent.com/28542789/159790381-c3e5c02c-ba65-4dec-94ec-c648f62e782f.png" width="700px" />
</div>

## Reference
Please, cite the following paper when using `urssa`:

> Poppiel, R.R.; Paiva, A.F.S.; Demattê, J.A.M. Bridging the gap between soil spectroscopy and tradition-al laboratory: insights for routine implementation. Geoderma, 2022. DOI: [https://](https://)
