---
title: 'CloudnetTools.jl: A Julia package to apply Cloudnet to data from the ARM facility'
tags:
  - Julia
  - cloud radar
  - lidar
  - cloud physics
  - ARM
authors:
  - name: Pablo Saavedra&nbsp;Garfias
    orcid: 0000-0002-4596-946X
    corresponding: true
    affiliation: 1
  - name: Heike Kalesse-Los
    orcid: 0000-0001-6699-7040
	corresponding: false
    affiliation: 1
affiliations:
 - name: Leipzig Institute for Meteorology, Faculty of Physics and Earth Science Systems. Leipzig University, Germany
   index: 1
date: 24 May 2024
bibliography: paper.bib
repo: github.com/pablosaa/CloudnetTools.jl
---

# Summary
Atmospheric clouds play an important role in Earth's climate system. The complexity and variability of cloud formation processess make the accurate representation of their development in weather and climate models a taunting task, clouds are one of the main sources of uncertainties in our current understanding of the atmosphere.

Remote sensing of clouds is the main tool atmospheric scientists have to study clouds in a wide range of temporal and spatial resolutions from space as well as ground-based stations. For the latter, passive and active remote sensing instruments such as lidars, radars, and microwave radiometers are the main tool for cloud monitoring. Those have the capability to resolve the atmospheric observations with high vertical and temporal resolution. In order to physically characterize the clouds it is necessary to utilize classification and retrieval algorithms to obtain the cloud thermodynamical phase states as well as the specific water content within each phase. One important processing-chain was developed by [@illingworth-cloudnetcontinuous-2007] in order to evaluate operational numerical weather prediction models with cloud profiles from ground-based observations. The `Cloudnet` classification algorithm was originally implemented as proprietary software coded in MATLAB.

The Aerosol, Clouds and Trace Gases Research Infrastructure (ACTRIS) [@ACTRIS-handbook] has developed an open source version of the original propietary Cloudnet classification algorithm, named `CloudnetPy` [@tukiainen-joss-2020] that covers the full processing chain for the standard remote sensing instrumentation of the ACTRIS network. Moreover, the ACTRIS repository provides numerical weather prediction (NWP) model data necessary to run `CloudnetPy`; this, however, only for locations of observatories that are ACTRIS members. 

The Atmospheric Radiation Measurement (ARM) program of the U.S. Department of Energy is one of the worldwide most important programs carrying valuable long-term remote sensing observations for the study of the Earth's atmosphere and climate [https://arm.gov/data](arm.gov/data). The ARM program operates three mobile facilities (AMF), and three long-term permanent facilities: the South Great Plains (SPG) in Oklahoma U.S.A, the North Slope of Alaska (NSA) in Utqiaǵvik Alaska, and the Eastern North Atlantic (ENA)located on Graciosa Island, Azores, Portugal. Some of the AMR sites have been monitoring the atmosphere for more than 30 years. The ARM mobile facilities have been deployed to locations from the Arctic to the tropics, on ground- and ship-based platforms, and for campaigns typically lasting about a year. Due to the quality of long-term observations and the continued plans for future campaigns, the ARM program is one of the most important sources for the study of the atmosphere in general and the physics of clouds in particular. 

# Statement of need
The ARM program is not an ACTRIS member, which implies that none of its facilities support the ACTRIS implementation of `Cloudnet`. Moreover ACTRIS has neither actual nor future plans to support any of the ARM instruments, therefore the vast openly available dataset from the ARM facilities currently cannot be utilized with `Cloudnet`. Nevertheless, ARM has expressed solid interest to support ACTRIS activities especially within the framework of the EarthCARE satellite observations for which all ARM facilities provide valuable opportunities for ground validation activities as explained by @ARM-earthcare.
In order to use data from the ARM facility with the `Cloudnet` algorithm an intermediate step needs to be performed. This consists on adapting the ARM products to be compliant with `CloudnetPy` input formats. The package [github.com/`CloudnetTools.jl`](https://github.com/pablosaa/CloudnetTools.jl) has been developed to accomplish that purpose. The package is written in Julia language [@Julia-2017] due to its portability, speed, and seamless interoperability with languages in which `Cloudnet` was developed e.g. Python or MATLAB. The package supports a wide range of ARM instrumentation and data product levels as listed in \autoref{TBL:instrument} which are standard instrumentation at most of ARM permanent and mobile facilities. Furthermore, the required information about the atmospheric state, usually provided by Numerical Weather Prediction models, is replaced and adapted using ARM's radiosonde data. The package workflow is visualized in \autoref{FIG:workflow} where the ARM data for radar, lidar, and microwave radiometer (MWR) is the main input. The package adapts the ARM data into netCDF files in a format required by `CloudnetPy` (\autoref{FIG:workflow} arrow 1) and calls the routines to produce the target categorization (\autoref{FIG:workflow} arrow 2). In addition, `CloudnetTools.jl` post-processes the standard output from `Cloudnet` to estimate higher-level cloud properties like cloud multi-layer detection, cloud top temperature, integrated water content for every detected cloud layer, and basic visualization (\autoref{FIG:workflow} arrow 4).

+---------------------+-------------------------+-----------------------+
| Instrument          | Products                | Level                 |
+:====================+:========================+======================:+
| Radar               | `KAZR ARSCR 1kollias`   | c0, c1                |
+---------------------+-------------------------+-----------------------+
|                     | `KAZR CRFGE`            | a1                    |
+---------------------+-------------------------+-----------------------+
|                     | `KAZR CRFCORGE`         | c0                    |
+---------------------+-------------------------+-----------------------+
|                     | `MWACR CFR`             | a1                    |
+---------------------+-------------------------+-----------------------+
| Lidar               | `CEIL10m`               | b1                    |
+---------------------+-------------------------+-----------------------+
|                     | `HSRL`                  | a1                    |
+---------------------+-------------------------+-----------------------+
| Microwave Radiometer| `LOS    `               | b1                    |
+---------------------+-------------------------+-----------------------+
|                     | `RET 1liljclou `        | c1                    |
+---------------------+-------------------------+-----------------------+
| Radiosonde          | `INTERPOLATESONDE`      | c1                    |
+---------------------+-------------------------+-----------------------+

Table: List of ARM instrumentation supported by the `CloudnetTools.jl` package.\label{TBL:instrument}

One limitation for `Cloudnet` algorithms is the inavility to detect liquid cloud beyond the lidar's total attenuation, i.e. the first liquid cloud layer. In order to overcome this limitation, [@schimmel-amt-2022] have developed a convolutional neural network algorithm (`VOODOO`) as suplement for `CloudnetPy`. `VOODOO` takes advantage of features from the radar Doppler spectra signal to identify cloud droplets at heights when `Cloudnet` is unable to resolve liquid clouds. In this regard `CloudnetTools.jl` is also able to process the spectrum data files from the ARM millimeter wavelength radars, listed in \autoref{TBL:instrument}, to produce a compatible input for the `VOODOO` algorithm (\autoref{FIG:workflow} arrow 3). The processing consists of the identification of the ARM Doppler spectrum noise level which is then used to normalize the spectrum between [0,1] to rearrange them into 4D spectrograms as required by `VOODOO`, which is then stored as `zarr` files compliant to the `VOODOO` algorithm [@schimmel-amt-2022]. In summary the Julia package `CloudnetTools.jl` provides the whole processing chain for `CloudnetPy` plus `Voodoo` for the ARM dataset.

![Package workflow: with ARM data as input, `CloudnetTools.jl` generates compliant inputs for `CloudnetPy` (1), then interfaces with `CloudnetPy` to use the generated inputs (2). Alternatively it generates inputs for `VOODOO` (3) and performs post-processing upon the `CloudnetPy` retrievals (4). The shape of the blocks represent: cloud storage (cylinder), processing (rectangles), and results storage (rounded block). \label{FIG:workflow}](CloudnetTools_jl.drawio.pdf){width=60%}

# Mentions
`CloudnetTools.jl` has been primary developed and applied to Arctic cloud observations during the MOSAiC expedition [@garfias-acp-2023; @garfias2023-datamosaic], to long-term analysis of cloud properties from observations at the ARM's North Slope of Alaska (NSA) site in Utqiaǵvik, Alaska [@garfias-acp-2023], and to characterize cloud properties from observations at the East River Watershed site near Crested Butte, Colorado, during the SAIL campaign [sail.lbl.gov](https://sail.lbl.gov).

# Acknowledgements
This work was supported by the _Deutsche Forschungsgemeinschaft_ (DFG) funded Transregio-project TR-172 "Arctic Amplification (AC)3" (grant 268020496), sub-project B07 (grant 437153667) Principal Investigator J.Prof. Heike Kalesse-Los. 
The Author would like to thank to ACTRIS and FMI for making the `CloudnetPy` code open source, and the atmospheric radiation measurement (ARM) user facility, for openly providing the valuable data. 

# References
