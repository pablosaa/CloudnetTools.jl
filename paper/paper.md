---
title: 'CloudnetTools.jl: A Julia package for running Cloudnet with ARM data'
tags:
  - Julia
  - cloud radar
  - lidar
  - microwave radiometer
  - remote sensing
  - Arctic
authors:
  - name: Pablo Saavedra~Garfias
    orcid: 0000-0002-4596-946X
    corresponding: true
    affiliation: 1
  - name: Heike Kalesse-Los
    orcid: 0000-0001-6699-7040
    affiliation: 1
affiliations:
 - name: Leipzig Institut für Meteorologie, Leipzig Universität, Leipzig, Germany
   index: 1
date: 23 December 2023
bibliography: paper.bib
---

# Summary
The Aerosol, Clouds and Trace Gases Research Infrastructure (ACTRIS) [@ACTRIS-handbook] has developed an open source version of the original propietary Cloudnet classification algorithm by [@illingworth-cloudnetcontinuous-2007] that covers the full processing chain for the ACTRIS network [@tukiainen-joss-2020]. 

# Statement of need
The Atmospheric Radiation Measurement (ARM) program of the Department of Energy in the United States of America is one of the most important program carring on long-term and valuable remote sensing observations for the study of the atmosphere and climate (https://arm.gov/data). The ARM program operates several observational facilities, some of them have been monitoring the atmosphere for more than 30 years. The ARM facilites comprises of permanent locations as well as deployments by its mobile facilitis (MF) to an unprecedented number of locations ranging from the Amazon to the Arctic.

Although the Cloudnetpy algorithm supports a rich number of remote sensing instrumentation of the ACTRIS member organizations, none of the ARM facility instrumentation is suported. Therefore, in order to use cloudnetpy with data from the ARM facility an intermediate step needs to be performed. This consists on the adaptation of the ARM data files to cloudnetpy supported inputs. Furthermore, ACTRIS provides the necessary numerical weather prediction model data for all observatories which are ACTRIS members. Therefore, for the usage with ARM data this needs to be replaced with radiosonde data which is similarly provided by ARM as an standard product for all its sites.

+---------------------+-------------------------+-----------------------+
| Instruement         | Products                | Level                 |
+:====================+:=======================:+:=====================:+
| RADAR               | `KAZR ARSCR 1kollias`   | c0, c1                |
+---------------------+-------------------------+-----------------------+
|                     | `KAZR CORGE`            | a1                    |
+---------------------+-------------------------+-----------------------+
|                     |  `KAZR CRF GE`          | a1                    |
+---------------------+-------------------------+-----------------------+
|                     |  `MWACR CFR`            | a1                    |
+---------------------+-------------------------+-----------------------+
| LIDAR               | `CEIL10m`               | b1                    |
+---------------------+-------------------------+-----------------------+
|                     | `HSRL`                  | a1                    |
+---------------------+-------------------------+-----------------------+
| MWR                 | `LOS    `               | b1                    |
+---------------------+-------------------------+-----------------------+
|                     | `RET 1liljclou `        | c1                    |
+---------------------+-------------------------+-----------------------+
| RADIOSONDE          | `INTERPOLATESONDE`      | c1                    |
+---------------------+-------------------------+-----------------------+
| MODEL               | `ECMWF`      `          | ~~not ARM ~~          |
+---------------------+-------------------------+-----------------------+

# Mentions
CloudnetTools.jl has been primary developed and applied to Arctic cloud observations during the MOSAiC expedition [garfias-acp-2023, garfias2023-datamosaic] as well as the long-term observations provided by the ARM's North Slope of Alaska (NSA) site in Utqiagvik, Alaska [garfias-acp-2023]. Furthermore, CloudnetTools.jl has been used to classify and retrieve cloud properties from observations at a montaneous region during the SAIL campaign (https://meteo.uni-leipzig.de).

# Acknowledgements
This work was supported by the _Deutsche Forschungsgemeinschaft_ (DFG) funded Transregio-project TR-172 "Arctic Amplification (AC)3" (grant 268020496), sub-project B07 (grand ). 
The Author would like to thank to ACTRIS for making the Cloudnetpy code open source, and the atmospheric radiation measurement (ARM) user facility, for providing the valuable data. 

# References
