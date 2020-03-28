---
title: 'Enhanced software and platform for the Town Energy Balance (TEB) model'
tags:
  - Meteorology
  - Fortran
  - Energy
  - Urban
  - Energy Balance
  - Python
  - HVAC

authors:
  - name: D. Meyer
    orcid: 0000-0002-7071-7547
    affiliation: "1, 2"
  - name: R. Schoetter
    orcid: 0000-0002-2284-4592
    affiliation: 3
  - name: V. Masson
    orcid: 0000-0001-8807-0545
    affiliation: 3
  - name: S. Grimmond
    orcid: 0000-0002-3166-9415
    affiliation: 1
affiliations:
 - name: Department of Meteorology, University of Reading, Reading, UK
   index: 1
 - name: Department of Civil and Environmental Engineering, Imperial College London, London, UK
   index: 2
 - name: CNRM UMR 3589, Université Fédérale de Toulouse, Météo-France/CNRS, Toulouse, France
   index: 3
date: 04 December 2019
bibliography: paper.bib
---

# Summary

The Town Energy Balance (TEB) model [@masson_2000_Boundary-LayerMeteorology] is a physically based single layer Urban Canopy Model (UCM) to calculate the urban surface energy balance at neighborhood scale assuming a simplified canyon geometry. It includes several capabilities (Table 1) that have been extensively evaluated offline with flux observations [@leroyer_2010_J.Appl.Meteor.Climatol; @pigeon_2008_MeteorolAtmosPhys; @lemonsu_2004_J.Appl.Meteorol; @masson_2002_J.Appl.Meteorol] and online coupled to atmospheric models such as ALARO [@gerard_2009_Mon.Wea.Rev] in ALARO-TEB [@hamdi_2012_Wea.Forecasting], the Global Environmental Multiscale (GEM; Côté et al. [-@cote_1998_Mon.WeatherRev]) in GEM-TEB [@lemonsu_2009_Boundary-LayerMeteorol], Meso-NH [@lafore_1998_Ann.Geophys; @lac_2018_Geosci.ModelDev] in TEB-MesoNH [@lemonsu_2002_Boundary-LayerMeteorology], the Regional Atmospheric Modeling System (RAMS; Pielke et al. [-@pielke_1992_Meteorl.Atmos.Phys]) in RAMS-TEB [@freitas_2007_Boundary-LayerMeteorol], the Advanced Regional Prediction System (ARPS; Xue et al., [-@Xue2000]) in ARPS-TEB [@rozoff_2003_J.Appl.Meteorol], and the Weather Research and Forecasting  (WRF; Skamarock et al. [-@skamarock_2019_NCARTech.NoteNCARTN-556STR]) in WRF-TEB [@meyer_n.a._Submitt.J.Adv.Model.EarthSyst].


Here, we present an enhanced software and platform for the TEB model to help scientists and practitioners wishing to use the TEB model in their research as a standalone software application or as a library in their own software. This includes several features such as cross-platform support for Windows, Linux, and macOS using CMake[@kitware_inc._cmake_2019], static and dynamic library generation for integration with other software/models, namelist-based configuration, integration with MinimalDX [@meyer_d_2019_3562311] and PsychroLib [@meyer_2019_JOSS] to improve the modelling of air conditioners (AC) and psychrometric calculations respectively, a thin interface used in the coupling with WRF-CMake [@riechert_2019_JOSS], helper functions for Python for pre- and post-processing inputs and outputs files, and a tutorial in Jupyter Notebook to allow users to quickly become familiar with the general TEB modeling workflow. In the new platform we implement testing at every code commit through continuous integration (CI) and automate the generation of API documentation. The project is developed as a free and open source and community-driven project on GitHub ([https://github.com/teb-model/teb](https://github.com/teb-model/teb)) and represents a significant improvement that supports existing and new model applications with enhanced functionality. We welcome contributions and encourage users to provide feedback, bug reports and feature requests, via GitHub's issue system at [https://github.com/teb-model/teb/issues](https://github.com/teb-model/teb/issues).


| Modeling capability                                                   | Reference                                                                 |
| --------------------------------------------------------------------- | ------------------------------------------------------------------------- |
| Urban Surface Energy Balance and Snow                                 | Masson [-@masson_2000_Boundary-LayerMeteorology]                          |
| Building Energy Model (BEM)                                           | Bueno et al. [-@bueno_2012_Geosci.ModelDev], Pigeon et al. [-@Pigeon2014] |
| In-canyon urban vegetation and variable road orientation              | Lemonsu et al. [-@lemonsu_2012_Geosci.ModelDev]                           |
| Green roofs, irrigation of green roofs, gardens and watering of roads | de Munck et al. [-@demunck_2013_Int.J.Climatol]                           |
| Solar panels for hot water and/or photo-voltaic (PV)                  | Masson et al.  [-@masson_2014_Front.Environ.Sci]                          |
| Human behavior related to building energy consumption*                | Schoetter et al. [-@schoetter_2017_Geosci.ModelDev]                       |
| Calculation of urban carbon dioxide fluxes*                           | Goret et al. [-@goret_2019_AtmosphericEnvironment:X]                      |
| Urban trees*                                                          | Redon et al. [-@Redon2017; -@Redon2020]                                   |

Table: Main capabilities available in the Town Energy Balance (TEB) model [@masson_2000_Boundary-LayerMeteorology]. The number of features available in TEB has increased since its first version published in 2000. (*) Capability not currently available in the TEB software.



# Acknowledgments

We acknowledge the researchers who contributed to the scientific development of the TEB code: from CNRM: Aude Lemonsu, Grégoire Pigeon, Cécile de Munck, Bruno Bueno, Marine Goret, Emilie Redon; from IFSTTAR Katia Chancibaul, Xenia Stavropulos-Laffaille; and from Environnement and Changement Climatique Canada (ECCC): Sylvie Leroyer.


# References
