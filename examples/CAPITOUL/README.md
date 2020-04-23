# CAPITOUL case

In this folder data are provided to run TEB (Masson, 2000) with sample input data from the CAPITOUL (Canopy and Aerosol Particles Interactions in TOulouse Urban Layer; Masson et al., 2008) campaign. Data are provided for demonstration purposes only and are not quality controlled (i.e. should not be used to conduct scientific experiments). Both data and tutorial are only meant to guide the user through setting up a simple experiment in TEB. If you are looking to conduct a scientific experiment, please download and process the data yourself (see the [software documentation page](../../docs/software-docs.md)) and amend the input namelist accordingly (see the [namelist options page](../../docs/namelist-options)). 

## Install Python packages

Several Python packages are required to run the tutorial. To install these, execute the following from the command line:

```
pip install -r requirements.txt
```

## Tutorial

To run the tutorial, execute the following from the command line:

```
jupyter notebook tutorial.ipynb
```

Alteratively, your can view its content in [`tutorial.ipynb`](tutorial.ipynb).


## References

> Masson, V., 2000: A Physically-Based Scheme For The Urban Energy Budget In Atmospheric Models. Boundary-Layer Meteorology, 94, 357–397, https://doi.org/10.1023/A:1002463829265.

> Masson, V., and Coauthors, 2008: The Canopy and Aerosol Particles Interactions in TOulouse Urban Layer (CAPITOUL) experiment. Meteorol Atmos Phys, 102, 135–157, https://doi.org/10.1007/s00703-008-0289-4.