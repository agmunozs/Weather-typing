# Flow-dependent Cross-timescale Model Diagnostics 
Several tools for flow-dependent model diagnostics and general weather typing, as part of our contribution to NOAA's Model Diagnostics Task Force (MDTF). More information:
A Weather-Type-Based Cross-Time-Scale Diagnostic Framework for Coupled Circulation Models 
Ángel G. Muñoz;  Xiaosong Yang;  Gabriel A. Vecchi;  Andrew W. Robertson;  William F. Cooke
J. Climate (2017) 30 (22): 8951–8972.
https://doi.org/10.1175/JCLI-D-17-0115.1

# Weather-typing
Several weather-typing codes, in Matlab and Python.

## PyWR
Python scripts and Jupyter notebooks to compute weather types/regimes diagnostics using K-means. Model weather types are projected (or not) into the observed ones in the EOF space.

Procrustes decomposition and callibration of computed weather types. Plotting functionality to view weather types comparison (reanalysis vs model data), weather types decomposition and callibration data for each computed weather type.

### Authors:
Ángel G. Muñoz (agmunoz@iri.columbia.edu), Drew M. Resnick (drewr@iri.columbia.edu)
### Collaborators:
Andrew W. Robertson (awr@iri.columbia.edu), James Doss-Gollin (james.doss-gollin@columbia.edu)

## Matlab scripts
Matlab scripts to compute weather types/regimes using K-means. Model weather types are not projected into observed ones.
### Authors:
Ángel G. Muñoz (agmunoz@iri.columbia.edu) modifications to the original code by Andrew W. Robertson (awr@iri.columbia.edu)

## Function documentation
In order to view the html function documentation, do the following in your terminal:
1. Make sure you have [sphinx](https://anaconda.org/conda-forge/sphinx) installed.
2. Navigate to the `Weather-typing/docs/` subdirectory within the repository.
3. Edit line 15 of `source/conf.py` to be the path of the top level of the repository on your computer.
4. In `docs/`, run the command:
  `make html`
5. The HTML pages are in `Weather-typing/docs/build/html/`; Navigate to the subdirectory in your computer's finder window and 
   double click on `index.html` to open in browser.

Once done, you can access the documentation any time by opening `index.html`. However, any time you pull updates from the main repository, it is recommended you remake the documentation in case it was updated. To do this, navigate to `Weather-typing/docs/` and:
1. Enter `make clean` to remove the outdated `index.html` from the `build/` folder.
2. Enter `make html` to rebuild the updated `index.html` in the `build/` folder.
