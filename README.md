# Coxmos: Cox MultiBlock Survival Analysis for high dimensional and multiomic data

**Coxmos** is still a beta-version. Work in progress. We strongly recommend to not use it yet.

* [Introduction](https://github.com/ConesaLab/Coxmos/edit/master/README.md#introduction)

* [Installation](https://github.com/ConesaLab/Coxmos/edit/master/README.md#installation)

* [Getting started](https://github.com/ConesaLab/Coxmos/edit/master/README.md#getting-started)

* [Contact](https://github.com/ConesaLab/Coxmos/edit/master/README.md#contact)

* [References](https://github.com/ConesaLab/Coxmos/edit/master/README.md#references)


### Introduction
The **coxmos** R package is an end-to-end pipeline designed for the study of survival analysis for 
high dimensional data. Updating classical methods and adding new ones based on sPLS technologies. 
Furthermore, includes multiblock functions to work with multiple sets of information to improve 
survival accuracy. 

The pipeline includes three basic analysis blocks:

1. Computing cross-validation functions and getting the models. 

2. Evaluating all the models to select the better one for multiple metrics.

3. Understanding the results in terms of the global model and the original variables.

*Coxmos* contains the necessary functions and documentation to obtain from raw data the final models
after compare them, evaluate with test data, study the performance individually and in terms of 
components and graph all the results to understand which variables are more relevant for each case 
of study.

![](images/logo.png)

### Installation

#### Dependencies requiring manual installation

Some of the metrics available in *coxmos* are optional based and will not be included in the 
standard *coxmos* installation. A list of all optional packages are shown below:

* nsROC:
* smoothROCtime:
* survivalROC:
* risksetROC:
* ggforce:
* RColorConesa:

#### Installing coxmos

The *coxmos* R package and all the remaining dependencies can be installed 
from GitHub using `devtools`:

```
devtools::install_github("ConesaLab/Coxmos")
```

To access vignettes, you will need to force building with
`devtools::install_github(build_vignettes = TRUE)`. Please note that this will
also install all suggested packages required for vignette build and might 
increase install time. Alternatively, an HTML version of the vignette is
available under the [vignettes](https://github.com/ConesaLab/Coxmos/tree/master/vignettes)
folder.


### Getting started

In order to use coxmos, you will need the following items:

- A explanatory X matrix.
- A response survival Y matrix (with two columns, "time" and "event").

Please note that two toy datasets are included in the package. Details to load and use them can be 
found in the package's vignette.


### Contact
If you encounter a problem, please 
[open an issue](https://github.com/ConesaLab/Coxmos/issues) via GitHub.

  
### References
If you use *coxmos* in your research, please cite the original publication:

