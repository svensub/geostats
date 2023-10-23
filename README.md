## Geostatistical Modelling

The repository contains code used for a project that employed geospatial methods to spatially predict treatment-priority areas for onchocerciasis (commonly known as river-blindness) in Tanzania. The data used for this project was obtained from Rapid Epidemiological Mapping of Onchocerciasis (REMO) data for Tanzania; the data contains patient information and cannot be made publically available but this code can be used with other spatial datasets. The primary package used in this project was the PrevMap R package (https://doi.org/10.18637/jss.v078.i08), which implements the fitting of geostatistical models for binomial data. The commented code used for the analysis can be found under the main_code.R document.

## Initial exploratory analysis

Under the section "Handling Spatial Data", an initial exploratory analysis is carried out, comprising of 1) a linear spline regression analysis of log-transformed precipitation on logit-transformed prevalence; 2) an assessment of residual spatial correlation. This section of the code first plots a point-map to demonstrate empirical prevalence of of river-blindness by sampled location. Following this, an empirical logit transformation is applied to explore the relationship between precipitation (the chosen spatially referenced candidate covariate) and onchocerciasis prevalence, and a linear splines are fitted. The code further estimates an empirical variogram to explore the presence of spatial correlation.

## Model formulation

Under the section "Analysis", an empirical variogram is estimated to explore the presence of spatial correlation. The code fits both a linear geostatistical model and a binomial geostatistical model of prevalence. A Monte Carlo maximum likelihood method is employed to estimate parameters. The code also carries out spatial prediction for prevalence and exceedance probabilities, as well as variogram-based validation.

