## General

This repository contains code and data for fitting the model and running the simulation study presented in "Improving ecological inference and uncertainty quantification from camera trap data through the fusion of AI confidences and manual annotations." Please contact Adira Cohen (acohen6@ncsu.edu) with any questions.

Computation was done using R version 4.5.2.

## Files

The following R files are included in the models folder.

```OrdinalCompModelFuncs.R```
* Functions to support model fitting, including the main MCMC algorithm, likelihoods, priors, and plotting functions.
* This is the only file necessary for fitting our model to data.

```OrdinalCompModelSim.R```
* Code to run the simulation study in parallel as presented in the manuscript.
* This generates the files found in the sim_results folder.
* The framework for the sim study is that for each epoch, a theoretical "complete" dataset is generated using ```gen_sim_data``` which is paired down (e.g., removing 80% of annotations) using ```convert_sim_data``` and fit using ```MH_alg``` for each setting.
* The JAGS package is only used for setting 1, the linear model, so users can skip this installation unless they wish to run setting 1. 

```OrdinalCompModelSetup.R```
* Handles settings for the simulation study, such as what percent of the data to annotate and the number of sites.

```OrdinalCompModelSimFuncs.R```
* Code to support the simulation study, mainly generating and converting simulated data.

```OrdinalCompModelAnalysis.R```
* Code to generate tables and figures for the simulation study.

## Running the model

The file ```example.R``` was included to demonstrate a simple example of running the model. To use your own data instead of simulated data, it is important to format your data as a ```list``` with the following elements. Matrices and vectors that correspond to images, sequences, or sites must be ordered consistently.

```L```
* The number of categories in the response (e.g., 5 in the paper)

```N```
* The number of sites

```R```
* The number of images

```r```
* The number of images in each sequences

```seq_ID```
* The sequence that each image belongs to

```p```
* The number of environmental covariates

```q```
* The number of image quality covariates

```A```
* The number of annotators
* If you do not wish to differentiate between annotators, set this to 1

```Z```
* Vector of integer manual annotations
* Non-annotated images should be ```NA```

```a```
* Vector of integers representing the annotator for each image
* If you do not wish to differentiate between annotators, this will be a vector of 1s
* Non-annotated images should be ```NA```

```Z_labelled```
* Logical vector indicating whether the image is annotated (```TRUE```) or not (```FALSE```)

```C```
* Matrix of AI confidences, where each row represents the confidence vector for that image
* Rows for non-AI classified images should be ```NA```

```C_labelled```
* Logical vector indicating whether the image is AI classified (```TRUE```) or not (```FALSE```)

```C_tilde```
* Rounded version of ```C``` to remove zeros, which cause errors, using the provided function ```fix_zeros(C)```

```X```
* Matrix of site covariates, where each row are the covariates for that site
* The provided plotting functions use the column names of ```X``` to improve readability

```U```
* Matrix of image quality covariates, where each row are the covariates for that image
* The provided plotting functions use the column names of ```U``` to improve readability
