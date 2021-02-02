# KenyaSerology

This repository contains the openly available data and code in order to reproduce the results of *Revealing the extent of the COVID-19 pandemic in Kenya based on serological and PCR-test data*, and/or, use the package with data from another setting.

## Prerequisites and recommended background knowledge:
* Access and basic familiarity with the [Julia programming language](https://julialang.org/).
* The underlying dynamical system representing the unobserved infection process is a modification to the basic SEIR model that is described in [this paper](https://journals.sagepub.com/doi/full/10.1177/0962280217747054), albeit we implement a continuous time rather than discrete version of the model.
* Solutions of the infection process are generated using the performant, and well documented, package [SciML/DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). Familiarisation with this package is desirable.
* Hamiltonian MCMC (HMC) is implemented using the [dynamicHMC.jl](https://github.com/tpapp/DynamicHMC.jl) package. The log-likelihood function for parameters is directly defined in KenyaSerology, and log-likelihood gradients (necessary for HMC) are calculated using forward-mode automatic differentiation. The combination of ODE solutions and log-likelihood function gradients in code was inspired by [DiffEqBayes.jl](https://github.com/SciML/DiffEqBayes.jl). A good conceptual introduction to HMC can be found [here](https://arxiv.org/abs/1701.02434).

## Data sets in this repository

`/opendatacsvs`folder contains .csv datafiles (see *supplementary information* for the main manuscript to read data file captions). This data is also present on the repository in the form of .jld2 datafiles.

## Tutorial notebooks

To aid reproducability and reusability of this model we are currently writing a series of notebooks that either elucidate our methodology, or directly reproduce results in the paper. These cover the basic ideas behind the model (how the contact rates are generated, the underlying transmission model etc) and the methods used to draw inferences from the data (the interaction between time-since-infection and detectability by different tests) and make predictions (as well as visualising these predictions).

The notebooks are given as `.jmd` files, and require [Weave.jl](https://github.com/JunoLab/Weave.jl) to create html output files. Usage of julia markdown with Weave is very similar to the popular R markdown `.rmd` with `knitr`.   
