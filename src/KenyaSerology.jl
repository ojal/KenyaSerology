__precompile__(true)

module KenyaSerology

using OrdinaryDiffEq,Distributions,Plots,
        Dates,JLD2,Parameters,
        TransformVariables, LogDensityProblems,
        DynamicHMC, DynamicHMC.Diagnostics,MCMCDiagnostics,
        Statistics, Random, MCMCChains,Optim,BlackBoxOptim
using Plots.PlotMeasures
using LinearAlgebra

import ForwardDiff
import AbstractMCMC: AbstractChains

export CoVAreaModel,
        inferparameters!,
        cornerplot,
        plotfittoPCRdata,
        population_plot,
        plot_incidence,
        basic_prior_contactrate_PCR_Peff_mombasa,
        basic_prior_contactrate_PCR_Peff_nairobi,
        basic_prior_contactrate_PCR_Peff_rural

include("datatypes.jl")
include("transmissionmodel.jl")
include("priors.jl")
include("loglikelihoods.jl")
include("observationmodel.jl")
include("plotting.jl")



#Load defaults for the PCR and serological sensitivity after infection
@load("data/default_sero_detection_after_infection.jld2")
@load("data/default_PCR_detection_after_infection.jld2")

end # module
