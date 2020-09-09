"""
struct MCMCResults

Struct for holding the MCMC chain (a vector of NamedTuples which represent the accepted parameter sets), and the MCMC tree statistics
"""
struct MCMCResults
    chain
    logL::Vector{Float64}
    treestatistics::Vector{DynamicHMC.TreeStatisticsNUTS}
end

"""
mutable struct CoVAreaModel

Struct for holding the fixed data of the area and MCMC results
"""
@with_kw mutable struct CoVAreaModel
    areaname::String #Name of areas
    PCR_cases::VecOrMat{Int64} #Daily positive (and maybe negative) PCR detections
    sero_cases::Matrix{Int64} #Daily positive and negative serology detections
    dates::Vector{Date}#Dates of each day
    N::Float64#Area's population size
    σ::Float64 = 1/3.1 # 1/mean latent period
    γ::Float64 = 1/2.4#Recovery rate
    contactrate_data
    relative_testing_rate::Vector{Float64} = ones(500)
    prob::ODEProblem
    PCR_array::Vector{Float64} = PCR_array #Relative Sensitivity of PCR test on days post-infection
    PCR_sensitivity::Float64 = 1. #Base sensitivity for RT-PCR test
    PCR_specificity::Float64 = 0.995 #Base specificity for RT_PCR test
    sero_array::Vector{Float64} = rel_sero_array_26days #Relative sensitivity of serology test on days post-infection
    sero_sensitivity::Float64 = 0.825 #Base sensitivity for serology test
    sero_specificity::Float64 = 0.992 #Base specificity for serology test
    M_BB::Float64 = 40.74 #Shrinkage factor for BetaBinomial distribution that incorporates the uncertainty in the sensitivity (fitted to serology group)
    log_priors::Function = θ -> 0.
    log_likelihood::Function = θ -> 0.
    MCMC_results::Union{MCMCResults,Nothing} = nothing #This field gets added after the MCMC run
end

"""
mutable struct CoVAreaForecastModel

Struct for holding the fixed data of the area, MCMC results and distributions important for forecasting
"""
@with_kw struct CoVAreaForecastModel
    areamodel::CoVAreaModel #Inference data and MCMC results
    forecastprob::ODEProblem #Potentially a different ODEProblem for forecasting compared to inference
    forecasts
end
