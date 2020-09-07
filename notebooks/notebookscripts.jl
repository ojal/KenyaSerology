push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaSerology/src"))
using OrdinaryDiffEq,Distributions,Plots,
        Dates,JLD2,Parameters,
        TransformVariables, LogDensityProblems,
        DynamicHMC, DynamicHMC.Diagnostics,MCMCDiagnostics,
        Statistics, Random,Optim

# using DiffEqCallbacks

import KenyaSerology


@load("data/case_data_by_area_21feb_to_6aug.jld2");
@load("data/serologydata_21feb_6aug.jld2");
@load("data/death_data_by_area_21feb_to_6aug.jld2");


push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaSerology/src"))#Add the source directory for KenyaSerology to your path
using OrdinaryDiffEq,Distributions,Plots,
        Dates,JLD2,TransformVariables
import KenyaSerology
