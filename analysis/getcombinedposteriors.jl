push!(LOAD_PATH, joinpath(homedir(),"GitHub/Kenya-Serology/src"))

using Distributions,Plots,Dates,JLD2,TransformVariables,Optim,FileIO,DataFrames,CSV,MAT,StatsPlots
import KenyaSerology

include("analysis_functions.jl");
models_semiurban = ["modelfits/kajiado_model_Peff_shortsi.jld2","modelfits/machakos_model_Peff_shortsi.jld2"]

@load("modelfits/kwale_model_Peff_shortsi.jld2")
models_rural = ["modelfits/tanariver_model_Peff_shortsi.jld2",
                "modelfits/wajir_model_Peff_shortsi.jld2",
                "modelfits/homabay_model_Peff_shortsi.jld2",
                "modelfits/turkana_model_Peff_shortsi.jld2"]

fitted_priors_semiurban = gather_and_create_posterior_fits(models_semiurban)
@save("analysis/fitted_priors_semiurban.jld2",fitted_priors_semiurban)
fitted_priors_rural = gather_and_create_posterior_fits(models_rural)
@save("analysis/fitted_priors_rural.jld2",fitted_priors_rural)
