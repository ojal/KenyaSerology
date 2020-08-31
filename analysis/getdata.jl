@load("data/p_ID.jld2")
@load("data/death_data_by_area_21feb_to_6aug.jld2")
countydata = DataFrame!(CSV.File("data/2019-population_census-report-per-county.csv"))
countynames = sort!(countydata.County)
countynames[5] = "Elgeyo Marakwet"
countynames[8] = "Homabay"
countynames[29] = "Muranga"
countynames[41] = "Tharaka Nithi"
## Create fitted priors --- defined on workspace ahead of loading model fits
@load("analysis/fitted_priors_rural.jld2")
@load("analysis/fitted_priors_semiurban.jld2")
function data_prior_rural(θ)
        E₀,I₀,R,α,p_test,P_eff = θ
        return logpdf(fitted_priors_rural.d_E₀,E₀) +
                logpdf(fitted_priors_rural.d_I₀,I₀) +
                logpdf(fitted_priors_rural.d_Peff,P_eff) +
                logpdf(fitted_priors_rural.d_R,R) +
                logpdf(fitted_priors_rural.d_ptest,p_test) +
                logpdf(fitted_priors_rural.d_α,α)
end
function data_prior_semiurban(θ)
        E₀,I₀,R,α,p_test,P_eff = θ
        return logpdf(fitted_priors_semiurban.d_E₀,E₀) +
                logpdf(fitted_priors_semiurban.d_I₀,I₀) +
                logpdf(fitted_priors_semiurban.d_Peff,P_eff) +
                logpdf(fitted_priors_semiurban.d_R,R) +
                logpdf(fitted_priors_semiurban.d_ptest,p_test) +
                logpdf(fitted_priors_semiurban.d_α,α)
end
