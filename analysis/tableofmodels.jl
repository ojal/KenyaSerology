push!(LOAD_PATH, joinpath(homedir(),"GitHub/Kenya-Serology/src"))

using Distributions,Plots,Dates,JLD2,TransformVariables,Optim,FileIO,DataFrames,CSV,MAT,StatsPlots
using Revise
import KenyaSerology
include("analysis_functions.jl");
include("getdata.jl");
nai_death = vec(death_data.deaths[:,death_data.areas .== "NAIROBI"])
mom_death = vec(death_data.deaths[:,death_data.areas .== "MOMBASA"])


@load("data/projected_contact_data_10082020.jld2")
@load("modelfits/nairobi_model_Peff_shortsi_var_testing.jld2")


##Doubling time calculations
SEIR_r(R₀,σ,γ) = 0.5*(sqrt(4*(R₀-1)*σ*γ + (σ+γ)^2) - (σ+γ))
SEIR_r(1.17,1/3.1,1/(5.5-3.1))
log(2)/SEIR_r(1.17,1/3.1,1/(5.5-3.1))
log(2)/SEIR_r(1.35,1/3.1,1/(5.5-3.1))

log(2)/SEIR_r(2.23*0.527,1/3.1,1/(5.5-3.1))
log(2)/SEIR_r(2.01*0.518,1/3.1,1/(5.5-3.1))
2.23*0.527
2.01*0.518
plot(projected_contactrate_mombasa.contactrate)
(mincontactrate, day) = findmin(projected_contactrate_nairobi.contactrate)
Date(2020,2,20) + Day(day[1])

@load("modelfits/nairobi_model_Peff_shortsi_var_testing.jld2")

KenyaSerology.modeldic(nairobi_model_Peff_shortsi_var_testing)
table_data_nairobi_model_Peff_shortsi_var_testing = datafortable(nairobi_model_Peff_shortsi_var_testing,nai_death)
round.(table_data_nairobi_model_Peff_shortsi_var_testing.IFR_CI,sigdigits = 3)

@load("modelfitsforDIC/nairobi_model_Peff_shortsi.jld2")

KenyaSerology.modeldic(nairobi_model_Peff_shortsi)
table_data_nairobi_model_Peff_shortsi = datafortable(nairobi_model_Peff_shortsi,nai_death)

@load("modelfitsforDIC/nairobi_model_Peff_var_testing.jld2")
KenyaSerology.modeldic(nairobi_model_Peff_var_testing)
table_data_nairobi_model_Peff_var_testing = datafortable(nairobi_model_Peff_var_testing,nai_death)




@load("modelfits/mombasa_model_Peff_shortsi.jld2")
mean(mombasa_model_Peff_shortsi.MCMC_results.chain[:,2,1]) + mean(mombasa_model_Peff_shortsi.MCMC_results.chain[:,3,1])

KenyaSerology.modeldic(mombasa_model_Peff_shortsi)

table_data_mombasa_model_Peff_shortsi = datafortable(mombasa_model_Peff_shortsi,mom_death)


@load("modelfitsforDIC/mombasa_model_Peff.jld2")
KenyaSerology.modeldic(mombasa_model_Peff)
table_data_mombasa_model_Peff= datafortable(mombasa_model_Peff,mom_death)
function basic_prior_twoR_nairobi(θ)
        E₀,I₀,R_eff1,R_eff2,α,p_test = θ
        return logpdf(Gamma(1,100/1),E₀ + I₀) +
                logpdf(Gamma(2,2.5/2),R_eff1) +
                logpdf(Gamma(2,1.5/2),R_eff2) +
                logpdf(Gamma(10,3000/4.3e7),p_test)
end

nairobi_model_twoR = KenyaSerology.CoVAreaModel(areaname = "Nairobi",#Name
                                PCR_cases = nairobi_cases,#The PCR tests determined as positive on each day as a vector, first element of PCR cases has to be Feb-21st
                                sero_cases = nairobi_model_Peff_shortsi_var_testing.sero_cases,#The sero-pos tests determined on each day as a Matrix, first row has to be Feb-21st. Col one is pos sero tests, col two is neg sero tests
                                dates = case_data.dates,#Dates as a vector (not actually used)
                                N = 4.3e6,#Population size
                                γ = 1/(5.5-3.1),#Infectious period
                                contactrate_data = contactrate_nairobi,
                                prob = KenyaSerology.make_odeproblemforinference(startdate = Date(2020,2,21),#Don't change from Feb 21st as start date!
                                                                                    changedate = Date(2020,4,1),
                                                                                        enddate = Date(2020,8,6)),#This caps when inference occurs until
                                sero_array = sero_array_80,#IMPORTANT This determines the sensitivity of the sero assay as a function of time since infection, there is also a field for PCR sensitivity but this defaults to Zhou et al if not used
                                log_priors = basic_prior_twoR_nairobi,#Put in your choice of log-prior here
                                log_likelihood = KenyaSerology.loglikelihoodBBtwoR)

trans_twoR = as((R_eff1 = asℝ₊,R_eff2 = asℝ₊,E₀ = asℝ₊,I₀ = asℝ₊,α = asℝ₊,p_test = as𝕀))

f(x) = -nairobi_model_twoR((R_eff1=x[1],R_eff2=x[2],E₀=x[3],I₀ = x[4],α=x[5],p_test=x[6]))
searchrange = [(1.,4.),(1.,2.),(0.,10.),(0.,10.),(0.,3.),(0.00000,0.001)]
using BlackBoxOptim
res = bboptimize(f; SearchRange = searchrange,PopulationSize=100)

KenyaSerology.inferparameters!(nairobi_model_twoR,10000,trans_twoR,best_candidate(res))
prob = KenyaSerology.make_odeproblemforinference(startdate = Date(2020,2,21),#Don't change from Feb 21st as start date!
                                                    changedate = Date(2020,4,1),
                                                        enddate = Date(2020,8,6))
