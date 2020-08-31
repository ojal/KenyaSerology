## Load PATH and relevant packages (most of these aren't required, but I run everything in Main module
# then put into KenyaSerology module)

push!(LOAD_PATH, joinpath(homedir(),"GitHub/Kenya-Serology/src"))
# using OrdinaryDiffEq,Distributions,Plots,
#         Dates,JLD2,Parameters,
#         TransformVariables, LogDensityProblems,
#         DynamicHMC, DynamicHMC.Diagnostics,MCMCDiagnostics,
#         Statistics, Random,Optim
#
# using BlackBoxOptim
# using Plots.PlotMeasures
# using MCMCChains
# import AbstractMCMC: AbstractChains
# using StatsPlots
# import ForwardDiff
# using Revise
# using LinearAlgebra
# using KenyaSerology

using Distributions,Plots,Dates,JLD2,TransformVariables,Optim,FileIO
import KenyaSerology
## Choose plotting backend. I currently like pyplot (which uses the matplotlib python backend) because
#it renders Greek symbols
pyplot()
# plotlyjs()


## Load data and analysis pipeline
@load("data/case_data_by_area_21feb_to_6aug.jld2")
@load("data/cleaned_sero_data_by_area_21feb_to24july.jld2")
province_sero_data = deepcopy(sero_data)
@load("data/cleaned_sero_data_by_area_21feb_to21july.jld2")
@load("data/sero_detection_after_infection_80.jld2")
@load("data/rel_sero_detection_after_infection.jld2")
@load("data/rel_sero_detection_after_infection_with_decay.jld2")
@load("data/PCR_detection_after_infection.jld2")
@load("data/relative_testing_rate.jld2")
@load("data/projected_contact_data_10082020.jld2")
@load("data/death_data_by_area_21feb_to_6aug.jld2")
@load("data/p_ID.jld2")
#Create zero inflated data for sero --- this is because the number of days has to match PCR tests otherwise an error occurs in the likelihood
matchedserodata = vcat(sero_data.serodata,zeros(Int64,length(case_data.dates) - size(sero_data.serodata,1),30,2)) #Match length to case data by inflating with zeros
province_matchedserodata = vcat(province_sero_data.serodata,zeros(Int64,length(case_data.dates) - size(province_sero_data.serodata,1),8,2)) #Match length to case data by inflating with zeros
#Put in the variable transformation
trans_Peff = as((R = asℝ₊,E₀ = asℝ₊,I₀ = asℝ₊,α = asℝ₊,p_test = as𝕀,P_eff = as𝕀))
include("analysis_functions.jl");

## Create fitted priors
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


## List of counties by test result
casesbycounty = vec(sum(case_data.cases,dims = 1))
print(case_data.areas[sortperm(casesbycounty,rev = true)])

serobycounty = vec(sum(sero_data.serodata,dims = [1,3]))
print(sero_data.areas[sortperm(serobycounty,rev = true)])
print(setdiff(case_data.areas,uppercase.(sero_data.areas)))


## Define a model for inference and fit as one pipeline
# The models with short serial intervals 5.5 days - 3.1 latency giving 2.4 days infectious seem to fit better

marsabit_model_Peff_shortsi,marsabit_deaths = generate_model_for_fit("Marsabit",4.6e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
marsabit_model_Peff_shortsi,marsabit_deaths = pipeline_for_fit("Marsabit",4.6e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
bungoma_model_Peff_shortsi,bungoma_deaths = pipeline_for_fit("Bungoma",1.67e6,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
migori_model_Peff_shortsi,migori_deaths = pipeline_for_fit("Migori",1.12e6,KenyaSerology.basic_prior_contactrate_PCR_Peff_semiurban)
wajir_model_Peff_shortsi,wajir_death = pipeline_for_fit("Wajir",7.81e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
tanariver_model_Peff_shortsi,tanariver_death = pipeline_for_fit("Tana River",3.16e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
turkana_model_Peff_shortsi,turkana_death = pipeline_for_fit("Turkana",9.27e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
homabay_model_Peff_shortsi,homabay_death = pipeline_for_fit("Homabay",1.13e6,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
nyeri_model_Peff_shortsi,nyeri_death = pipeline_for_fit("Nyeri",7.6e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
@save("modelfits/nyeri_model_Peff_shortsi.jld2",nyeri_model_Peff_shortsi)
siaya_model_Peff_shortsi,siaya_death = pipeline_for_fit("Siaya",9.93e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
@save("modelfits/siaya_model_Peff_shortsi.jld2",siaya_model_Peff_shortsi)
taitataveta_model_Peff_shortsi,taitataveta_death = pipeline_for_fit("Taita Taveta",3.41e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
@save("modelfits/taitataveta_model_Peff_shortsi.jld2",taitataveta_model_Peff_shortsi)
machakos_model_Peff_shortsi,machakos_death = pipeline_for_fit("Machakos",1.42e6,KenyaSerology.basic_prior_contactrate_PCR_Peff_semiurban,Val(:var_testing))
@save("modelfits/machakos_model_Peff_shortsi.jld2",machakos_model_Peff_shortsi)
kwale_model_Peff_shortsi,kwale_death = pipeline_for_fit("Kwale",8.67e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)

kwale_model_Peff_shortsi,kwale_death = generate_model_for_fit("Kwale",8.67e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
KenyaSerology.inferparameters!(kwale_model_Peff_shortsi,10000,trans_Peff,[1.5,1.,1.,1,0.00001,0.5])
@save("modelfits/kwale_model_Peff_shortsi.jld2",kwale_model_Peff_shortsi)
ugish_model_Peff_shortsi,ugish_death = pipeline_for_fit("Uasin Gishu",1.16e6,KenyaSerology.basic_prior_contactrate_PCR_Peff_semiurban)
@save("modelfits/ugish_model_Peff_shortsi.jld2",ugish_model_Peff_shortsi)
kilifi_model_Peff_shortsi,kilifi_death = pipeline_for_fit("kilifi",1.45e6,KenyaSerology.basic_prior_contactrate_PCR_Peff_semiurban)


## Counties that need strong priors are all rural apart from Nakuru and Isiolo
isiolo_model_Peff_shortsi,nakuru_death = pipeline_for_fit("Isiolo",2.68e5,data_prior_semiurban)
@save("modelfits_inferred/isiolo_model_Peff_shortsi.jld2",isiolo_model_Peff_shortsi)
muranga_model_Peff_shortsi,muranga_death = pipeline_for_fit("Muranga",1.06e6,data_prior_rural)
@save("modelfits_inferred/muranga_model_Peff_shortsi.jld2",muranga_model_Peff_shortsi)
baringo_model_Peff_shortsi,baringo_death = pipeline_for_fit("Baringo",6.67e5,data_prior_rural)
@save("modelfits_inferred/baringo_model_Peff_shortsi.jld2",baringo_model_Peff_shortsi)
embu_model_Peff_shortsi,embu_death = pipeline_for_fit("Embu",6.09e5,data_prior_rural)
@save("modelfits_inferred/embu_model_Peff_shortsi.jld2",embu_model_Peff_shortsi)
elgeyomarakwet_model_Peff_shortsi,elgeyomarakwet_death = pipeline_for_fit("Elgeyo Marakwet",4.54e5,data_prior_rural)
@save("modelfits_inferred/elgeyomarakwet_model_Peff_shortsi.jld2",elgeyomarakwet_model_Peff_shortsi)
laikipia_model_Peff_shortsi,laikipia_death = pipeline_for_fit("Laikipia",5.19e5,data_prior_rural)
@save("modelfits_inferred/laikipia_model_Peff_shortsi.jld2",laikipia_model_Peff_shortsi)
meru_model_Peff_shortsi,meru_death = pipeline_for_fit("Meru",1.55e6,data_prior_rural)
@save("modelfits_inferred/meru_model_Peff_shortsi.jld2",meru_model_Peff_shortsi)
samburu_model_Peff_shortsi,samburu_death = pipeline_for_fit("Samburu",3.1e5,data_prior_rural)
@save("modelfits_inferred/samburu_model_Peff_shortsi.jld2",samburu_model_Peff_shortsi)




## Quick look at parameter posteriors

nairobi_model_Peff_var_testing.MCMC_results.chain
kiambu_model_Peff_shortsi_var_testing.MCMC_results.chain
narok_model_Peff_shortsi.MCMC_results.chain
nyeri_model_Peff_shortsi.MCMC_results.chain





## Saving model fits
@save("modelfits/nairobi_model_Peff_shortsi_var_testing.jld2",nairobi_model_Peff_shortsi_var_testing)
@save("modelfits/nairobi_model_Peff_var_testing.jld2",nairobi_model_Peff_var_testing)
@save("modelfits/nairobi_model_Peff_shortsi.jld2",nairobi_model_Peff_shortsi)
@save("modelfits/mombasa_model_Peff_shortsi.jld2",mombasa_model_Peff_shortsi)
@save("modelfits/mombasa_model_Peff.jld2",mombasa_model_Peff)
@save("modelfits/kiambu_model_Peff_shortsi_var_testing.jld2",kiambu_model_Peff_shortsi_var_testing)
@save("modelfits/narok_model_Peff_shortsi.jld2",narok_model_Peff_shortsi)
@save("modelfits/nyamira_model_Peff_shortsi.jld2",nyamira_model_Peff_shortsi)
@save("modelfits/kericho_model_Peff_shortsi.jld2",kericho_model_Peff_shortsi)
@save("modelfits/garissa_model_Peff_shortsi.jld2",garissa_model_Peff_shortsi)
@save("modelfits/marsabit_model_Peff_shortsi.jld2",marsabit_model_Peff_shortsi)
@save("modelfits/bungoma_model_Peff_shortsi.jld2",bungoma_model_Peff_shortsi)
@save("modelfits/migori_model_Peff_shortsi.jld2",migori_model_Peff_shortsi)
@save("modelfits/wajir_model_Peff_shortsi.jld2",wajir_model_Peff_shortsi)
@save("modelfits/tanariver_model_Peff_shortsi.jld2",tanariver_model_Peff_shortsi)
@save("modelfits/turkana_model_Peff_shortsi.jld2",turkana_model_Peff_shortsi)
@save("modelfits/homabay_model_Peff_shortsi.jld2",homabay_model_Peff_shortsi)



## Plots -- also good for visual confirmation


## These use the MCMCChains plotting methods for Chains type structs
plt_corner_nai = KenyaSerology.cornerplot(nairobi_model_Peff)
plt_corner_mom = KenyaSerology.cornerplot(mombasa_model_Peff)
plt_corner_ugish = KenyaSerology.cornerplot(ugish_model_Peff_shortsi)
# savefig(plt_corner_nai,"draft_corner_plot.png")
plot(nairobi_model.MCMC_results.chain)


## These I've "hand crafted" and are detailed in plotting.jl
plt_swab = KenyaSerology.plotfittoPCRdata(mombasa_model_Peff)
plt_swab = KenyaSerology.plotfittoPCRdata(mombasa_model_Peff_shortsi)
plt_swab = KenyaSerology.plotfittoPCRdata(nyeri_model_Peff_shortsi)
plot!(plt_swab,ylims = (0,20.))
savefig("plotsforpaper/mombasa_model_Peff_shortsi_swab.png")
savefig("plotsforpaper/ugish_model_Peff_shortsi_swab.png")


plt_swab = KenyaSerology.plotfittoPCRdata(nairobi_model_Peff)
plt_swab = KenyaSerology.plotfittoPCRdata(nairobi_model_Peff_shortsi)
plt_swab = KenyaSerology.plotfittoPCRdata(nairobi_model_Peff_shortsi_var_testing)
plt_swab = KenyaSerology.plotfittoPCRdata(kericho_model_Peff_shortsi)

plot!(plt_swab,titlefont = (18,"helvetica"),guidefont = (14,"helvetica"),tickfont = (9,"helvetica"))
plot!(plt_swab,ylims = (0,50.))

savefig("plotsforpaper/nairobi_model_Peff_shortsi_var_testing_swab.png")




plt_pop = KenyaSerology.population_plot(mombasa_model_Peff_shortsi,Val(:monthlyserology))
plt_pop = KenyaSerology.population_plot(turkana_model_Peff_shortsi,Val(:monthlyserology))

savefig("plotsforpaper/mombasa_model_Peff_shortsi_sero.png")
plt_pop = KenyaSerology.population_plot(ugish_model_Peff_shortsi,Val(:monthlyserology))


plt_pop = KenyaSerology.population_plot(nairobi_model_Peff)
plt_pop = KenyaSerology.population_plot(nairobi_model_Peff_shortsi,Val(:monthlyserology))
plt_pop = KenyaSerology.population_plot(marsabit_model_Peff_shortsi,Val(:monthlyserology))
plot!(plt_pop,titlefont = (18,"helvetica"),guidefont = (14,"helvetica"),tickfont = (9,"helvetica"))

savefig("plotsforpaper/nairobi_model_Peff_shortsi_var_testing_sero.png")

plt_pop_weekav = KenyaSerology.population_plot(nairobi_model_Peff,Val(:monthlyserology))

plot!(title = "Nairobi with 'cleaned' serology")
savefig("newnairobi_seroplot.png")


plt_incidence = KenyaSerology.plot_incidence(nairobi_model_Peff_shortsi,nairobi_deaths[1:(end-7)],p_ID)
plt_incidence = KenyaSerology.plot_incidence(nairobi_model_Peff_shortsi_var_testing,nairobi_deaths[1:(end-7)],p_ID)
plot!(plt_incidence,titlefont = (18,"helvetica"),guidefont = (14,"helvetica"),tickfont = (9,"helvetica"))
savefig("plotsforpaper/nairobi_model_Peff_shortsi_var_testing_inc.png")
plt_incidence = KenyaSerology.plot_incidence(mombasa_model_Peff_shortsi,mombasa_deaths[1:(end-3)],p_ID)
savefig("plotsforpaper/mombasa_model_Peff_shortsi_inc.png")
plt_incidence = KenyaSerology.plot_incidence(kiambu_model_Peff_shortsi_var_testing,kiambu_deaths[1:(end-3)],p_ID)
