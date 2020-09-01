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

using Distributions,Plots,Dates,JLD2,TransformVariables
import KenyaSerology
## Choose plotting backend. I currently like pyplot (which uses the matplotlib python backend) because
#it renders Greek symbols
pyplot()
# plotlyjs()


## Load data
@load("data/case_data_by_area_21feb_to_6aug.jld2")
@load("data/cleaned_sero_data_by_area_21feb_to24july.jld2")
province_sero_data = deepcopy(sero_data)
@load("data/cleaned_sero_data_by_area_21feb_to21july.jld2")
@load("data/sero_detection_after_infection_80.jld2")
@load("data/rel_sero_detection_after_infection.jld2")
@load("data/rel_sero_detection_after_infection_with_decay.jld2")
@load("data/PCR_detection_after_infection.jld2")
@load("data/relative_testing_rate.jld2")
@load("data/projected_contact_data_upto30july.jld2")
@load("data/death_data_by_area_21feb_to_24july.jld2")
@load("data/p_ID.jld2")

#Create zero inflated data for sero --- this is because the number of days has to match PCR tests otherwise an error occurs in the likelihood
matchedserodata = vcat(sero_data.serodata,zeros(Int64,length(case_data.dates) - size(sero_data.serodata,1),30,2)) #Match length to case data by inflating with zeros
province_matchedserodata = vcat(province_sero_data.serodata,zeros(Int64,length(case_data.dates) - size(province_sero_data.serodata,1),8,2)) #Match length to case data by inflating with zeros



## Get data for counties to fit ---ADD TO THIS AS REQUIRED

#Swab tests
nairobi_cases = vec(case_data.cases[:,case_data.areas .== "NAIROBI"])
mombasa_cases = vec(case_data.cases[:,case_data.areas .== "MOMBASA"])
ugish_cases = vec(case_data.cases[:,case_data.areas .== "UASIN GISHU"])
kilifi_cases = vec(case_data.cases[:,case_data.areas .== "KILIFI"])
kwale_cases = vec(case_data.cases[:,case_data.areas .== "KWALE"])
kisumu_cases = vec(case_data.cases[:,case_data.areas .== "KISUMU"])
baringo_cases = vec(case_data.cases[:,case_data.areas .== "BARINGO"])

machakos_cases = vec(case_data.cases[:,case_data.areas .== "MACHAKOS"])
taitataveta_cases = vec(case_data.cases[:,case_data.areas .== "TAITA TAVETA"])
kiambu_cases = vec(case_data.cases[:,case_data.areas .== "KIAMBU"])
kajiado_cases = vec(case_data.cases[:,case_data.areas .== "KAJIADO"])
busia_cases = vec(case_data.cases[:,case_data.areas .== "BUSIA"])

#Sero tests
# nairobi_sero = Matrix(sero_data.serodata[:,sero_data.areas.=="nairobi",:][:,1,:])
nairobi_sero = Matrix(province_matchedserodata[:,province_sero_data.areas.=="nairobi",:][:,1,:])
mombasa_sero = Matrix(matchedserodata[:,sero_data.areas.=="mombasa",:][:,1,:])
ugish_sero = Matrix(matchedserodata[:,sero_data.areas.=="uasin gishu",:][:,1,:])
kwale_sero = Matrix(matchedserodata[:,sero_data.areas.=="kwale",:][:,1,:])
kisumu_sero = Matrix(matchedserodata[:,sero_data.areas.=="kisumu",:][:,1,:])
machakos_sero = Matrix(matchedserodata[:,sero_data.areas.=="machakos",:][:,1,:])
kiambu_sero = Matrix(matchedserodata[:,sero_data.areas.=="kiambu",:][:,1,:])
scatter(kiambu_sero)


#Deaths
nairobi_deaths = vec(death_data.deaths[:,death_data.areas .== "NAIROBI"])
mombasa_deaths = vec(death_data.deaths[:,death_data.areas .== "MOMBASA"])
ugish_deaths = vec(death_data.deaths[:,death_data.areas .== "UASIN GISHU"])
kilifi_deaths = vec(death_data.deaths[:,death_data.areas .== "KILIFI"])
kwale_deaths = vec(death_data.deaths[:,death_data.areas .== "KWALE"])
kisumu_deaths = vec(death_data.deaths[:,death_data.areas .== "KISUMU"])
machakos_deaths = vec(death_data.deaths[:,death_data.areas .== "MACHAKOS"])
taitataveta_deaths = vec(death_data.deaths[:,death_data.areas .== "TAITA TAVETA"])
kiambu_deaths = vec(death_data.deaths[:,death_data.areas .== "KIAMBU"])
kajiado_deaths = vec(death_data.deaths[:,death_data.areas .== "KAJIADO"])
busia_deaths = vec(death_data.deaths[:,death_data.areas .== "BUSIA"])


## List of counties by test result
casesbycounty = vec(sum(case_data.cases,dims = 1))
print(case_data.areas[sortperm(casesbycounty,rev = true)])




## Define a model for inference
# The models with short serial intervals 5.5 days - 3.1 latency giving 2.4 days infectious seem to fit better

#This model assumes testing is flat from mid-April
nairobi_model_Peff_shortsi = KenyaSerology.CoVAreaModel(areaname = "Nairobi",
                                PCR_cases = nairobi_cases,
                                sero_cases = nairobi_sero,
                                dates = case_data.dates,
                                N = 4.3e6,
                                γ = 1/(5.5 - 3.1),
                                contactrate_data = projected_contactrate_nairobi,
                                prob = KenyaSerology.make_odeproblemforinference(projected_contactrate_mombasa,
                                                                                        startdate = Date(2020,2,21),
                                                                                        enddate = Date(2020,8,6)),
                                sero_array = rel_sero_array_26days,
                                log_priors = basic_prior_contactrate_PCR_Peff_nairobi,
                                log_likelihood = KenyaSerology.loglikelihood_contactratemodelBB_Peff)

#This model assumes testing varies relative to the whole Kenyan increase in testing rate (which seems to be driven by increased testing in Nairobi)
nairobi_model_Peff_shortsi_var_testing = CoVAreaModel(areaname = "Nairobi",
                                PCR_cases = nairobi_cases,
                                sero_cases = nairobi_sero,
                                dates = case_data.dates,
                                N = 4.3e6,
                                γ = 1/(5.5 - 3.1),
                                contactrate_data = projected_contactrate_nairobi,
                                relative_testing_rate = relative_testing_rate_nairobi.relative_testing_rate, #This is where relative testing rate based on relative to mean testing in Kenya goes in
                                prob = KenyaSerology.make_odeproblemforinference(projected_contactrate_mombasa,
                                                                                        startdate = Date(2020,2,21),
                                                                                        enddate = Date(2020,8,6)),
                                sero_array = rel_sero_array_26days,
                                log_priors = KenyaSerology.basic_prior_contactrate_PCR_Peff_nairobi,
                                log_likelihood = KenyaSerology.loglikelihood_contactratemodelBB_Peff)

nairobi_model_Peff_var_testing = CoVAreaModel(areaname = "Nairobi",
                                PCR_cases = nairobi_cases,
                                sero_cases = nairobi_sero,
                                dates = case_data.dates,
                                N = 4.3e6,
                                γ = 1/9,
                                contactrate_data = projected_contactrate_nairobi,
                                relative_testing_rate = relative_testing_rate_nairobi.relative_testing_rate, #This is where relative testing rate based on relative to mean testing in Kenya goes in
                                prob = KenyaSerology.make_odeproblemforinference(projected_contactrate_mombasa,
                                                                                        startdate = Date(2020,2,21),
                                                                                        enddate = Date(2020,8,6)),
                                sero_array = rel_sero_array_26days,
                                log_priors = KenyaSerology.basic_prior_contactrate_PCR_Peff_nairobi,
                                log_likelihood = KenyaSerology.loglikelihood_contactratemodelBB_Peff)

mombasa_model_Peff_shortsi = KenyaSerology.CoVAreaModel(areaname = "Mombasa",
                                PCR_cases = mombasa_cases,
                                sero_cases = mombasa_sero,
                                dates = case_data.dates,
                                N = 1.2e6,
                                γ = 1/(5.5 - 3.1),
                                contactrate_data = projected_contactrate_mombasa,
                                prob = KenyaSerology.make_odeproblemforinference(projected_contactrate_mombasa,
                                                                                        startdate = Date(2020,2,21),
                                                                                        enddate = Date(2020,8,6)),
                                sero_array = rel_sero_array_26days,
                                log_priors = KenyaSerology.basic_prior_contactrate_PCR_Peff_mombasa,
                                log_likelihood = KenyaSerology.loglikelihood_contactratemodelBB_Peff)

mombasa_model_Peff = KenyaSerology.CoVAreaModel(areaname = "Mombasa",
                                PCR_cases = mombasa_cases,
                                sero_cases = mombasa_sero,
                                dates = case_data.dates,
                                N = 1.2e6,
                                γ = 1/9,
                                contactrate_data = projected_contactrate_mombasa,
                                prob = KenyaSerology.make_odeproblemforinference(projected_contactrate_mombasa,
                                                                                        startdate = Date(2020,2,21),
                                                                                        enddate = Date(2020,8,6)),
                                sero_array = rel_sero_array_26days,
                                log_priors = KenyaSerology.basic_prior_contactrate_PCR_Peff_mombasa,
                                log_likelihood = KenyaSerology.loglikelihood_contactratemodelBB_Peff)


kiambu_model_Peff_shortsi_var_testing = KenyaSerology.CoVAreaModel(areaname = "Kiambu",
                                PCR_cases = ugish_cases,
                                sero_cases = ugish_sero,
                                dates = case_data.dates,
                                N = 2.42e6,
                                γ = 1/(5.5 - 3.1),
                                contactrate_data = projected_contactrate_kenya,
                                relative_testing_rate = relative_testing_rate_nairobi.relative_testing_rate,
                                prob = KenyaSerology.make_odeproblemforinference(projected_contactrate_mombasa,
                                                                                        startdate = Date(2020,2,21),
                                                                                        enddate = Date(2020,8,6)),
                                sero_array = rel_sero_array_26days,
                                log_priors = KenyaSerology.basic_prior_contactrate_PCR_Peff_semiurban,
                                log_likelihood = KenyaSerology.loglikelihood_contactratemodelBB_Peff)



## Inference
#Put in the variable transformation
trans_Peff = as((R = asℝ₊,E₀ = asℝ₊,I₀ = asℝ₊,α = asℝ₊,p_test = as𝕀,P_eff = as𝕀))

#Inference method is in-place. Before inference the MCMC_results field is 'nothing'.
#After inference the MCMC_results field contains a 'Chain' struct (see MCMCChains package) and the tree statistics of the HMC routine
#The Chain object contains draws from the posterior distribution and can used in the TuringLang probabilistic programming ecosystem (worth looking this up)
KenyaSerology.inferparameters!(nairobi_model_Peff_shortsi_var_testing,10000,trans_Peff)
KenyaSerology.inferparameters!(nairobi_model_Peff_var_testing,10000,trans_Peff)
KenyaSerology.inferparameters!(nairobi_model_Peff_shortsi,10000,trans_Peff)
KenyaSerology.inferparameters!(mombasa_model_Peff_shortsi,10000,trans_Peff)
KenyaSerology.inferparameters!(mombasa_model_Peff,10000,trans_Peff)

KenyaSerology.inferparameters!(kiambu_model_Peff_shortsi_var_testing,10000,trans_Peff)


KenyaSerology.modeldic(nairobi_model_Peff_shortsi_var_testing)
KenyaSerology.modeldic(nairobi_model_Peff_var_testing)
KenyaSerology.modeldic(nairobi_model_Peff_shortsi)
KenyaSerology.modeldic(mombasa_model_Peff_shortsi)
KenyaSerology.modeldic(mombasa_model_Peff)
KenyaSerology.modeldic(kiambu_model_Peff_shortsi_var_testing)





#There are now two inferparameters! methods. As before where an initial guess is supplied...
KenyaSerology.inferparameters!(nairobi_model_Peff_shortsi,10000,trans_Peff,[2.5,1.,1.,3.,0.0001,1])
#... and a new method which attempts to find a good initial guess using a "blackbox" global optimiser
KenyaSerology.inferparameters!(mombasa_model_Peff_shortsi,10000,trans_Peff)
#The blackbox search can help, I have also tweaked the initial kinetic energy of the HMC proposal method. But we will still need multiple attempts I think.

## Saving model fits
@save("modelfits/nairobi_model_Peff_shortsi_var_testing.jld2",nairobi_model_Peff_shortsi_var_testing)
@save("modelfits/nairobi_model_Peff_var_testing.jld2",nairobi_model_Peff_var_testing)
@save("modelfits/nairobi_model_Peff_shortsi.jld2",nairobi_model_Peff_shortsi)
@save("modelfits/mombasa_model_Peff_shortsi.jld2",mombasa_model_Peff_shortsi)
@save("modelfits/mombasa_model_Peff.jld2",mombasa_model_Peff)





## Quick look at parameter posteriors

nairobi_model_Peff_var_testing.MCMC_results.chain
kiambu_model_Peff_shortsi_var_testing.MCMC_results.chain

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
plt_swab = KenyaSerology.plotfittoPCRdata(ugish_model_Peff_shortsi)
plot!(plt_swab,ylims = (0,40.))
savefig("plotsforpaper/mombasa_model_Peff_shortsi_swab.png")
savefig("plotsforpaper/ugish_model_Peff_shortsi_swab.png")


plt_swab = KenyaSerology.plotfittoPCRdata(nairobi_model_Peff)
plt_swab = KenyaSerology.plotfittoPCRdata(nairobi_model_Peff_shortsi)
plt_swab = KenyaSerology.plotfittoPCRdata(nairobi_model_Peff_shortsi_var_testing)
plt_swab = KenyaSerology.plotfittoPCRdata(kiambu_model_Peff_shortsi_var_testing)

plot!(plt_swab,titlefont = (18,"helvetica"),guidefont = (14,"helvetica"),tickfont = (9,"helvetica"))
plot!(plt_swab,ylims = (0,1000.))

savefig("plotsforpaper/nairobi_model_Peff_shortsi_var_testing_swab.png")




plt_pop = KenyaSerology.population_plot(mombasa_model_Peff_shortsi,Val(:monthlyserology))
savefig("plotsforpaper/mombasa_model_Peff_shortsi_sero.png")
plt_pop = KenyaSerology.population_plot(ugish_model_Peff_shortsi,Val(:monthlyserology))


plt_pop = KenyaSerology.population_plot(nairobi_model_Peff)
plt_pop = KenyaSerology.population_plot(nairobi_model_Peff_shortsi,Val(:monthlyserology))
plt_pop = KenyaSerology.population_plot(nairobi_model_Peff_shortsi_var_testing,Val(:monthlyserology))
plot!(plt_pop,titlefont = (18,"helvetica"),guidefont = (14,"helvetica"),tickfont = (9,"helvetica"))

savefig("plotsforpaper/nairobi_model_Peff_shortsi_var_testing_sero.png")

plt_pop_weekav = KenyaSerology.population_plot(nairobi_model_Peff,Val(:monthlyserology))

plot!(title = "Nairobi with 'cleaned' serology")
savefig("newnairobi_seroplot.png")


plt_incidence = KenyaSerology.plot_incidence(nairobi_model_Peff_shortsi,nairobi_deaths[1:(end-7)],p_ID)
plt_incidence = KenyaSerology.plot_incidence(nairobi_model_Peff_shortsi_var_testing,nairobi_deaths[1:(end-7)],p_ID)
plot!(plt_incidence,titlefont = (18,"helvetica"),guidefont = (14,"helvetica"),tickfont = (9,"helvetica"))
savefig("plotsforpaper/nairobi_model_Peff_shortsi_var_testing_inc.png")
plt_incidence = KenyaSerology.plot_incidence(mombasa_model_Peff_shortsi,mombasa_deaths[1:(end-7)],p_ID)
savefig("plotsforpaper/mombasa_model_Peff_shortsi_inc.png")
plt_incidence = KenyaSerology.plot_incidence(ugish_model_Peff_shortsi,ugish_deaths[1:(end-7)],p_ID)
