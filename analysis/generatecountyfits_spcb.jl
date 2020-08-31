## Load PATH and relevant packages (most of these aren't required, but I run everything in Main module
# then put into KenyaSerology module)

push!(LOAD_PATH, joinpath(homedir(),"GitHub/Kenya-Serology/src"))

using Distributions,Plots,Dates,JLD2,TransformVariables,Optim,FileIO
import KenyaSerology
## Choose plotting backend. I currently like pyplot (which uses the matplotlib python backend) because
#it renders Greek symbols
pyplot()
# plotlyjs()



## Load data and analysis pipeline
@load("data/case_data_by_area_21feb_to_6aug.jld2")
@load("data/serologydata_21feb_6aug.jld2")
province_sero_data = deepcopy(sero_data)
@load("data/sero_detection_after_infection_80.jld2")
@load("data/rel_sero_detection_after_infection.jld2")
@load("data/rel_sero_detection_after_infection_with_decay.jld2")
@load("data/PCR_detection_after_infection.jld2")
@load("data/relative_testing_rate.jld2")
@load("data/projected_contact_data_10082020.jld2")
@load("data/death_data_by_area_21feb_to_6aug.jld2")
@load("data/p_ID.jld2")


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
#Create zero inflated data for sero --- this is because the number of days has to match PCR tests otherwise an error occurs in the likelihood
matchedserodata = vcat(sero_data.serodata,zeros(Int64,length(case_data.dates) - size(sero_data.serodata,1),30,2)) #Match length to case data by inflating with zeros
province_matchedserodata = vcat(province_sero_data.serodata,zeros(Int64,length(case_data.dates) - size(province_sero_data.serodata,1),8,2)) #Match length to case data by inflating with zeros
#Put in the variable transformation
trans_Peff = as((R = asℝ₊,E₀ = asℝ₊,I₀ = asℝ₊,α = asℝ₊,p_test = as𝕀,P_eff = as𝕀))
include("analysis_functions.jl");


scatter(case_data.cases[:,case_data.areas.=="NAKURU"])


## List of counties by test result
casesbycounty = vec(sum(case_data.cases,dims = 1))
print(case_data.areas[sortperm(casesbycounty,rev = true)])

serobycounty = vec(sum(sero_data.serodata,dims = [1,3]))
print(sero_data.areas[sortperm(serobycounty,rev = true)])
print(setdiff(case_data.areas,uppercase.(sero_data.areas)))

## Define a model for inference and fit as one pipeline
# The models with short serial intervals 5.5 days - 3.1 latency giving 2.4 days infectious seem to fit better

marsabit_model_Peff_shortsi,marsabit_deaths = pipeline_for_fit("Marsabit",4.6e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
wpokot_model_Peff_shortsi,wpokot_deaths = pipeline_for_fit("West Pokot",6.21e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
busia_model_Peff_shortsi,busia_death = pipeline_for_fit("Busia",8.94e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
vihiga_model_Peff_shortsi,vihiga_deaths = pipeline_for_fit("Vihiga",5.9e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_semiurban)
makueni_model_Peff_shortsi,makueni_deaths = pipeline_for_fit("Makueni",9.88e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
transnzoia_model_Peff_shortsi,transnzoia_deaths = pipeline_for_fit("Trans Nzoia",9.90e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
kisii_model_Peff_shortsi,kisii_deaths = pipeline_for_fit("Kisii",1.27e6,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
lamu_model_Peff_shortsi,lamu_deaths = pipeline_for_fit("Lamu",1.4e5,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
@save("modelfits/lamu_model_Peff_shortsi.jld2",lamu_model_Peff_shortsi)
kisumu_model_Peff_shortsi,kisumu_deaths = pipeline_for_fit("Kisumu",1.56e6,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)

kisumu_model_Peff_shortsi,kisumu_deaths = generate_model_for_fit("Kisumu",1.56e6,KenyaSerology.basic_prior_contactrate_PCR_Peff_semiurban)
KenyaSerology.inferparameters!(kisumu_model_Peff_shortsi,10000,trans_Peff,[1.5,1.,1.,1.,0.00001,0.6])
@save("modelfits/kisumu_model_Peff_shortsi.jld2",kisumu_model_Peff_shortsi)
kilifi_model_Peff_shortsi,kilifi_death = pipeline_for_fit("kilifi",1.45e6,KenyaSerology.basic_prior_contactrate_PCR_Peff_rural)
@save("modelfits/kilifi_model_Peff_shortsi.jld2",kilifi_model_Peff_shortsi)
kajiado_model_Peff_shortsi,kajiado_death = pipeline_for_fit("Kajiado",1.12e6,KenyaSerology.basic_prior_contactrate_PCR_Peff_semiurban,Val(:var_testing))
@save("modelfits/kajiado_model_Peff_shortsi.jld2",kajiado_model_Peff_shortsi)

## Counties that need strong priors are all rural apart from Nakuru and Isiolo
nakuru_model_Peff_shortsi,nakuru_death = pipeline_for_fit("Nakuru",2.16e6,data_prior_semiurban,Val(:var_testing))
@save("modelfits_inferred/nakuru_model_Peff_shortsi.jld2",nakuru_model_Peff_shortsi)
kakamega_model_Peff_shortsi,kakamega_deaths = pipeline_for_fit("Kakamega",1.87e6,data_prior_rural)
@save("modelfits_inferred/kakamega_model_Peff_shortsi.jld2",kakamega_model_Peff_shortsi)
nyandarua_model_Peff_shortsi,nyandarua_deaths = pipeline_for_fit("Nyandarua",6.38e5,data_prior_rural)
@save("modelfits_inferred/nyandarua_model_Peff_shortsi.jld2",nyandarua_model_Peff_shortsi)
bomet_model_Peff_shortsi,bomet_deaths = pipeline_for_fit("Bomet",8.76e5,data_prior_rural)
@save("modelfits_inferred/bomet_model_Peff_shortsi.jld2",bomet_model_Peff_shortsi)
kirinyaga_model_Peff_shortsi,kirinyaga_deaths = pipeline_for_fit("Kirinyaga",6.1e5,data_prior_rural)
@save("modelfits_inferred/kirinyaga_model_Peff_shortsi.jld2",kirinyaga_model_Peff_shortsi)
tharakanithi_model_Peff_shortsi,tharakanithi_deaths = pipeline_for_fit("Tharaka Nithi",3.93e5,data_prior_rural)
@save("modelfits_inferred/tharakanithi_model_Peff_shortsi.jld2",tharakanithi_model_Peff_shortsi)
mandera_model_Peff_shortsi,mandera_deaths = pipeline_for_fit("Mandera",8.67e5,data_prior_rural)
@save("modelfits_inferred/mandera_model_Peff_shortsi.jld2",mandera_model_Peff_shortsi)
nandi_model_Peff_shortsi,nandi_deaths = pipeline_for_fit("Nandi",8.86e5,data_prior_rural)
@save("modelfits_inferred/nandi_model_Peff_shortsi.jld2",nandi_model_Peff_shortsi)
wpokot_model_Peff_shortsi,wpokot_deaths = pipeline_for_fit("West Pokot",6.21e5,data_prior_rural)
@save("modelfits_inferred/wpokot_model_Peff_shortsi.jld2",wpokot_model_Peff_shortsi)

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
@save("modelfits/wpokot_model_Peff_shortsi.jld2",wpokot_model_Peff_shortsi)
@save("modelfits/busia_model_Peff_shortsi.jld2",busia_model_Peff_shortsi)
@save("modelfits/vihiga_model_Peff_shortsi.jld2",vihiga_model_Peff_shortsi)
@save("modelfits/makueni_model_Peff_shortsi.jld2",makueni_model_Peff_shortsi)
@save("modelfits/transnzoia_model_Peff_shortsi.jld2",transnzoia_model_Peff_shortsi)
@save("modelfits/kisii_model_Peff_shortsi.jld2",kisii_model_Peff_shortsi)

kisii_model_Peff_shortsi

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
plt_swab = KenyaSerology.plotfittoPCRdata(kajiado_model_Peff_shortsi)
plot!(plt_swab,ylims = (0,120.))
savefig("plotsforpaper/mombasa_model_Peff_shortsi_swab.png")
savefig("plotsforpaper/ugish_model_Peff_shortsi_swab.png")


plt_swab = KenyaSerology.plotfittoPCRdata(nairobi_model_Peff)
plt_swab = KenyaSerology.plotfittoPCRdata(nairobi_model_Peff_shortsi)
plt_swab = KenyaSerology.plotfittoPCRdata(nairobi_model_Peff_shortsi_var_testing)
plt_swab = KenyaSerology.plotfittoPCRdata(wpokot_model_Peff_shortsi)

plot!(plt_swab,titlefont = (18,"helvetica"),guidefont = (14,"helvetica"),tickfont = (9,"helvetica"))
plot!(plt_swab,ylims = (0,20.))

savefig("plotsforpaper/nairobi_model_Peff_shortsi_var_testing_swab.png")




plt_pop = KenyaSerology.population_plot(mombasa_model_Peff_shortsi,Val(:monthlyserology))
plt_pop = KenyaSerology.population_plot(makueni_model_Peff_shortsi,Val(:monthlyserology))

savefig("plotsforpaper/mombasa_model_Peff_shortsi_sero.png")
plt_pop = KenyaSerology.population_plot(ugish_model_Peff_shortsi,Val(:monthlyserology))


plt_pop = KenyaSerology.population_plot(nairobi_model_Peff)
plt_pop = KenyaSerology.population_plot(nairobi_model_Peff_shortsi,Val(:monthlyserology))
plt_pop = KenyaSerology.population_plot(transnzoia_model_Peff_shortsi,Val(:monthlyserology))
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
