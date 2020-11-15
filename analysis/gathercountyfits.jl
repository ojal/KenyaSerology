## Test script for analysis functions
push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaSerologyPrivate/src"))
using Distributions,Plots,Dates,JLD2,TransformVariables,Optim,FileIO,DataFrames,CSV,MAT
using Revise
import KenyaSerology





## Choose plotting backend. I currently like pyplot (which uses the matplotlib python package) because
#it renders Greek symbols
# pyplot()
# plotlyjs()
gr()

## Load data and analysis pipeline

@load("data/death_data_by_area_21feb_to_17oct.jld2")
@load("data/p_ID.jld2")
@load("data/p_IH.jld2")
@load("data/hosp_data_admission_symptoms_by_area_21feb_to_21aug.jld2")


collecteddata_latesept_not_opt = KenyaSerology.gatherdatafrommodels("modelfits_pos_neg/",death_data,hosp_data,0.00264,p_IH,p_ID)
@save("analysis/collecteddata_lateseptbayesianfit_posneg_not_opt.jld2",collecteddata_latesept_not_opt)

n = (Date(2020,9,30) - Date(2020,2,20)).value


function create_percentage_infected_day_n(collecteddata,n)
	mean_total_infected_on_day_n = [ data.mean_true_infected[n] for data in collecteddata]
	std_total_infected_on_day_n = [ data.std_true_infected[n] for data in collecteddata]
	countynames = [ data.area for data in collecteddata]
	return countynames,mean_total_infected_on_day_n,std_total_infected_on_day_n
end

countynames,mean_total_infected_30sept,std_total_infected_30sept = create_percentage_infected_day_n(collecteddata_latesept_not_opt,n)

df_kenya = DataFrame!(CSV.File("data/2019-population_census-report-per-county.csv"))
df_kenya.County
converted_countynames = copy(df_kenya.County)
converted_countynames[converted_countynames.=="Elgeyo-Marakwet"] .= "Elgeyo Marakwet"
converted_countynames[converted_countynames.=="Murang'a"] .= "Muranga"
converted_countynames[converted_countynames.=="Tharaka-Nithi"] .= "Tharaka Nithi"
converted_countynames[converted_countynames.=="Homabay"] .= "HomaBay"

N = df_kenya.Total_Population19

mean_prop_total_infected_30sept = mean_total_infected_30sept./N
std_prop_total_infected_30sept = std_total_infected_30sept./N

matwrite("prop_infect_by_county.mat", Dict(
	"countynames" => countynames,
	"mean_perc_fitted" => mean_prop_total_infected_30sept,
	"std_perc_fitted" => std_prop_total_infected_30sept
); compress = true)



## Output the headline fits to a CSV file



modelnames = readdir("modelfits_pos_neg_opt/")
modelpaths = [joinpath("modelfits_pos_neg_opt/",name) for name in modelnames ]

df1 = KenyaSerology.parameterinferenceovercollection(modelpaths[1:23],death_data,0.00264,p_ID;num_sims = 1000)
df2 = KenyaSerology.parameterinferenceovercollection(modelpaths[23:end],death_data,0.00264,p_ID;num_sims = 1000)

# @load("analysis/collecteddata_lateseptbayesianfit_posneg.jld2")
CSV.write("countyparameters_fittedCt_1.csv",df1)
CSV.write("countyparameters_fittedCt_2.csv",df2)

modelnames = readdir("modelfits_pos_neg/")
modelpaths = [joinpath("modelfits_pos_neg/",name) for name in modelnames ]
df1 = KenyaSerology.parameterinferenceovercollection(modelpaths[1:23],death_data,0.00264,p_ID;num_sims = 1000)
df2 = KenyaSerology.parameterinferenceovercollection(modelpaths[24:end],death_data,0.00264,p_ID;num_sims = 1000)

# @load("analysis/collecteddata_lateseptbayesianfit_posneg.jld2")
CSV.write("countyparameters_googleCt_1.csv",df1)
CSV.write("countyparameters_googleCt_2.csv",df2)

## Get DIC for county fits

modelnames = readdir("modelfits_pos_neg/")
modelpaths = [joinpath("modelfits_pos_neg/",name) for name in modelnames ]
countynames = String[]
DIC = Float64[]
for (i,path) in enumerate(modelpaths)
	D = load(path)
	model = D[first(keys(D))]
	push!(countynames,model.areaname)
	push!(DIC,KenyaSerology.modeldic(model))
end
modelnames = readdir("modelfits_pos_neg_opt/")
modelpaths = [joinpath("modelfits_pos_neg_opt/",name) for name in modelnames ]
countynames_opt = String[]
DIC_opt = Float64[]
for (i,path) in enumerate(modelpaths)
	D = load(path)
	model = D[first(keys(D))]
	push!(countynames_opt,model.areaname)
	push!(DIC_opt,KenyaSerology.modeldic(model))
end

df = DataFrame(county=countynames,DIC_googleCt=DIC,DIC_fittedCt=DIC_opt)
CSV.write("DIC_for_models.csv",df)

## LPD for each county

function get_LPD(areamodel::KenyaSerology.CoVAreaModel,deaths::Vector{Int64},p_ID)
    inc = KenyaSerology.incidence_across_samples(areamodel,315);
	x = KenyaSerology.simple_conv(inc.true_incidence[:,1],p_ID)
	unscaleddeaths = zeros(length(x),size(inc.true_incidence,2))
	for j = 1:size(unscaleddeaths,2)
		unscaleddeaths[:,j] .= KenyaSerology.simple_conv(inc.true_incidence[:,j],p_ID)
	end
    D = sum(deaths)
    μ_ests = [Erlang(D+1,1/((1/0.00264) + sum(unscaleddeaths[1:length(deaths),n]))) for n = 1:size(inc.true_incidence,2)]
	return KenyaSerology.logdensity_deaths(deaths,mean.(μ_ests),unscaleddeaths)
end

@load("data/death_data_by_area_21feb_to_17oct.jld2")
@load("data/p_ID.jld2")
modelnames = readdir("modelfits_pos_neg/")
modelpaths = [joinpath("modelfits_pos_neg/",name) for name in modelnames ]
countynames = String[]
LPD_opt = Float64[]
for (i,path) in enumerate(modelpaths)
	D = load(path)
	model = D[first(keys(D))]
	deaths = vec(death_data.deaths[:,death_data.areas .== uppercase(model.areaname)][1:223])
	println("Starting on county $(model.areaname)")
	push!(countynames_opt,model.areaname)
	push!(LPD_opt,get_LPD(model,deaths,p_ID))
end
LPD = deepcopy(LPD_opt)

modelnames = readdir("modelfits_pos_neg_opt/")
modelpaths = [joinpath("modelfits_pos_neg_opt/",name) for name in modelnames ]
countynames_opt = String[]
LPD_opt = Float64[]
for (i,path) in enumerate(modelpaths)
	D = load(path)
	model = D[first(keys(D))]
	deaths = vec(death_data.deaths[:,death_data.areas .== uppercase(model.areaname)][1:223])
	println("Starting on county $(model.areaname)")
	push!(countynames_opt,model.areaname)
	push!(LPD_opt,get_LPD(model,deaths,p_ID))
end
df = DataFrame(county=countynames,DIC_googleCt=DIC,DIC_fittedCt=DIC_opt,LPD_googleCt = LPD,LPD_fittedCt = LPD_opt)
CSV.write("fit_scores_for_models.csv",df)


##
@load("modelfits_pos_neg/Nairobi_model.jld2")
inc = KenyaSerology.incidence_across_samples(model,300)

true_inf = KenyaSerology.create_credible_intervals(inc.true_infecteds)

plot(100*true_inf.mean_pred./model.N)



# ## Gather fitted data and save or load as required
# collecteddata_fitted = gatherdatafrommodels("modelfits/",death_data,0.00264,p_ID)
# collecteddata_fittedlateaug = gatherdatafrommodels("modelfitslateaug/",death_data,hosp_data,0.00264,p_IH,p_ID)
#
# # @save("analysis/collecteddata_fitted_bayesianfit.jld2",collecteddata_fitted)
# # @save("analysis/collecteddata_fittedlateaug_bayesianfit.jld2",collecteddata_fittedlateaug)
#
# @load("analysis/collecteddata_fittedlateaug_bayesianfit.jld2",collecteddata_fittedlateaug)
# @load("analysis/collecteddata_data_inferredlateaug_bayesianfit.jld2",collecteddata_datalateaug_inferred)
#
# collecteddata_data_inferred = gatherdatafrommodels("modelfits_inferred/",death_data,0.00264,p_ID)
# collecteddata_datalateaug_inferred = gatherdatafrommodels("modelfitslateaug_inferred/",death_data,hosp_data,0.00264,p_IH,p_ID)
# # @save("analysis/collecteddata_data_inferred_bayesianfit.jld2",collecteddata_data_inferred)
# # @save("analysis/collecteddata_data_inferredlateaug_bayesianfit.jld2",collecteddata_datalateaug_inferred)

## Collect data on peak timing for export to MATLAB peak plotting by county



peaktimes_fitted = createpeaktimes(collecteddata_fitted,countynames)

matwrite("plotsforpaper/fittedpeaktimesbycounty.mat", Dict(
	"countynames" => countynames,
	"peaktimes_fitted" => peaktimes_fitted
); compress = true)

peaktimes_data_inferred = createpeaktimes(collecteddata_data_inferred,countynames)


matwrite("plotsforpaper/datainferredpeaktimesbycounty.mat", Dict(
	"countynames" => countynames,
	"peaktimes_data_inferred" => peaktimes_data_inferred
); compress = true)

i = 100
f_R = findfirst(keys(model.MCMC_results.chain) .== :R)
y = model.MCMC_results.chain[100,f_R,1]
model.MCMC_results.chain.name_map.parameters
## Analysis for mean R outside of Nairobi and Mombasa

f = findall([(d.area != "Mombasa")&&(d.area != "Nairobi") for d in collecteddata_fitted])
mean_R_over_counties_array = [d.mean_R for d in collecteddata_fitted[f]]
histogram(mean_R_over_counties_array,bins = 20)
mean_R_over_counties = mean(mean_R_over_counties_array)

## Create CSV file with fits for each county for each parameter, peak, IFR_crude

df = parameterinferenceovercollection("modelfits/",death_data,0.00264,p_ID)
df.weakpriors = trues(30)

df_sps = parameterinferenceovercollection("modelfits_inferred/",death_data,0.00264,p_ID)
df_sps.weakpriors = falses(17)

df_total = vcat(df,df_sps)


CSV.write("analysis/parameterfitsforeachcounty_vs2.csv",df_total)


## Respond to information requests
pyplot()
dates = [Date(2020,2,20) + Day(k) for k =1:588]

xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,20)).value for k = 1:20 ]
xticklabs2020 = [(monthname(k)[1:3])*"-20" for k = 3:12]
xticklabs2021 = [(monthname(k)[1:3])*"-21" for k = 1:10]
xticklabs = vcat(xticklabs2020,xticklabs2021)

projectionsfittedmodels = incidence_severe_death_projectionovercollection("modelfitslateaug/",hosp_data,death_data,0.00264,p_IH,p_ID)
projectionsinferredmodels = incidence_severe_death_projectionovercollection("modelfitslateaug_inferred/",hosp_data,death_data,0.00264,p_IH,p_ID)
# @save("analysis/projectionsfittedmodels.jld2",projectionsfittedmodels)
# @save("analysis/projectionsinferredmodels.jld2",projectionsinferredmodels)
@load("analysis/projectionsfittedmodels.jld2",projectionsfittedmodels)
@load("analysis/projectionsinferredmodels.jld2",projectionsinferredmodels)



kenyan_inc = plot(projectionsfittedmodels.totalkenyanincidence .+ projectionsinferredmodels.totalkenyanincidence.+1,
		yscale = :log10,
		xticks = (xticktimes,xticklabs),
		yticks = ([0,1,10,100,1000,10000,100000,1000000].+1,[0,1,10,100,1000,10000,100000,1000000]),
		tickfont = 6,
		title = "Total Kenya incidence",
		lab = "forecast incidence",
		ylabel = "Daily new infections",
		size = (700,500),dpi=250)
daylessthan100 = findfirst((projectionsfittedmodels.totalkenyanincidence .+ projectionsinferredmodels.totalkenyanincidence)[31:end] .<= 100)+30
plot!(kenyan_inc,[daylessthan100,daylessthan100],[1,10000],lw=2,
		ls = :dash,lab = "Less than 100 daily new infections on $(Date(2020,2,20)+Day(daylessthan100))")
savefig("plotsforreport/kenyan_incidence_to_end.pdf")

kenyan_hosp = plot(projectionsfittedmodels.totalkenyanhosp .+ projectionsinferredmodels.totalkenyanhosp,
		xticks = (xticktimes,xticklabs),
		tickfont = 6,
		title = "Total Kenya Hospitalisations --- fitted upto 17th July",
		lab = "forecast hospitalisation",
		ylabel = "Daily new hospitalisation",
		size = (700,500),dpi=250)
scatter!(sum(hosp_data.hosps,dims = 2),ms=5,lab = "Hospitalisation data")
fittingdate = (Date(2020,7,17)-Date(2020,2,20)).value
plot!([fittingdate,fittingdate],[0,20],lw = 2,ls = :dash,lab = "Fitting window")
savefig("plotsforreport/kenyan_hospitalisation_to_end.pdf")



@load("analysis/projectionsfittedmodels.jld2",projectionsfittedmodels)
@load("analysis/projectionsinferredmodels.jld2",projectionsinferredmodels)




## Edwine's questions

kenyadata = combinedata(vcat(collecteddata_fittedlateaug,collecteddata_datalateaug_inferred))
kenyahosps = projectionsfittedmodels.totalkenyanhosp .+ projectionsinferredmodels.totalkenyanhosp

# Duration of the pandemic (at the end of the epi curve- no more susceptible)
# Cumulative number of people infected (at the end of the epi curve- no more susceptible)
cumsum(kenyadata.kenyaincidence)[end],3*sqrt.(sum(kenyadata.kenyaincidence_std.^2))

# Cumulative number of people hospitalized (severe disease) (at the end of the epi curve- no more susceptible)
cumsum(kenyahosps)[end],invlogcdf(Poisson(cumsum(kenyahosps)[end]),log(0.025)),invlogcdf(Poisson(cumsum(kenyahosps)[end]),log(0.975))

# Cumulative number of deaths (at the end of the epi curve- no more susceptible)
cumsum(kenyadata.kenyadeaths)[end],invlogcdf(Poisson(cumsum(kenyadata.kenyadeaths)[end]),log(0.025)),invlogcdf(Poisson(cumsum(kenyadata.kenyadeaths)[end]),log(0.975))

## Forecasting requirements

hosp_kenya = projectionsfittedmodels.totalkenyanhosp .+ projectionsinferredmodels.totalkenyanhosp

#Discrete distribution of infection to symptoms
d_incubation = LogNormal(1.644,0.363)#Lauer estimate
d_duration_in_hosp = Weibull(1.572,17.819)#Haw et al
d_ICUstay = Gamma(2.66,3.42)#Fit to IQR
d_recovery = Exponential(2.4)
#lag distributions
p_IS = [cdf(d_incubation,t) - cdf(d_incubation,t-1) for t = 1:100]#infection to symptoms
p_ICU = [cdf(d_ICUstay,t) - cdf(d_ICUstay,t-1) for t = 1:100]#ICU to leave ICU
p_HR = [cdf(d_duration_in_hosp,t) - cdf(d_duration_in_hosp,t-1) for t = 1:100]#Hospital to discharge
p_ICUR = KenyaSerology.simple_conv(p_ICU,p_HR)#ICU to discharge assuming its the sum of the two
p_SH = [0.2 for t = 1:5] #Moghadas et al uniform 1-5 day estimate
p_R = [cdf(d_recovery,t) - cdf(d_recovery,t-1) for t = 1:1000]

#Upper tail functions
Q_HR = vcat([1.],[1 - cumsum(p_HR)[t] for t = 1:100])
Q_ICUH = vcat([1.],[1 - cumsum(p_ICU)[t] for t = 1:100])
Q_ICUR = vcat([1.],[1 - cumsum(p_ICUR)[t] for t = 1:100])
Q_R = vcat([1.],[1 - cumsum(p_R)[t] for t = 1:1000])
F_R = 1 .- Q_R
#Parameters
symptomaticrate = 0.0657 #From linelist
prop_critical = 0.21

#Incidence results
asymp_incidence = (1-symptomaticrate)*KenyaSerology.simple_conv(kenyadata.kenyaincidence,p_IS)
symp_incidence = symptomaticrate*KenyaSerology.simple_conv(kenyadata.kenyaincidence,p_IS)
severe_incidence = (1-prop_critical)*hosp_kenya
crit_incidence = prop_critical*hosp_kenya
death_incidence = kenyadata.kenyadeaths
#prevalence results
asymp_prev = KenyaSerology.simple_conv(asymp_incidence,Q_R)
critical_prev = KenyaSerology.simple_conv(crit_incidence,Q_ICUH)
severe_prev = KenyaSerology.simple_conv(severe_incidence,Q_HR) .+ KenyaSerology.simple_conv(crit_incidence,Q_ICUR) .- critical_prev
mild_prev = KenyaSerology.simple_conv(symp_incidence,Q_R) .- severe_prev .- critical_prev
recovered_prev = KenyaSerology.simple_conv(asymp_incidence.+symp_incidence,F_R)


"""

# Arguments
- `n::Integer`: the number of elements to compute.
- `dim::Integer=1`: the dimensions along which to perform the computation.
"""
function predict_incidence_and_prev(;incidence,
										hospitalisations,deaths,
										symptomaticrate,
										prop_critical,
										p_IS,
										Q_R,
										Q_ICUH,
										F_R)
	asymp_incidence = (1-symptomaticrate)*KenyaSerology.simple_conv(incidence,p_IS)
	symp_incidence = symptomaticrate*KenyaSerology.simple_conv(incidence,p_IS)
	severe_incidence = (1-prop_critical)*hospitalisations
	crit_incidence = prop_critical*hospitalisations
	death_incidence = deaths
	asymp_prev = KenyaSerology.simple_conv(asymp_incidence,Q_R)
	critical_prev = KenyaSerology.simple_conv(crit_incidence,Q_ICUH)
	severe_prev = KenyaSerology.simple_conv(severe_incidence,Q_HR) .+ KenyaSerology.simple_conv(crit_incidence,Q_ICUR) .- critical_prev
	mild_prev = KenyaSerology.simple_conv(symp_incidence,Q_R) .- severe_prev .- critical_prev
	recovered_prev = KenyaSerology.simple_conv(asymp_incidence.+symp_incidence,F_R)
	return (asymp_incidence = asymp_incidence,
	symp_incidence = symp_incidence,
	severe_incidence = severe_incidence,
	crit_incidence = crit_incidence,
	death_incidence = death_incidence,
	asymp_prev = asymp_prev,
	critical_prev = critical_prev,
	severe_prev = severe_prev,
	mild_prev = mild_prev,
	recovered_prev = recovered_prev)
end

function predict_incidence_and_prev(modelfit,symptomaticrate,
										prop_critical,
										p_IS,
										Q_R,
										Q_ICUH,
										F_R)
	incidence = modelfit.mean_pred_incidence
	hospitalisations = modelfit.pred_hosps_ftc
	deaths = modelfit.pred_deaths_ftc
	name = modelfit.area
	return predict_incidence_and_prev(incidence,hospitalisations,deaths,
										symptomaticrate,
										prop_critical,
										p_IS,
										Q_R,
										Q_ICUH,
										F_R),name
end

predictions,name = predict_incidence_and_prev(collecteddata_datalateaug_inferred[1],symptomaticrate,
										prop_critical,
										p_IS,
										Q_R,
										Q_ICUH,
										F_R)


function outputcountypredictions(outputdirname,modelfitarray,symptomaticrate,
										prop_critical,
										p_IS,
										Q_R,
										Q_ICUH,
										F_R)
	for modelfit in modelfitarray
		predictions,name = predict_incidence_and_prev(modelfit,symptomaticrate,
										prop_critical,
										p_IS,
										Q_R,
										Q_ICUH,
										F_R)
		dates = string.(collect(Date(2020,3,13):Day(1):Date(2021,9,30)))
		df = DataFrame(Date = dates,
				inc_infection = (predictions.asymp_incidence.+predictions.symp_incidence)[22:end],
 				prev_infection = (predictions.asymp_prev.+predictions.critical_prev.+predictions.severe_prev.+predictions.mild_prev)[22:end],
				inc_asymptomatic = predictions.asymp_incidence[22:end],
				prev_asymptomatic = predictions.asymp_prev[22:end],
				inc_mild =predictions.symp_incidence[22:end],
				prev_mild = predictions.mild_prev[22:end],
				inc_severe = predictions.severe_incidence[22:end],
				prev_severe = predictions.severe_prev[22:end],
				inc_critical = predictions.crit_incidence[22:end],
				prev_critical = predictions.critical_prev[22:end],
				inc_death = predictions.death_incidence[22:end],
				prev_recovereds =  predictions.recovered_prev[22:end])

		CSV.write(joinpath(outputdirname,"outputforreport_$(name).csv"),df)
	end
	return nothing
end

outputcountypredictions("outputforreport",collecteddata_fittedlateaug,symptomaticrate,
										prop_critical,
										p_IS,
										Q_R,
										Q_ICUH,
										F_R)

outputcountypredictions("outputforreport",collecteddata_datalateaug_inferred,symptomaticrate,
										prop_critical,
										p_IS,
										Q_R,
										Q_ICUH,
										F_R)
## Output to CSV file for C&P
df = DataFrame(inc_infection = (asymp_incidence.+symp_incidence)[22:end],
 				prev_infection = (asymp_prev.+critical_prev.+severe_prev.+mild_prev)[22:end],
				inc_asymptomatic = asymp_incidence[22:end],
				prev_asymptomatic = asymp_prev[22:end],
				inc_mild =symp_incidence[22:end],
				prev_mild = mild_prev[22:end],
				inc_severe = severe_incidence[22:end],
				prev_severe = severe_prev[22:end],
				inc_critical = crit_incidence[22:end],
				prev_critical = critical_prev[22:end],
				inc_death = death_incidence[22:end],
				prev_recovereds =  recovered_prev[22:end])

CSV.write("plotsforreport/outputforreport.csv",df)

kenya_hosp_prev = plot(severe_prev,
		xticks = (xticktimes,xticklabs),
		xtickfont = 6,
		ytickfont = 10,
		lw=2,
		title = "Total Kenya hospital prevalence",
		lab = "Numbers in hospital (not ICU)",
		ylabel = "Numbers on day",
		size = (700,500),dpi=250)
plot!(kenya_hosp_prev,critical_prev,
		lw=2,
		lab = "Numbers in ICU")
savefig("plotsforreport/numbersinhospital.pdf")

#Columns needed for data entry
#DATE --- 13/03/2020 (day 22) 30/09/2021 (day 588)
# incident infection
# prevalent infection
# incident asymptomatic
# prevalent asymptomatic
# incident mild
# prevalent mild
# incident severe
# prevalent severe
# incident critical
# prevalent critical
# incident death
# prevalent recovered
# total (check)
# total infectious (check)
