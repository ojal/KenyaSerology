push!(LOAD_PATH, joinpath(homedir(),"GitHub/Kenya-Serology/src"))
using Distributions,Plots,Dates,JLD2,TransformVariables,Optim,FileIO,DataFrames,CSV,MAT
import KenyaSerology

include("analysis_functions.jl");
include("getdata.jl");

collecteddata_fitted = gatherdatafrommodels("modelfits/",death_data,0.00264,p_ID)
@save("analysis/collecteddata_fitted.jld2",collecteddata_fitted)
# @load("analysis/collecteddata_fitted.jld2",collecteddata_fitted)

peaktimes_fitted = createpeaktimes(collecteddata_fitted,countynames)

matwrite("plotsforpaper/fittedpeaktimesbycounty.mat", Dict(
	"countynames" => countynames,
	"peaktimes_fitted" => peaktimes_fitted
); compress = true)

collecteddata_data_inferred = gatherdatafrommodels("modelfits_inferred/",death_data,0.00264,p_ID)
peaktimes_data_inferred = createpeaktimes(collecteddata_data_inferred,countynames)
@save("analysis/collecteddata_data_inferred.jld2",collecteddata_data_inferred)

matwrite("plotsforpaper/datainferredpeaktimesbycounty.mat", Dict(
	"countynames" => countynames,
	"peaktimes_data_inferred" => peaktimes_data_inferred
); compress = true)


## Analysis for mean R outside of Nairobi and Mombasa
f = findall([(d.area != "Mombasa")&&(d.area != "Nairobi") for d in collecteddata_fitted])
mean_R_over_counties_array = [d.mean_R for d in collecteddata_fitted[f]]
histogram(mean_R_over_counties_array,bins = 20)
mean_R_over_counties = mean(mean_R_over_counties_array)

## Create CSV file with fits for each county for each parameter, peak, IFR_crude, posterior predictive p-value using test statistic T(X) = ∑(X_n+1 - X_n)^2

@load("modelfits/nairobi_model_Peff_shortsi_var_testing.jld2")
plot(diff(death_data.deaths[:,death_data.areas .=="NAIROBI"],dims=1))

nairobi_model_Peff_shortsi_var_testing.MCMC_results.chain
inc_nai = KenyaSerology.incidence_across_samples(nairobi_model_Peff_shortsi_var_testing,315)
deaths_nai = vec(death_data.deaths[:,death_data.areas .=="NAIROBI"])
T_data_nai = sum(diff(deaths_nai).^2) #This is the test statistic


datafornai = getdataforCSVfile(nairobi_model_Peff_shortsi_var_testing,deaths_nai,0.00264,p_ID)

df = DataFrame(Countyname = String[],
				R0 = String[],
				E0 = String[],
				I0 = String[],
				clusteringfactor = String[],
				p_test = String[],
				EffPopSize = String[],
				Infectionratepeak = String[],
				IFR = String[])

# function create				







D = sum(deaths_nai[1:165])
IFR_ests = mean.([Erlang(D+1,1/((1/0.00265) + sum(KenyaSerology.simple_conv(inc_nai.true_incidence[:,n],p_ID)[1:165]))) for n = 1: 10000])
v = IFR_ests[1]*KenyaSerology.simple_conv(inc_nai.true_incidence[:,1],p_ID)[1:165]

T_array = [sum([rand(Poisson(μ)) for μ in IFR_ests[k]*KenyaSerology.simple_conv(inc_nai.true_incidence[:,k],p_ID)[1:165]]) for k in 1:length(IFR_ests) ]
mean(T_array .> T_nai)
T = sumsqchange(sample)
T_nai = sum(deaths_nai[1:165])




string(round(fit.posteriormeans[2],sigdigits = 3))*" "*string(round.(fit.CIs[2],sigdigits = 3))
x = (1,(1,2))
string((1,2))

round(π/10000,sigdigits = 3)
