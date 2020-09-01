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

## Create CSV file with fits for each county for each parameter, peak, IFR_crude

df = parameterinferenceovercollection("modelfits/",death_data,0.00264,p_ID)
df.weakpriors = trues(30)

df_sps = parameterinferenceovercollection("modelfits_inferred/",death_data,0.00264,p_ID)
df_sps.weakpriors = falses(17)

df_total = vcat(df,df_sps)


CSV.write("analysis/parameterfitsforeachcounty_vs2.csv",df_total)
