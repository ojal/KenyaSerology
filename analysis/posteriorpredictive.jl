push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaSerology/src"))

using Distributions,Plots,Dates,JLD2,TransformVariables,Optim,FileIO,DataFrames,CSV,MAT,StatsPlots
import KenyaSerology
include("analysis_functions.jl");
include("getdata.jl");

@load("modelfits/nairobi_model_Peff_shortsi_var_testing.jld2")
deaths_nai = death_data.deaths[:,death_data.areas.=="NAIROBI"][1:165]

inc_nai = KenyaSerology.incidence_across_samples(nairobi_model_Peff_shortsi_var_testing,315)
D = sum(deaths_nai)
IFR_array = [Erlang(D+1,1/((1/0.00264) + sum(KenyaSerology.simple_conv(inc_nai.true_incidence[:,n],p_ID)[1:165]))) for n = 1: 10000]


@time lpd_actual = death_mean_pred(deaths_nai,inc_nai.true_incidence,mean.(IFR_array),p_ID)



lpd_array_nai = generate_simulated_death_lpds(nairobi_model_Peff_shortsi_var_testing,deaths_nai,p_ID)
histogram(lpd_array_nai,norm = :pdf,bins = 50)
plot!([lpd_actual,lpd_actual],[0.,0.02])
sum(lpd_actual .< lpd_array_nai)/1000
