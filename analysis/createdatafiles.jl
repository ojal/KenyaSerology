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


## load fitted priors
@load("analysis/fitted_priors_rural.jld2")
@load("analysis/fitted_priors_semiurban.jld2")


matchedserodata = vcat(sero_data.serodata,zeros(Int64,length(case_data.dates) - size(sero_data.serodata,1),30,2)) #Match length to case data by inflating with zeros
province_matchedserodata = vcat(province_sero_data.serodata,zeros(Int64,length(case_data.dates) - size(province_sero_data.serodata,1),8,2)) #Match length to case data by inflating with zeros
province_sero_data.areas
sero_data.areas
sero_data

nai_prov_data = province_matchedserodata[:,5,:]
sum(nai_prov_data)
sum(matchedserodata[:,7,:])
scatter(province_matchedserodata[:,5,:].-matchedserodata[:,7,:])
scatter!(matchedserodata[:,7,:])


matchedserodata[:,7,:] .= province_matchedserodata[:,5,:]
sero_data.dates = [Date(2020,2,20) + Day(k) for k = 1:168]
sero_data2 = (serodata = matchedserodata,areas = sero_data.areas,dates = [Date(2020,2,20) + Day(k) for k = 1:168])


sero_data = deepcopy(sero_data2)

@save("data/serologydata_21feb_6aug.jld2",sero_data)

function createserodatadataframe(sero_data)
    df_sero = DataFrame(dates = string.(sero_data.dates))
    for (i,county) in enumerate(sero_data.areas)
        str_pos = Symbol(county*"_pos")
        str_neg = Symbol(county*"_neg")
        df_sero[!,str_pos] = sero_data.serodata[:,i,1]
        df_sero[!,str_neg] = sero_data.serodata[:,i,2]
    end
    return df_sero
end
function createPCRpositivedataframe(case_data)
    df_PCRpos = DataFrame(sampledates = string.(case_data.dates))
    for (i,county) in enumerate(case_data.areas)
        str_pos = Symbol(county*"_posswabtest")
        df_PCRpos[!,str_pos] = case_data.cases[:,i]
    end
    return df_PCRpos
end
function createdeathsdataframe(death_data)
    df = DataFrame(dateofdeath = string.(death_data.dates))
    for (i,county) in enumerate(death_data.areas)
        str = Symbol(county*"_deaths")
        df[!,str] = death_data.deaths[:,i]
    end
    return df
end
df_sero = createserodatadataframe(sero_data)
df_PCRpos = createPCRpositivedataframe(case_data)
df_deaths = createdeathsdataframe(death_data)
CSV.write("data/serologicaldata_21feb_6thaug.csv",df_sero)
CSV.write("data/PCRpostiveswabs_21feb_6aug.csv",df_PCRpos)
CSV.write("deaths_21feb_30thsept.csv",df_deaths)


function createPCRpositive_and_negativedataframe(case_data)
    df_PCRpos = DataFrame(sampledates = string.(case_data.dates))
    for (i,county) in enumerate(case_data.areas)
        str_pos = Symbol(county*"_posswabtest")
        df_PCRpos[!,str_pos] = max.(case_data.cases[:,i,1],0)
        str_neg = Symbol(county*"_negswabtest")
        df_PCRpos[!,str_neg] = case_data.cases[:,i,2]
    end
    return df_PCRpos
end
df_PCRposandneg = createPCRpositive_and_negativedataframe(case_data_with_pos_neg)
CSV.write("PCRpostiveandnegswabs_21feb_30thsept.csv",df_PCRposandneg)
