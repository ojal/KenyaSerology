push!(LOAD_PATH, joinpath(homedir(),"GitHub/Kenya-Serology/src"))

using Distributions,Plots,Dates,JLD2,TransformVariables,Optim,FileIO,DataFrames,CSV,MAT,StatsPlots
import KenyaSerology
include("analysis_functions.jl");
include("getdata.jl");

worldometerdeaths =[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,2,1,0,0,2,0,0,1,0,0,1,1,0,1,1,0,1,2,0,0,0,0,0,0,0,0,0,1,2,4,1,2,0,0,2,3,0,1,2,1,3,4,2,
        3,5,0,0,0,0,0,0,0,1,1,0,3,3,4,1,1,5,2,3,4,1,4,1,1,3,1,3,4,4,3,1,1,2,10,2,2,2,2,3,2,2,5,4,2,1,4,1,3,2,5,1,4,
        3,2,4,8,3,1,12,5,7,8,5,3,9,4,12,10,3,11,4,2,5,14,12,14,16,23,5,13,6,3,8,14,5,2,3,15,18]
worldometerdeaths[7:end]

@load("analysis/collecteddata_fitted.jld2")
@load("analysis/collecteddata_data_inferred.jld2")
@load("data/case_data_by_area_21feb_to_6aug.jld2")
@load("data/death_data_by_area_21feb_to_6aug.jld2")
kenya_case_data = vec(sum(case_data.cases,dims = 2))
kenya_death_data = vec(sum(death_data.deaths,dims = 2))

kenyadata = combinedata(vcat(collecteddata_fitted,collecteddata_data_inferred))

## Plot Kenya data
xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:10 ]
xticklabs = [(monthname(k)[1:3]) for k = 3:12]
pyplot()

plt_kenya_cases = scatter(kenya_case_data,ms = 7,
                        color=:blue,
                        legend=:topleft,
                        label = "Positive swab tests",
                        xticks = (xticktimes,xticklabs),
                        ylabel = "Daily incidence",
                        size = (700,500),dpi = 250,
                        title = "Kenya PCR positive test prediction and data",
                        guidefont = 20,
                        tickfont = 14,
                        legendfont = 12,
                        titlefont = 18)
plot!(plt_kenya_cases,kenyadata.kenyaPCR,
        ribbon = (min.(3*kenyadata.kenyaPCR_std,kenyadata.kenyaPCR),3*kenyadata.kenyaPCR_std),
        lw =3,
        color = :red,
        fillalpha = 0.3,
        lab = "Predicted trend")
savefig(plt_kenya_cases,"plotsforpaper/kenyan_swab_test_plot.pdf")

plt_kenya_deaths = scatter(kenya_death_data,
                        ms = 10,
                        color=:black,
                        legend=:topleft,
                        label = "Reported deaths",
                        xticks = (xticktimes,xticklabs),
                        ylabel = "Daily deaths",
                        size = (700,500),dpi = 250,
                        title = "Kenya daily deaths prediction and data",
                        guidefont = 20,
                        tickfont = 14,
                        legendfont = 12,
                        titlefont = 18)
plot!(plt_kenya_deaths,kenyadata.kenyadeaths,
        ribbon = (kenyadata.kenyadeaths_lpred,kenyadata.kenyadeaths_upred),
        lw =3,
        color = :black,
        fillalpha = 0.2,
        lab = "Predicted trend")

plot!(cumsum(kenyadata.kenyadeaths),lw = 1,
        grid = nothing,
        inset = (1,bbox(0.25, -0.2, 0.3, 0.3, :center)),
        xticks = (xticktimes[1:5:end],xticklabs[1:5:end]),
        yticks = 0:200:1000,
        ylims = (0.,1000.),
        lab="",
        color = :black,
        subplot = 2,
        titlefont = 10,
        title = "Cumulative observed deaths",
        bg_inside = nothing)
scatter!(cumsum(kenya_death_data),
        ms = 5,
        color = :black,
        lab = "",
        subplot = 2)

savefig(plt_kenya_deaths,"plotsforpaper/kenyan_deaths_plot.pdf")

plt_kenya_incidence = plot(kenyadata.kenyaincidence.+1,yscale = :log10,
        ribbon = (min.(3*kenyadata.kenyaincidence_std,kenyadata.kenyaincidence),3*kenyadata.kenyaincidence_std),
        yticks = ([1,2,11,101,1001,10001,100001],[0,1,10,100,1000,10000,100000]),
        xticks = (xticktimes,xticklabs),
        lw =3,
        lab = "",
        color = :red,
        fillalpha = 0.3,
        title = "Predicted Kenyan true daily incidence",
        titlefont = (18,"helvetica"),
        guidefont = (18,"helvetica"),
        tickfont = (10,"helvetica"),
        legendfont = (12,"helvetica"),
        size = (700,500),
        ylabel = "Daily incidence")
savefig(plt_kenya_incidence,"plotsforpaper/kenyan_incidence_plot.pdf")
