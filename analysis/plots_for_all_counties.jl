## Activate the Project environment, make sure working directory is the KenyaSerology folder
using Pkg;
Pkg.activate(".");
Pkg.precompile();
using Revise, Suppressor

##Use/import relevant packages

using Distributions, Plots, Dates, JLD2, TransformVariables, Optim, FileIO, CSV, DataFrames, OrdinaryDiffEq
using Plots.PlotMeasures
import KenyaSerology
gr()

##Load death data

@load("data/death_data_by_area_21feb_to_17oct.jld2")
@load("data/p_ID.jld2")
## Loop over all fitted counties --- using the fitted Ct rates

fitfiles = readdir("modelfits_pos_neg_fittedCt/", join = true)
oct_first = (Date(2020, 10, 1) - Date(2020, 2, 20)).value

#The saved ODEProblems are from an older version of SciML, suppressed the warning about type reconstruction
#and reconstruct the ODEProblem
for filename in fitfiles
    fit_dict = @suppress_err load(filename)
    model = @suppress_err fit_dict[first(keys(fit_dict))]
    model.prob = KenyaSerology.make_odeproblemforinference(model.contactrate_data;
        startdate = Date(2020, 2, 21),
        enddate = Date(2020, 10, 1))
    name = model.areaname
    uprname = uppercase(name)
    deaths = death_data.deaths[:, death_data.areas.==uprname][1:(end-14)]
    println("Plotting for county $(name)")
    plt_PCR = KenyaSerology.weeklyplotfittoPCRdata(model)
    plt_pop = KenyaSerology.population_plot(model, Val(:monthlyserology))
    plot!(plt_pop, xlims = (0.0, oct_first),
        ylims = (-5, 120),
        guidefontsize = 20,
        tickfontsize = 14,
        legendfontsize = 9,
        legend = :left,
        titlefontsize = 18,
        left_margin = 5mm, right_margin = 10mm,
        yticks = [0, 20, 40, 60, 80, 100])

    plt_deaths = KenyaSerology.plot_deaths(model, deaths, p_ID)

    plot!(plt_deaths, title = "$(name): Observed and predicted deaths")
    savefig(plt_PCR, "countyplot/PCR/PCR_plot_$(name).png")
    savefig(plt_pop, "countyplot/population_exposure/pop_exposure_plot_$(name).png")
    savefig(plt_deaths, "countyplot/deaths/deaths_plot_$(name).png")

end

##

filename = fitfiles[30]
fit_dict = @suppress_err load(filename)
model = @suppress_err fit_dict[first(keys(fit_dict))]

model.prob = KenyaSerology.make_odeproblemforinference(model.contactrate_data;
    startdate = Date(2020, 2, 21),
    enddate = Date(2020, 10, 1))
name = model.areaname
plt_pop = KenyaSerology.population_plot(model, Val(:monthlyserology))
plot!(plt_pop, xlims = (0.0, oct_first),
    ylims = (-5, 120),
    guidefontsize = 20,
    tickfontsize = 14,
    legendfontsize = 9,
    legend = :left,
    titlefontsize = 18,
    left_margin = 5mm, right_margin = 10mm,
    yticks = [0, 20, 40, 60, 80, 100])

uprname = uppercase(name)
deaths = death_data.deaths[:, death_data.areas.==uprname][1:(end-14)]
plt_deaths = KenyaSerology.plot_deaths(model, deaths, p_ID)
