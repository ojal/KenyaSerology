push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaSerologyPrivate/src"))

using Distributions,Plots,Dates,JLD2,TransformVariables,Optim,FileIO,DataFrames,CSV,MAT,StatsPlots
import KenyaSerology
# include("analysis_functions.jl");
# include("getdata.jl");
@load("data/case_data_with_pos_neg_21feb_to_30sept.jld2")
@load("data/death_data_by_area_21feb_to_17oct.jld2")



@load("analysis/collecteddata_lateseptbayesianfit_posneg_opt.jld2")


kenya_case_data = vec(sum(case_data_with_pos_neg.cases[:,:,1],dims = 2))
weekly_kenya_case_data = KenyaSerology.group_PCR_pos_by_week(kenya_case_data)
scatter(weekly_kenya_case_data)

kenya_death_data = vec(sum(death_data.deaths[death_data.dates .<= Date(2020,9,30),:],dims = 2))
weekly_kenya_death_data = KenyaSerology.group_PCR_pos_by_week(kenya_death_data)


kenya_dates = [Date(2020,2,20) + Day(k) for k = 1:224]
weekdates_on_scale = collect(0:(length(weekly_kenya_death_data)-1)).*7 .- 4
scatter(weekdates_on_scale,weekly_kenya_death_data)


function combinedata(combindeddata,case_data_with_pos_neg)
        kenyaincidence = zeros(size(combindeddata[1].mean_pred_incidence))
        kenyaincidence_std = zeros(size(combindeddata[1].mean_pred_incidence))
        kenyaPCR = zeros(size(combindeddata[1].mean_pred_PCR_NB))
        kenyaPCR_std = zeros(size(combindeddata[1].mean_pred_PCR_NB))
        kenyadeaths = zeros(size(combindeddata[1].pred_deaths))
        kenyahosps = zeros(size(combindeddata[1].pred_hosps))

        for data in combindeddata
			days_with_neg_tests = vec(case_data_with_pos_neg.cases[:,case_data_with_pos_neg.areas .== uppercase(data.area),2]) .>= 0
			total_tests = vec(case_data_with_pos_neg.cases[:,case_data_with_pos_neg.areas .== uppercase(data.area),1] .+ case_data_with_pos_neg.cases[:,case_data_with_pos_neg.areas .== uppercase(data.area),2])

            kenyaincidence .+= data.mean_pred_incidence
            kenyaincidence_std .+= data.std_pred_incidence.^2
			for t = 1:size(case_data_with_pos_neg.cases,1)
				if days_with_neg_tests[t]
					kenyaPCR[t] += data.mean_pred_P_PCR_BB[t]*total_tests[t]
					kenyaPCR_std[t] += (data.std_pred_P_PCR_BB[t]*total_tests[t])^2
				else
					kenyaPCR[t] += data.mean_pred_PCR_NB[t]
					kenyaPCR_std[t] += (data.std_pred_PCR_NB[t])^2
				end
			end
            kenyadeaths .+= data.pred_deaths
            kenyahosps .+= data.pred_hosps

        end

        kenyaincidence_std = sqrt.(kenyaincidence_std)
        kenyaPCR_std = sqrt.(kenyaPCR_std)
        kenyadeaths_upred = [invlogcdf(Poisson(d),log(0.975)) for d in kenyadeaths]
        kenyadeaths_lpred = [invlogcdf(Poisson(d),log(0.025)) for d in kenyadeaths]
        kenyahosps_upred = [invlogcdf(Poisson(d),log(0.975)) for d in kenyahosps]
        kenyahosps_lpred = [invlogcdf(Poisson(d),log(0.025)) for d in kenyahosps]

        return (kenyaincidence=kenyaincidence,
                kenyaincidence_std=kenyaincidence_std,
                kenyaPCR = kenyaPCR,
                kenyaPCR_std = kenyaPCR_std,
                kenyadeaths = kenyadeaths,
                kenyadeaths_upred = kenyadeaths_upred .- kenyadeaths,
                kenyadeaths_lpred = kenyadeaths .- kenyadeaths_lpred,
                kenyahosps = kenyahosps,
                kenyahosps_upred = kenyahosps_upred .- kenyahosps,
                kenyahosps_lpred = kenyahosps .- kenyahosps_lpred)
end
sept_30_day = (Date(2020,9,30) - Date(2020,2,20)).value
# kenyadata = combinedata(vcat(collecteddata_fitted,collecteddata_data_inferred))
kenyadata = combinedata(collecteddata_latesept_opt,case_data_with_pos_neg)

plot(kenyadata.kenyaincidence[1:223])
## Plot Kenya data
# xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:10 ]
# xticklabs = [(monthname(k)[1:3]) for k = 3:12]

xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,20)).value for k = 1:8 ]
xticklabs2020 = [(monthname(k)[1:3]) for k = 3:10]
# xticklabs2021 = [(monthname(k)[1:3])*"-21" for k = 1:10]
# xticklabs = vcat(xticklabs2020,xticklabs2021)
xticklabs = xticklabs2020
pyplot()
# gr()
weekly_kenyaPCR_pred = KenyaSerology.group_PCR_pos_by_week(kenyadata.kenyaPCR[1:223])
weekly_kenyaPCR_std = KenyaSerology.group_PCR_pos_by_week(kenyadata.kenyaPCR_std[1:223])

plt_kenya_cases = plot(weekdates_on_scale[1:(end-1)],weekly_kenyaPCR_pred[1:(end-1)],
        ribbon = (min.(3*weekly_kenyaPCR_std[1:(end-1)],weekly_kenyaPCR_pred[1:(end-1)]),3*weekly_kenyaPCR_std[1:(end-1)]),
        lw =3,
        color = :red,
        fillalpha = 0.3,
        lab = "Predicted trend")
scatter!(plt_kenya_cases,weekdates_on_scale[1:(end-1)],weekly_kenya_case_data[1:(end-1)],ms = 5,
                        color=:blue,
                        legend=:topleft,
                        label = "Positive swab tests",
                        xticks = (xticktimes,xticklabs),
                        ylabel = "Weekly incidence",
						xlims = (0,234),
                        size = (700,500),dpi = 250,
                        title = "Kenya PCR positive test prediction and data",
                        guidefont = 20,
                        tickfont = 14,
                        legendfont = 12,
                        titlefont = 14)


savefig(plt_kenya_cases,"plotsforpaper/revisedkenyan_swab_test_plot_fitted_up_to_30th_sept.pdf")

##
@load("data/case_data_with_symptoms_by_area_21feb_to_5oct.jld2")
kenya_symptomatic_cases = vec(sum(case_data_with_symptoms.cases[:,:,1], dims = 2))
kenya_asymptomatic_cases = vec(sum(case_data_with_symptoms.cases[:,:,2], dims = 2))

plt_kenya_cases_symptoms = scatter(kenya_symptomatic_cases.+1,ms = 4,yscale = :log10,
                        color=:red,
                        legend=:topright,
                        label = "Symptomatic cases",
                        xticks = (xticktimes,xticklabs),
                        ylabel = "Daily incidence",
						xlims = (0,340),
                        size = (700,500),dpi = 250,
                        title = "",
                        guidefont = 20,
                        tickfont = 8,
                        legendfont = 12,
                        titlefont = 14)
scatter!(plt_kenya_cases_symptoms,kenya_asymptomatic_cases.+1,ms=4,
 		lab = "Asymptomatic cases",
		color = :blue)


## Plot Deaths

weekly_kenyadeath_pred = KenyaSerology.group_PCR_pos_by_week(kenyadata.kenyadeaths[1:223])
weekly_kenyaPCR_upred = [quantile(Poisson(d),0.975) - d  for d in weekly_kenyadeath_pred]
weekly_kenyaPCR_lpred = [d - quantile(Poisson(d),0.025)  for d in weekly_kenyadeath_pred]
plot(weekly_kenyaPCR_lpred)
plot!(weekly_kenyadeath_pred)
plot!(weekly_kenyaPCR_upred)
plt_kenya_deaths = plot(weekdates_on_scale[1:(end-1)],weekly_kenyadeath_pred[1:(end-1)],
        ribbon = (weekly_kenyaPCR_lpred[1:(end-1)],weekly_kenyaPCR_upred[1:(end-1)]),
        lw =3,
        color = :black,
        fillalpha = 0.2,
        lab = "Predicted trend")
scatter!(plt_kenya_deaths,weekdates_on_scale[1:(end-1)],weekly_kenya_death_data[1:(end-1)],
                        ms = 5,
                        color=:black,
                        legend=:topleft,
                        label = "Reported deaths",
                        xticks = (xticktimes,xticklabs),
                        ylabel = "Weekly deaths",
                        size = (700,500),dpi = 250,
                        title = "Kenya daily deaths prediction and data",
                        guidefont = 20,
                        tickfont = 12,
                        legendfont = 12,
                        titlefont = 18,
						xlims = (0,225))


plot!(weekdates_on_scale[1:(end-1)],cumsum(weekly_kenyadeath_pred[1:(end-1)]),lw = 1,
        grid = nothing,
        inset = (1,bbox(-0.25, -0, 0.3, 0.3, :center)),
        xticks = (xticktimes[1:2:end],xticklabs[1:2:end]),
        yticks = 0:200:1200,
        ylims = (0.,1200.),
        lab="",
        color = :black,
        subplot = 2,
        titlefont = 10,
        title = "Cumulative observed deaths",
        bg_inside = nothing)
scatter!(weekdates_on_scale[1:(end-1)],cumsum(weekly_kenya_death_data[1:(end-1)]),
        ms = 5,
        color = :black,
        lab = "",
        subplot = 2)

savefig(plt_kenya_deaths,"plotsforpaper/revisedkenyan_deaths_up_to_30sept.pdf")

plt_kenya_incidence = plot(kenyadata.kenyaincidence.+1,yscale = :log10,
        ribbon = (min.(3*kenyadata.kenyaincidence_std,kenyadata.kenyaincidence),3*kenyadata.kenyaincidence_std),
        yticks = ([1,2,11,101,1001,10001,100001],[0,1,10,100,1000,10000,100000]),
        xticks = (xticktimes,xticklabs),
        lw =3,
        lab = "",
        color = :red,
        fillalpha = 0.3,
        title = "Predicted Kenyan true daily incidence",
        titlefont = 18,
        guidefont = 18,
        tickfont = 6,
        legendfont = 8,
        size = (700,500),dpi=250,
        ylabel = "Daily incidence")
daylessthan10000 = findfirst((kenyadata.kenyaincidence)[91:end] .<= 10000)+90
daylessthan100 = findfirst((kenyadata.kenyaincidence)[91:end] .<= 100)+90

plot!(plt_kenya_incidence,[daylessthan100,daylessthan100],[1,1000000],lw=2,
	ls = :dash,
	lab = "Less than 100 daily new infections on $(Date(2020,2,20)+Day(daylessthan100))",
	legend = :topright)
plot!(plt_kenya_incidence,[daylessthan10000,daylessthan10000],[1,1000000],lw=2,
	ls = :dash,
	lab = "Less than 10,000 daily new infections on $(Date(2020,2,20)+Day(daylessthan10000))",
	legend = :topright)
savefig(plt_kenya_incidence,"plotsforreport/kenyan_incidence_plotlateaugfit.pdf")
