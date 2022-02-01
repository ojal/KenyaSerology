#Get instantaneous and effective R(t) for each county

## Activate the Project environment, make sure working directory is the KenyaSerology folder
using Pkg;
Pkg.activate(".");
Pkg.precompile();
using Revise, Suppressor

##
using Distributions, Plots, Dates, JLD2, TransformVariables, Optim, FileIO, CSV, DataFrames
using Revise
import KenyaSerology

##

## Load the gathered Rt data, keeping Nairobi and Mombasa separate
@load("analysis/gathered_Rt_data.jld2", gathered_Rt_data)

f_notnaiormom = findall([data.name != "Nairobi" && data.name != "Mombasa" for data in gathered_Rt_data])
Rt_data_outsidenaiandmom = gathered_Rt_data[f_notnaiormom]

"""
        function IQR_Rt(Rt_data)

Returns the IQR over the Rt estimates                
"""
function IQR_Rt(Rt_data)
        Rt_array = zeros(length(Rt_data[1].mean_Rt), length(Rt_data))
        for (j, data) in enumerate(Rt_data)
                Rt_array[:, j] = data.mean_Rt
        end
        R_med = [median(Rt_array[t, :]) for t = 1:size(Rt_array, 1)]
        R_LQ = [quantile(Rt_array[t, :], 0.25) for t = 1:size(Rt_array, 1)]
        R_UQ = [quantile(Rt_array[t, :], 0.75) for t = 1:size(Rt_array, 1)]
        R_mean = [mean(Rt_array[t, :]) for t = 1:size(Rt_array, 1)]

        return (R_med = R_med,
                R_mean = R_mean,
                R_LQ = R_LQ,
                R_UQ = R_UQ)
end
Rt_summary_outside_naiandmom = IQR_Rt(Rt_data_outsidenaiandmom)

## Get the eff Rt for Nairobi and Mombasa

#The saved ODEProblems are from an older version of SciML, suppressed the warning about type reconstruction
#and reconstruct the ODEProblem
@load("modelfits_pos_neg_fittedCt/Mombasa_model_opt.jld2")
model_mombasa = deepcopy(model_opt)
model_mombasa.prob = KenyaSerology.make_odeproblemforinference(model_mombasa.contactrate_data;
        startdate = Date(2020, 2, 21),
        enddate = Date(2020, 10, 1))
@load("modelfits_pos_neg_fittedCt/Nairobi_model_opt.jld2")
model_nairobi = deepcopy(model_opt)
model_nairobi.prob = KenyaSerology.make_odeproblemforinference(model_nairobi.contactrate_data;
        startdate = Date(2020, 2, 21),
        enddate = Date(2020, 10, 1))


inc_nai = KenyaSerology.incidence_across_samples(model_nairobi, 224)
inc_mom = KenyaSerology.incidence_across_samples(model_mombasa, 224)

## Calculate Eff. Rt averaged over MCMC ensemble
eff_trans_nai = mean(inc_nai.eff_transmission[1:220], dims = 2)
effRt_nai = eff_trans_nai .* model_nairobi.contactrate_data.contactrate

eff_trans_mom = mean(inc_mom.eff_transmission[1:220], dims = 2)
effRt_mom = eff_trans_mom .* model_mombasa.contactrate_data.contactrate

##Plot Rt and eff. Rt

f_nai = findfirst([data.name == "Nairobi" for data in gathered_Rt_data])
f_mom = findfirst([data.name == "Mombasa" for data in gathered_Rt_data])
mean_R_nai = gathered_Rt_data[f_nai].mean_Rt
mean_R_mom = gathered_Rt_data[f_mom].mean_Rt
date_for_plot = [Date(2020, 2, 20) + Day(k) for k = 1:length(mean_R_nai)]


xticktimes = [((Date(2020, 2, 1) + Month(k)) - Date(2020, 2, 21)).value for k = 1:8]
xticklabs = [(monthname(k)[1:3]) for k = 3:10]

## Dates of restrictions
#Plot dates of 
#group 1
firstcase = (Date(2020, 3, 12) - Date(2020, 2, 20)).value
suspendpublicgathering = (Date(2020, 3, 13) - Date(2020, 2, 20)).value
schoolsclosed = (Date(2020, 3, 15) - Date(2020, 2, 20)).value
retroactivequarantine = (Date(2020, 3, 17) - Date(2020, 2, 20)).value

#group 2
inboundflightsuspended = (Date(2020, 3, 25) - Date(2020, 2, 20)).value
nationalcurfew = (Date(2020, 3, 27) - Date(2020, 2, 20)).value
# testingquarantined = (Date(2020,3,27)-Date(2020,2,20)).value
regionallockdowns = (Date(2020, 4, 6) - Date(2020, 2, 20)).value

#group 3
localnaiandmomlockdowns = (Date(2020, 5, 6) - Date(2020, 2, 20)).value
tazandsomborderclosed = (Date(2020, 5, 16) - Date(2020, 2, 20)).value

#Easing restrictions
endnaiandmomlockdowns = (Date(2020, 7, 6) - Date(2020, 2, 20)).value
resumeairtravel = (Date(2020, 8, 1) - Date(2020, 2, 20)).value


##
gr()
plt_Rt = plot(mean_R_nai,
        lw = 2.5,
        color = :red,
        ylims = (0.5, 3.55),
        xticks = (xticktimes, xticklabs),
        lab = "Nairobi: R(t)",
        size = (700, 500), dpi = 250,
        guidefontsize = 20,
        tickfontsize = 12,
        legendfontsize = 12,
        titlefontsize = 18,
        title = "Inferred Kenyan Reproductive Ratios")

plot!(plt_Rt, mean_R_mom,
        lw = 2.5,
        color = :green,
        lab = "Mombasa: R(t)")
plot!(plt_Rt, Rt_summary_outside_naiandmom.R_mean, lab = "Other counties: R(t) IQR",
        ribbon = (Rt_summary_outside_naiandmom.R_mean .- Rt_summary_outside_naiandmom.R_LQ,
                Rt_summary_outside_naiandmom.R_UQ .- Rt_summary_outside_naiandmom.R_mean),
        color = :blue,
        lw = 2.5,
        fillalpha = 0.2)
plot!(plt_Rt, effRt_nai,
        lw = 2.5, ls = :dot,
        color = :red,
        lab = "Nairobi: Eff. R(t)")
plot!(plt_Rt, effRt_mom,
        lw = 2.5, ls = :dot,
        color = :green,
        lab = "Mombasa: Eff. R(t)")



plot!(plt_Rt, [firstcase, firstcase], [0.5, 3.5],
        ls = :dash, lw = 0.5, color = :black, lab = "Restrictions")
plot!(plt_Rt, [suspendpublicgathering, suspendpublicgathering], [0.5, 3.5],
        ls = :dash, lw = 0.5, color = :black, lab = "")
plot!(plt_Rt, [schoolsclosed, schoolsclosed], [0.5, 3.5],
        ls = :dash, lw = 0.5, color = :black, lab = "")
plot!(plt_Rt, [retroactivequarantine, retroactivequarantine], [0.5, 3.5],
        ls = :dash, lw = 0.5, color = :black, lab = "")
plot!(plt_Rt, [inboundflightsuspended, inboundflightsuspended], [0.5, 3.5],
        ls = :dash, lw = 0.5, color = :black, lab = "")
plot!(plt_Rt, [nationalcurfew, nationalcurfew], [0.5, 3.5],
        ls = :dash, lw = 0.5, color = :black, lab = "")
plot!(plt_Rt, [regionallockdowns, regionallockdowns], [0.5, 3.5],
        ls = :dash, lw = 0.5, color = :black, lab = "")
plot!(plt_Rt, [localnaiandmomlockdowns, localnaiandmomlockdowns], [0.5, 3.5],
        ls = :dash, lw = 0.5, color = :black, lab = "")
plot!(plt_Rt, [tazandsomborderclosed, tazandsomborderclosed], [0.5, 3.5],
        ls = :dash, lw = 0.5, color = :black, lab = "")
vline!(plt_Rt, [endnaiandmomlockdowns], ls = :dash, lw = 2.5, color = 2, lab = "End movement restriction to/from main cities")
vline!(plt_Rt, [resumeairtravel], ls = :dash, lw = 2.5, color = 3, lab = "Resumed international air travel")

ann1 = (19, 0.76, Plots.text("1", :left))
ann2 = (39, 0.76, Plots.text("2", :left))
ann3 = (79, 0.76, Plots.text("3", :left))
plot!(plt_Rt, annotations = [ann1, ann2, ann3])
plot!(plt_Rt, ylabel = "R(t)")
savefig(plt_Rt, "Rt_plot.pdf")
