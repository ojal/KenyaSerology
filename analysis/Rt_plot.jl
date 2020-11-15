#Get instantaneous and effective R(t) for each county

push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaSerologyPrivate/src"))

using Distributions,Plots,Dates,JLD2,TransformVariables,Optim,FileIO,CSV,DataFrames
using Revise
import KenyaSerology
@load("data/projected_contact_data_09102020.jld2")
@load("modelfits_pos_neg_opt/Nairobi_model_opt.jld2")

@load("analysis/collecteddata_lateseptbayesianfit_posneg_opt.jld2")
data = collecteddata_latesept_opt[1]

dir_name = "modelfits_pos_neg_opt"
model_paths = [ joinpath(dir_name,filename) for filename in readdir(dir_name)]

nai_path = "modelfits_pos_neg_opt/Nairobi_model_opt.jld2"
mom_path = "modelfits_pos_neg_opt/Mombasa_model_opt.jld2"
other_paths = setdiff(model_paths,[nai_path,mom_path])

function gather_smoothed_ct_fits(paths,projected_contactrate)
    gathered_smoothed_ct = []
    for path in paths
        model = load(path)["model_opt"]
        smoothed_ct = KenyaSerology.optimise_ct_model(model,projected_contactrate)
        push!(gathered_smoothed_ct,(model.areaname,smoothed_ct))
    end
    return gathered_smoothed_ct
end

nai_smoothed_ct = gather_smoothed_ct_fits([nai_path],projected_contactrate_nairobi)

mom_smoothed_ct = gather_smoothed_ct_fits([mom_path],projected_contactrate_mombasa)

other_smoothed_ct_1 = gather_smoothed_ct_fits(other_paths[1:22],projected_contactrate_kenya)

other_smoothed_ct_2 = gather_smoothed_ct_fits(other_paths[23:34],projected_contactrate_kenya)
@save("analysis/smoothed_ct_2.jld2",other_smoothed_ct_2)


@load("analysis/gathered_ct_fits.jld2")
@load("modelfits_pos_neg_opt/Nairobi_model_opt.jld2")
projected_contactrate_kenya.area

function add_contact_rate_and_calculate_Rt(path,gathered_ct_data)
    model_opt = load(path)["model_opt"]
    f = findfirst([data[1] == model_opt.areaname for data in gathered_ct_data])
    ct = gathered_ct_data[f][2]
    fitted_contact_data = (date = [Date(2020,2,20) + Day(k) for k = 1:length(ct)],
                    area = model_opt.contactrate_data.area,
                    contactrate = ct)
    model_opt.contactrate_data = fitted_contact_data
    @save("modelfits_pos_neg_opt_for_paper/$(model_opt.areaname)_model_opt.jld2",model_opt)
    mean_R = mean(model_opt.MCMC_results.chain[:,1,1])
    R_lb = mean_R - quantile(model_opt.MCMC_results.chain[:,1,1],0.025)
    R_ub =  quantile(model_opt.MCMC_results.chain[:,1,1],0.975) - mean_R

    return (mean_Rt = mean_R.*ct,Rt_lpred = R_lb.*ct, Rt_upred = R_ub.*ct,name = model_opt.areaname)
end


gathered_Rt_data = [add_contact_rate_and_calculate_Rt(path,gathered_smoothed_cts) for path in model_paths]
@save("analysis/gathered_Rt_data.jld2",gathered_Rt_data)

## Plot Rt and eff Rt for Nai and mom
@load("analysis/gathered_Rt_data.jld2",gathered_Rt_data)

f_notnaiormom = findall([data.name != "Nairobi" && data.name != "Mombasa" for data in gathered_Rt_data])
Rt_data_outsidenaiandmom = gathered_Rt_data[f_notnaiormom]

function IQR_Rt(Rt_data)
    Rt_array = zeros(length(Rt_data[1].mean_Rt),length(Rt_data))
    for (j,data) in enumerate(Rt_data)
         Rt_array[:,j] = data.mean_Rt
    end
    R_med = [median(Rt_array[t,:]) for t = 1:size(Rt_array,1)]
    R_LQ = [quantile(Rt_array[t,:],0.25) for t = 1:size(Rt_array,1)]
    R_UQ = [quantile(Rt_array[t,:],0.75) for t = 1:size(Rt_array,1)]
    R_mean = [mean(Rt_array[t,:]) for t = 1:size(Rt_array,1)]

    return (R_med=R_med,
            R_mean=R_mean,
            R_LQ=R_LQ,
            R_UQ=R_UQ)
end


Rt_summary_outside_naiandmom = IQR_Rt(Rt_data_outsidenaiandmom)

# Get the eff Rt for Nairobi and Mombasa
@load("modelfits_pos_neg_opt_for_paper/Mombasa_model_opt.jld2")
model_mombasa = deepcopy(model_opt)
@load("modelfits_pos_neg_opt_for_paper/Nairobi_model_opt.jld2")
model_nairobi = deepcopy(model_opt)

inc_nai = KenyaSerology.incidence_across_samples(model_nairobi,224)
inc_mom = KenyaSerology.incidence_across_samples(model_mombasa,224)

eff_trans_nai = mean(inc_nai.eff_transmission[1:220],dims = 2)
effRt_nai = eff_trans_nai.*model_nairobi.contactrate_data.contactrate

eff_trans_mom = mean(inc_mom.eff_transmission[1:220],dims = 2)
effRt_mom = eff_trans_mom.*model_mombasa.contactrate_data.contactrate

plot(effRt_nai)
plot!(effRt_mom)

f_nai = findfirst([data.name == "Nairobi" for data in gathered_Rt_data])
f_mom = findfirst([data.name == "Mombasa" for data in gathered_Rt_data])
mean_R_nai = gathered_Rt_data[f_nai].mean_Rt
mean_R_mom = gathered_Rt_data[f_mom].mean_Rt
date_for_plot = [Date(2020,2,20) + Day(k) for k = 1:length(mean_R_nai)]


xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:8 ]
xticklabs = [(monthname(k)[1:3]) for k = 3:10]


plt_Rt = plot(mean_R_nai,
    lw = 2.5,
    color = :red,
    ylims = (0.5,3.55),
    xticks = (xticktimes,xticklabs),
    lab = "Nairobi: R(t)",
    size = (700,500),dpi = 250,
    guidefont = 20,
    tickfont = 12,
    legendfont = 12,
    titlefont = 18,
    title = "Inferred Kenyan Reproductive Ratios")
plot!(plt_Rt,mean_R_mom,
    lw = 2.5,
    color = :green,
    lab = "Mombasa: R(t)")
plot!(plt_Rt,Rt_summary_outside_naiandmom.R_mean,lab = "Other counties: R(t) IQR",
    ribbon = (Rt_summary_outside_naiandmom.R_mean .-Rt_summary_outside_naiandmom.R_LQ ,
                Rt_summary_outside_naiandmom.R_UQ .- Rt_summary_outside_naiandmom.R_mean),
    color = :blue,
    lw = 2.5,
    fillalpha = 0.2)
plot!(plt_Rt,effRt_nai,
    lw = 2.5,ls = :dot,
    color = :red,
    lab = "Nairobi: Eff. R(t)")
plot!(plt_Rt,effRt_mom,
    lw = 2.5,ls = :dot,
    color = :green,
    lab = "Mombasa: Eff. R(t)")

#Plot dates
#group 1
firstcase = (Date(2020,3,12)-Date(2020,2,20)).value
suspendpublicgathering = (Date(2020,3,13)-Date(2020,2,20)).value
schoolsclosed = (Date(2020,3,15)-Date(2020,2,20)).value
retroactivequarantine = (Date(2020,3,17)-Date(2020,2,20)).value

#group 2
inboundflightsuspended = (Date(2020,3,25)-Date(2020,2,20)).value
nationalcurfew = (Date(2020,3,27)-Date(2020,2,20)).value
# testingquarantined = (Date(2020,3,27)-Date(2020,2,20)).value
regionallockdowns = (Date(2020,4,6)-Date(2020,2,20)).value

#group 3
localnaiandmomlockdowns = (Date(2020,5,6)-Date(2020,2,20)).value
tazandsomborderclosed = (Date(2020,5,16)-Date(2020,2,20)).value

plot!(plt_Rt,[firstcase,firstcase],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "Restrictions")
plot!(plt_Rt,[suspendpublicgathering,suspendpublicgathering],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")
plot!(plt_Rt,[schoolsclosed,schoolsclosed],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")
plot!(plt_Rt,[retroactivequarantine,retroactivequarantine],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")
plot!(plt_Rt,[inboundflightsuspended,inboundflightsuspended],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")
plot!(plt_Rt,[nationalcurfew,nationalcurfew],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")
plot!(plt_Rt,[regionallockdowns,regionallockdowns],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")
plot!(plt_Rt,[localnaiandmomlockdowns,localnaiandmomlockdowns],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")
plot!(plt_Rt,[tazandsomborderclosed,tazandsomborderclosed],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")

ann1 = (19, 0.76, Plots.text("1", :left))
ann2 = (39, 0.76, Plots.text("2", :left))
ann3 = (79, 0.76, Plots.text("3", :left))
plot!(plt_Rt,annotations = [ann1,ann2,ann3])
plot!(plt_Rt,ylabel = "R(t)")
savefig(plt_Rt,"plotsforpaper/revised_Rt_plot.pdf")
