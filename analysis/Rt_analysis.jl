#Get instantaneous and effective R(t) for each county

## Activate the Project environment, make sure working directory is the KenyaSerology folder
using Pkg;
Pkg.activate(".");
Pkg.precompile();
using Revise, Suppressor

## Load packages
using Distributions, Plots, Dates, JLD2, TransformVariables, Optim, FileIO, CSV, DataFrames
using Revise
import KenyaSerology

## Load data and gather fits

@load("data/projected_contact_data_09102020.jld2")
@load("modelfits_pos_neg_fittedCt/Nairobi_model_opt.jld2")
@load("analysis/collecteddata_lateseptbayesianfit_posneg_opt.jld2")
data = collecteddata_latesept_opt[1]

dir_name = "modelfits_pos_neg_fittedCt"
model_paths = [joinpath(dir_name, filename) for filename in readdir(dir_name)]

nai_path = "modelfits_pos_neg_fittedCt/Nairobi_model_opt.jld2"
mom_path = "modelfits_pos_neg_fittedCt/Mombasa_model_opt.jld2"
other_paths = setdiff(model_paths, [nai_path, mom_path])

function gather_smoothed_ct_fits(paths, projected_contactrate)
        gathered_smoothed_ct = []
        for path in paths
                model = load(path)["model_opt"]
                smoothed_ct = KenyaSerology.optimise_ct_model(model, projected_contactrate)
                push!(gathered_smoothed_ct, (model.areaname, smoothed_ct))
        end
        return gathered_smoothed_ct
end

nai_smoothed_ct = gather_smoothed_ct_fits([nai_path], projected_contactrate_nairobi)
mom_smoothed_ct = gather_smoothed_ct_fits([mom_path], projected_contactrate_mombasa)
other_smoothed_ct_1 = gather_smoothed_ct_fits(other_paths[1:22], projected_contactrate_kenya)
other_smoothed_ct_2 = gather_smoothed_ct_fits(other_paths[23:34], projected_contactrate_kenya)
@save("analysis/smoothed_ct_2.jld2", other_smoothed_ct_2)


@load("analysis/gathered_ct_fits.jld2")
@load("modelfits_pos_neg_opt/Nairobi_model_opt.jld2")
projected_contactrate_kenya.area

function add_contact_rate_and_calculate_Rt(path, gathered_ct_data)
        model_opt = load(path)["model_opt"]
        f = findfirst([data[1] == model_opt.areaname for data in gathered_ct_data])
        ct = gathered_ct_data[f][2]
        fitted_contact_data = (date = [Date(2020, 2, 20) + Day(k) for k = 1:length(ct)],
                area = model_opt.contactrate_data.area,
                contactrate = ct)
        model_opt.contactrate_data = fitted_contact_data
        # @save("modelfits_pos_neg_opt_for_paper/$(model_opt.areaname)_model_opt.jld2", model_opt)
        mean_R = mean(model_opt.MCMC_results.chain[:, 1, 1])
        R_lb = mean_R - quantile(model_opt.MCMC_results.chain[:, 1, 1], 0.025)
        R_ub = quantile(model_opt.MCMC_results.chain[:, 1, 1], 0.975) - mean_R

        return (mean_Rt = mean_R .* ct, Rt_lpred = R_lb .* ct, Rt_upred = R_ub .* ct, name = model_opt.areaname)
end


gathered_Rt_data = [add_contact_rate_and_calculate_Rt(path, gathered_smoothed_cts) for path in model_paths]
@save("analysis/gathered_Rt_data.jld2", gathered_Rt_data)
