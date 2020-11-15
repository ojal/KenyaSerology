## Load PATH and relevant packages (most of these aren't required, but I run everything in Main module
# then put into KenyaSerology module)

push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaSerologyPrivate/src"))

using Distributions,Plots,Dates,JLD2,TransformVariables,Optim,FileIO,CSV,DataFrames
using Revise
import KenyaSerology
## Choose plotting backend. I currently like pyplot (which uses the matplotlib python backend) because
#it renders Greek symbols
pyplot()
gr()
plotlyjs()

## Revised submission data up to 30th September

#With R(t) fitting
@load("data/death_data_by_area_21feb_to_17oct.jld2")
@load("data/p_ID.jld2")
nai_deaths = vec(death_data.deaths[death_data.dates .<= Date(2020,10,1),death_data.areas .== "NAIROBI"])
mom_deaths = vec(death_data.deaths[death_data.dates .<= Date(2020,10,1),death_data.areas .== "MOMBASA"])
totaldeaths = sum(death_data.deaths[death_data.dates .< Date(2020,10,1),:])

@load("data/case_data_with_pos_neg_21feb_to_30sept.jld2")
sum(case_data_with_pos_neg.cases[:,case_data_with_pos_neg.areas .== "NAIROBI",1])
sum(case_data_with_pos_neg.cases[:,case_data_with_pos_neg.areas .== "MOMBASA",1])

sum(case_data_with_pos_neg.cases[:,:,:])

@load("modelfits_pos_neg/Nairobi_model.jld2")
nairobi_model = deepcopy(model)
@load("modelfits_pos_neg_opt/Nairobi_model_opt.jld2")
nairobi_model_opt = deepcopy(model_opt)
# @load("modelfits_pos_neg_opt_vs2/Nairobi_model_opt.jld2")
# nairobi_model_opt2 = deepcopy(model_opt)





@load("modelfits_pos_neg/Mombasa_model.jld2")
mombasa_model = deepcopy(model)

@load("modelfits_pos_neg_opt/Mombasa_model_opt.jld2")
mombasa_model_opt = deepcopy(model_opt)
# @load("modelfits_pos_neg_opt_vs2/Mombasa_model_opt.jld2")
# mombasa_model_opt2 = deepcopy(model_opt)

oct_first = (Date(2020,10,1) - Date(2020,2,20)).value



plt_PCR_nai = KenyaSerology.weeklyplotfittoPCRdata(nairobi_model)
# plot!(plt_PCR_mom_opt,xlims = (0,10))

plt_PCR_nai_opt = KenyaSerology.weeklyplotfittoPCRdata(nairobi_model_opt)

plt_PCR_nai_opt2 = KenyaSerology.weeklyplotfittoPCRdata(nairobi_model_opt2)

plt_PCR_mom = KenyaSerology.weeklyplotfittoPCRdata(mombasa_model)

plt_PCR_mom_opt = KenyaSerology.weeklyplotfittoPCRdata(mombasa_model_opt)
plt_PCR_mom_opt2 = KenyaSerology.weeklyplotfittoPCRdata(mombasa_model_opt2)

plot!(plt_PCR_mom,ylims = (-10,450))
plot!(plt_PCR_mom_opt,ylims = (-10,650))

savefig(plt_PCR_nai,"plotsforpaper/revisedweekly_PCR_fit_nairobi_google.pdf")
savefig(plt_PCR_nai_opt,"plotsforpaper/weekly_PCR_fit_nairobi_opt.pdf")
savefig(plt_PCR_mom,"plotsforpaper/revisedweekly_PCR_fit_mombasa_google.pdf")
savefig(plt_PCR_mom_opt,"plotsforpaper/weekly_PCR_fit_mombasa_opt.pdf")

## population plots


plt_pop_nai = KenyaSerology.population_plot(nairobi_model,Val(:monthlyserology))
plot!(plt_pop_nai,xlims = (0.,oct_first),
        ylims = (-5,120),
        guidefont = 20,
        tickfont = 14,
        legendfont = 12,
        titlefont = 18,
        yticks = [0,20,40,60,80,100])
plt_pop_nai_opt = KenyaSerology.population_plot(nairobi_model_opt,Val(:monthlyserology))
plot!(plt_pop_nai_opt,xlims = (0.,oct_first),
        ylims = (-5,120),
        guidefont = 20,
        tickfont = 14,
        legendfont = 12,
        titlefont = 18,
        yticks = [0,20,40,60,80,100])

plt_pop_mom = KenyaSerology.population_plot(mombasa_model,Val(:monthlyserology))
plot!(plt_pop_mom,xlims = (0.,oct_first),
        ylims = (-5,120),
        guidefont = 20,
        tickfont = 14,
        legendfont = 12,
        titlefont = 18,
        yticks = [0,20,40,60,80,100])
plt_pop_mom_opt = KenyaSerology.population_plot(mombasa_model_opt,Val(:monthlyserology))
plot!(plt_pop_mom_opt,xlims = (0.,oct_first),
        ylims = (-5,120),
        guidefont = 20,
        tickfont = 14,
        legendfont = 12,
        titlefont = 18,
        yticks = [0,20,40,60,80,100])
savefig(plt_pop_nai,"plotsforpaper/revised_sero_pop_nai_google.pdf")
savefig(plt_pop_nai_opt,"plotsforpaper/revised_sero_pop_nai_opt.pdf")
savefig(plt_pop_mom,"plotsforpaper/revised_sero_pop_mom_google.pdf")
savefig(plt_pop_mom_opt,"plotsforpaper/revised_sero_pop_mom_opt.pdf")

## Incidence
#
# plt_PCR = weeklyplotfittoPCRdata(nairobi_model)
# plt_PCR_nai_opt = weeklyplotfittoPCRdata(nairobi_model_opt)
# plt_PCR_mom = weeklyplotfittoPCRdata(mombasa_model)
# plt_PCR_mom_opt = weeklyplotfittoPCRdata(mombasa_model_opt)

plt_inc_mom_opt = KenyaSerology.plot_incidenceonly(mombasa_model_opt)
plot!(plt_inc_mom_opt,xlims = (0.,oct_first),
        yscale = :linear,
        guidefont = 20,
        tickfont = 14,
        legendfont = 12,
        titlefont = 18,
        title = "Mombasa: Predicted infection incidence")

savefig(plt_PCR,"weekly_PCR_fit_nairobi.pdf")
savefig(plt_PCR_mom,"weekly_PCR_fit_mombasa.pdf")
savefig(plt_PCR_nai_opt,"weekly_PCR_fit_nairobi_opt.pdf")
savefig(plt_PCR_mom_opt,"weekly_PCR_fit_mombasa_opt.pdf")


## Deaths

plt_deaths_nai = KenyaSerology.plot_deaths(nairobi_model,nai_deaths,p_ID)

plot!(plt_deaths_nai,title = "Nairobi: Observed and predicted deaths")
savefig("plotsforpaper/nairobi_deaths_vs2.pdf")

plt_deaths_mom = KenyaSerology.plot_deaths(mombasa_model_opt2,mom_deaths,p_ID)
plot!(plt_deaths_mom,title = "Mombasa: Observed and predicted deaths")
savefig("plotsforpaper/mombasa_deaths_vs2.pdf")

## Original pre-print data up to 6th August
@load("modelfits/nairobi_model_Peff_shortsi_var_testing.jld2")
@load("modelfits/mombasa_model_Peff_shortsi.jld2")

plt_case_nai = KenyaSerology.plotfittoPCRdata(nairobi_model_Peff_shortsi_var_testing)
plot!(plt_case_nai,ylims = (0,1400),
        guidefont = 20,
        tickfont = 14,
        legendfont = 12,
        titlefont = 19,
        title = "Nairobi: PCR positive test prediction and data")
savefig(plt_case_nai,"plotsforpaper/nairobi_cases_vs2.pdf")
plt_case_mom = KenyaSerology.plotfittoPCRdata(mombasa_model_Peff_shortsi)
plot!(plt_case_mom,ylims = (0,125),
        guidefont = 20,
        tickfont = 14,
        legendfont = 12,
        titlefont = 18,
        title = "Mombasa: PCR positive test prediction and data")
savefig(plt_case_mom,"plotsforpaper/mombasa_cases_vs2.pdf")

plt_pop_nai = KenyaSerology.population_plot(nairobi_model_Peff_shortsi_var_testing,Val(:monthlyserology))
plot!(plt_pop_nai,ylims = (-5,120),
        guidefont = 20,
        tickfont = 14,
        legendfont = 12,
        titlefont = 18,
        yticks = [0,20,40,60,80,100])
savefig(plt_pop_nai,"plotsforpaper/nairobi_sero_vs2.pdf")

plt_pop_mom = KenyaSerology.population_plot(mombasa_model_Peff_shortsi,Val(:monthlyserology))
plot!(plt_pop_mom,ylims = (-5,120),
        guidefont = 20,
        tickfont = 14,
        legendfont = 12,
        titlefont = 18,
        yticks = [0,20,40,60,80,100])
savefig(plt_pop_mom,"plotsforpaper/mombasa_sero_vs2.pdf")

plt_inc_nai = KenyaSerology.plot_incidenceonly(nairobi_model_Peff_shortsi_var_testing)
plot!(plt_inc_nai,
        guidefont = 20,
        tickfont = 14,
        legendfont = 12,
        titlefont = 18,
        title = "Nairobi: Predicted infection incidence")
savefig("plotsforpaper/nairobi_inc_only_vs2.pdf")

plt_inc_mom = KenyaSerology.plot_incidenceonly(mombasa_model_Peff_shortsi)
plot!(plt_inc_mom,
        guidefont = 20,
        tickfont = 14,
        legendfont = 12,
        titlefont = 18,
        title = "Mombasa: Predicted infection incidence")
savefig("plotsforpaper/mombasa_inc_only_vs2.pdf")

plt_deaths_nai = KenyaSerology.plot_deaths(nairobi_model_Peff_shortsi_var_testing,nai_death,p_ID)
plot!(plt_deaths_nai,title = "Nairobi: Observed and predicted deaths")
savefig("plotsforpaper/nairobi_deaths_vs2.pdf")

plt_deaths_mom = KenyaSerology.plot_deaths(mombasa_model_Peff_shortsi,mom_death,p_ID)
plot!(plt_deaths_mom,title = "Mombasa: Observed and predicted deaths")
savefig("plotsforpaper/mombasa_deaths_vs2.pdf")

## Daily deaths
plt_nai_deaths = KenyaSerology.plot_deaths_weekly(nairobi_model,nai_deaths,p_ID,55)
savefig(plt_nai_deaths,"plotsforpaper/revisednairobi_deaths_only_google.pdf")
plt_mom_deaths = KenyaSerology.plot_deaths_weekly(mombasa_model,mom_deaths,p_ID,20)
savefig(plt_mom_deaths,"plotsforpaper/revisedmombasa_deaths_only_google.pdf")


plt_nai_deaths_opt = KenyaSerology.plot_deaths_weekly(nairobi_model_opt,nai_deaths,p_ID)
savefig("plotsforpaper/nairobi_deaths_only_vs2.pdf")
plt_mom_deaths_opt = KenyaSerology.plot_deaths_weekly(mombasa_model_opt,mom_deaths,p_ID)
savefig("plotsforpaper/mombasa_deaths_only_vs2.pdf")

## Daily incidence
plt_nai_inc = KenyaSerology.plot_incidenceonly(nairobi_model)
savefig(plt_nai_inc,"plotsforpaper/revisednairobi_inc_only_google.pdf")
plt_mom_inc = KenyaSerology.plot_incidenceonly(mombasa_model)
savefig(plt_mom_inc,"plotsforpaper/revisedmombasa_inc_only_google.pdf")

plt_nai_inc_opt = KenyaSerology.plot_incidenceonly(nairobi_model_opt)
savefig("plotsforpaper/nairobi_inc_only_vs2.pdf")
plt_mom_inc_opt = KenyaSerology.plot_incidenceonly(mombasa_model_opt)
savefig("plotsforpaper/mombasa_inc_only_vs2.pdf")



## initial conditions
init_inf_nai = [nairobi_model_Peff_shortsi_var_testing.MCMC_results.chain[k,2,1] + nairobi_model_Peff_shortsi_var_testing.MCMC_results.chain[k,3,1] for k = 1:10000 ]
init_inf_mom = [mombasa_model_Peff_shortsi.MCMC_results.chain[k,2,1] + nairobi_model_Peff_shortsi_var_testing.MCMC_results.chain[k,3,1] for k = 1:10000 ]



mean(init_inf_nai), (quantile(init_inf_nai,0.025),quantile(init_inf_nai,0.975))
mean(init_inf_mom), (quantile(init_inf_mom,0.025),quantile(init_inf_mom,0.975))

## True peaks
inc_nai = KenyaSerology.incidence_across_samples(nairobi_model_opt,325)
peak_nai = [argmax(inc_nai.true_incidence[:,k]) for k = 1:10000]


inc_mom = KenyaSerology.incidence_across_samples(mombasa_model_opt,325)
peak_mom = [argmax(inc_mom.true_incidence[:,k]) for k = 1:10000]

mean(peak_nai), (quantile(peak_nai,0.025),quantile(peak_nai,0.975))
mean(peak_mom), (quantile(peak_mom,0.025),quantile(peak_mom,0.975))
mombasa_model_opt.MCMC_results.chain

display(nairobi_model_opt.MCMC_results.chain)
mean(nairobi_model_opt.MCMC_results.chain[:,end,1])
quantile(nairobi_model_opt.MCMC_results.chain[:,end,1],0.025)
quantile(nairobi_model_opt.MCMC_results.chain[:,end,1],0.975)
mean(mombasa_model_opt.MCMC_results.chain[:,end,1])
quantile(mombasa_model_opt.MCMC_results.chain[:,end,1],0.025)
quantile(mombasa_model_opt.MCMC_results.chain[:,end,1],0.975)
display(mombasa_model_opt.MCMC_results.chain)

function doubling_time(R₀,σ,γ )
        r = 0.5*( sqrt((4*(R₀ - 1)*σ*γ) + (σ + γ)^2) - (σ+γ))
        return log(2)/r
end

doubling_time(2,1/3.1,1/2.4)


## True infections

inc_nai = KenyaSerology.incidence_across_samples(nairobi_model_opt,325)
eff_trans_nai = KenyaSerology.create_credible_intervals(inc_nai.eff_transmission)
eff_sus_nai = KenyaSerology.create_credible_intervals(inc_nai.eff_susceptible)

true_inf_nai = KenyaSerology.create_credible_intervals(inc_nai.true_infecteds)
sero_nai = KenyaSerology.create_credible_intervals(inc_nai.sero_converted)
plot(sero_nai.mean_pred)


c_t_nai_array = inc_nai.true_incidence./(nairobi_model_opt.γ.*inc_nai.eff_transmission.*inc_nai.true_infectious)
mean_ct_nai = mean(c_t_nai_array,dims = 2)
mean_ct_nai[1] = mean_ct_nai[2]

@load("data/projected_contact_data_09102020.jld2")
plot(mean_ct_nai)
plot!(projected_contactrate_nairobi.contactrate)

dayaug1 = Date(2020,8,1) - Date(2020,2,20)
sept_30 = Date(2020,9,30) - Date(2020,2,20)

100*true_inf_nai.mean_pred[223]./nairobi_model.N
((true_inf_nai.mean_pred[223] - true_inf_nai.lb_pred[223])*100/nairobi_model.N,
        (true_inf_nai.mean_pred[223] + true_inf_nai.ub_pred[223])*100/nairobi_model.N)

100*sero_nai.mean_pred[223]
((sero_nai.mean_pred[223] - sero_nai.lb_pred[223])*100,
        (sero_nai.mean_pred[223] + sero_nai.ub_pred[223])*100)



inc_mom = KenyaSerology.incidence_across_samples(mombasa_model_opt,325);
eff_trans_mom = KenyaSerology.create_credible_intervals(inc_mom.eff_transmission)

true_inf_mom = KenyaSerology.create_credible_intervals(inc_mom.true_infecteds)

eff_sus_mom = KenyaSerology.create_credible_intervals(inc_mom.eff_susceptible)

c_t_mom_array = inc_mom.true_incidence./(nairobi_model_opt.γ.*inc_mom.eff_transmission.*inc_mom.true_infectious)
mean_ct_mom = mean(c_t_mom_array,dims = 2)
mean_ct_mom[1] = mean_ct_mom[2]

plot(mean_ct_mom)

sero_mom = KenyaSerology.create_credible_intervals(inc_mom.sero_converted)


dayaug1 = Date(2020,8,1) - Date(2020,2,20)
100*true_inf_mom.mean_pred[223]./mombasa_model_opt.N
((true_inf_mom.mean_pred[223] - true_inf_mom.lb_pred[223])*100/mombasa_model_opt.N,
        (true_inf_mom.mean_pred[223] + true_inf_mom.ub_pred[223])*100/mombasa_model_opt.N)

100*sero_mom.mean_pred[223]
((sero_mom.mean_pred[223] - sero_mom.lb_pred[223])*100,
        (sero_mom.mean_pred[223] + sero_mom.ub_pred[223])*100)

## R value plot
dayaug21 = Date(2020,8,21) - Date(2020,2,20)

@load("data/projected_contact_data_09102020.jld2")

xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:8 ]
xticklabs = [(monthname(k)[1:3]) for k = 3:10]

mean_R_nai = mean(nairobi_model_opt.MCMC_results.chain[:,1,1])
lb_R_nai = mean_R_nai - quantile(nairobi_model_opt.MCMC_results.chain[:,1,1],0.025)
ub_R_nai = quantile(nairobi_model_opt.MCMC_results.chain[:,1,1],0.975) - mean_R_nai

mean_R_mom = mean(mombasa_model_opt.MCMC_results.chain[:,1,1])


plt_Rt = plot(mean_R_nai*projected_contactrate_nairobi.contactrate[1:224],
        yticks = 0.5:0.25:3.5,
        ylims = (0.5,2.5),
        # ribbon = (lb_R_nai*projected_contactrate_nairobi.contactrate[1:224],ub_R_nai*projected_contactrate_nairobi.contactrate[1:224]),
        xticks = (xticktimes,xticklabs),
        lab = "Google estimate: Nairobi",
        xlims = (1,224),
        color = :red,
        ls = :dash,
        guidefont = 20,
        tickfont = 14,
        titlefont = 16,
        legendfont = 10,
        size = (700,500),
        dpi = 250,
        title = "Inferred Kenyan reproductive ratios",
        ylabel = "R(t)")
plot!(plt_Rt,mean_R_mom*projected_contactrate_mombasa.contactrate[1:224],
        lab = "Google estimate: Mombasa",
        ls = :dash,
        color = :green)
plot!(plt_Rt,mean_R_nai*mean_ct_nai[1:224],
        lab = "fitted R(t): Nairobi",lw = 2.5,
        color = :red)
plot!(plt_Rt,mean_R_mom*mean_ct_mom[1:224],
        lab = "fitted R(t): Mombasa",lw = 2.5,
        color = :green)

plt_effRt = plot(eff_trans_nai.mean_pred[1:224].*projected_contactrate_nairobi.contactrate[1:224],
        yticks = 0.5:0.25:3.5,
        ylims = (0.5,2.5),
        # ribbon = (eff_trans_nai.lb_pred[1:224].*projected_contactrate_nairobi.contactrate[1:224],eff_trans_nai.ub_pred[1:224].*projected_contactrate_nairobi.contactrate[1:224]),
        xticks = (xticktimes,xticklabs),
        lab = "Google estimate: Nairobi",
        xlims = (1,224),
        color = :red,
        ls = :dash,
        guidefont = 20,
        tickfont = 14,
        titlefont = 16,
        legendfont = 10,
        size = (700,500),
        dpi = 250,
        title = "Effective reproductive ratios",
        ylabel = "Eff. R(t)")
plot!(plt_effRt,eff_trans_mom.mean_pred[1:224].*projected_contactrate_mombasa.contactrate[1:224],
        # ribbon = (eff_trans_mom.lb_pred[1:224].*projected_contactrate_mombasa.contactrate[1:224],eff_trans_mom.ub_pred[1:224].*projected_contactrate_mombasa.contactrate[1:224]),
        lab = "Google estimate: Mombasa",
        ls = :dash,
        color = :green)
plt_effRt = plot(eff_trans_nai.mean_pred[1:224].*mean_ct_nai[1:224],
        ribbon = (eff_trans_nai.lb_pred[1:224].*mean_ct_nai[1:224],eff_trans_nai.ub_pred[1:224].*mean_ct_nai[1:224]),
        lab = "fitted R(t): Nairobi",lw = 2.5,
        yticks = 0.5:0.25:3.5,
        ylims = (0.5,2.5),
        # ribbon = (eff_trans_nai.lb_pred[1:224].*projected_contactrate_nairobi.contactrate[1:224],eff_trans_nai.ub_pred[1:224].*projected_contactrate_nairobi.contactrate[1:224]),
        xticks = (xticktimes,xticklabs),
        xlims = (1,224),
        color = :red,
        guidefont = 20,
        tickfont = 14,
        titlefont = 16,
        legendfont = 10,
        size = (700,500),
        dpi = 250,
        title = "Effective reproductive ratios",
        ylabel = "Eff. R(t)")
plot!(plt_effRt,eff_trans_mom.mean_pred[1:224].*mean_ct_mom[1:224],
        ribbon = (eff_trans_mom.lb_pred[1:224].*mean_ct_mom[1:224],eff_trans_mom.ub_pred[1:224].*mean_ct_mom[1:224]),
        lab = "fitted R(t): Mombasa",lw = 2.5,
        color = :green)

mean_R_restofkenya = 1.95
# plot!(plt_Rt,mean_R_restofkenya*projected_contactrate_kenya.contactrate,
#         lab = "",
#         ls = :dash,
#         color = :blue)
plot!(plt_effRt, 2.01*projected_contactrate_nairobi.contactrate[1:183],
        lw = 2,
        lab = "Nairobi: basic R(t)",
        color = :red)
plot!(plt_effRt, projected_contactrate_nairobi.contactrate.*eff_trans_nai.mean_pred[1:length(projected_contactrate_nairobi.contactrate)],
        lw = 1,
        ls = :dot,
        lab = "Nairobi: eff. R(t)",
        color = :red)
plot!(plt_effRt, 2.23*projected_contactrate_mombasa.contactrate[1:183],
        lw = 2,
        lab = "Mombasa: basic R(t)",
        color = :green)
plot!(plt_effRt,
        projected_contactrate_mombasa.contactrate.*eff_trans_mom.mean_pred[1:length(projected_contactrate_mombasa.contactrate)],
        lw = 1,
        ls = :dot,
        lab = "Mombasa: eff. R(t)",
        color = :green)
plot!(plt_effRt, 1.95*projected_contactrate_kenya.contactrate[1:183],
        lw = 2,
        lab = "Rest of Kenya: mean basic R(t)",
        color = :blue)

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

plot!(plt_effRt,[firstcase,firstcase],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "Restrictions")
plot!(plt_effRt,[suspendpublicgathering,suspendpublicgathering],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")
plot!(plt_effRt,[schoolsclosed,schoolsclosed],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")
plot!(plt_effRt,[retroactivequarantine,retroactivequarantine],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")
plot!(plt_effRt,[inboundflightsuspended,inboundflightsuspended],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")
plot!(plt_effRt,[nationalcurfew,nationalcurfew],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")
plot!(plt_effRt,[regionallockdowns,regionallockdowns],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")
plot!(plt_effRt,[localnaiandmomlockdowns,localnaiandmomlockdowns],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")
plot!(plt_effRt,[tazandsomborderclosed,tazandsomborderclosed],[0.5,3.5],
        ls = :dash,lw=0.5,color = :black,lab = "")

ann1 = (19, 0.76, Plots.text("1", :left))
ann2 = (39, 0.76, Plots.text("2", :left))
ann3 = (79, 0.76, Plots.text("3", :left))

plot!(plt_effRt,annotations = [ann1,ann2,ann3])

savefig(plt_effRt,"plotsforpaper/revisedRt_plot.pdf")

## Effective transmission rate

plt_effRt = plot(projected_contactrate_nairobi.contactrate.*eff_trans_nai.mean_pred[1:length(projected_contactrate_nairobi.contactrate)])

eff_trans_nai.mean_pred
