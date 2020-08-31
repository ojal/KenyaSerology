## Load PATH and relevant packages (most of these aren't required, but I run everything in Main module
# then put into KenyaSerology module)

push!(LOAD_PATH, joinpath(homedir(),"GitHub/Kenya-Serology/src"))

using Distributions,Plots,Dates,JLD2,TransformVariables,Optim,FileIO,CSV,DataFrames
using Revise
import KenyaSerology
## Choose plotting backend. I currently like pyplot (which uses the matplotlib python backend) because
#it renders Greek symbols
pyplot()
gr()
plotlyjs()

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
## initial conditions
init_inf_nai = [nairobi_model_Peff_shortsi_var_testing.MCMC_results.chain[k,2,1] + nairobi_model_Peff_shortsi_var_testing.MCMC_results.chain[k,3,1] for k = 1:10000 ]
init_inf_mom = [mombasa_model_Peff_shortsi.MCMC_results.chain[k,2,1] + nairobi_model_Peff_shortsi_var_testing.MCMC_results.chain[k,3,1] for k = 1:10000 ]



mean(init_inf_nai), (quantile(init_inf_nai,0.025),quantile(init_inf_nai,0.975))
mean(init_inf_mom), (quantile(init_inf_mom,0.025),quantile(init_inf_mom,0.975))

## True peaks
inc_nai = KenyaSerology.incidence_across_samples(nairobi_model_Peff_shortsi_var_testing,325)
peak_nai = [argmax(inc_nai.true_incidence[:,k]) for k = 1:10000]
inc_mom = KenyaSerology.incidence_across_samples(mombasa_model_Peff_shortsi,325)
peak_mom = [argmax(inc_mom.true_incidence[:,k]) for k = 1:10000]

mean(peak_nai), (quantile(peak_nai,0.025),quantile(peak_nai,0.975))
mean(peak_mom), (quantile(peak_mom,0.025),quantile(peak_mom,0.975))



## True infections

inc_nai = KenyaSerology.incidence_across_samples(nairobi_model_Peff_shortsi_var_testing,325)
eff_trans_nai = KenyaSerology.create_credible_intervals(inc_nai.eff_transmission)

true_inf_nai = KenyaSerology.create_credible_intervals(inc_nai.true_infecteds)
sero_nai = KenyaSerology.create_credible_intervals(inc_nai.sero_converted_samples)
plot(sero_nai.mean_pred)

dayaug1 = Date(2020,8,1) - Date(2020,2,20)
100*true_inf_nai.mean_pred[163]./nairobi_model_Peff_shortsi_var_testing.N
((true_inf_nai.mean_pred[163] - true_inf_nai.lb_pred[163])*100/nairobi_model_Peff_shortsi_var_testing.N,
        (true_inf_nai.mean_pred[163] + true_inf_nai.ub_pred[163])*100/nairobi_model_Peff_shortsi_var_testing.N)

100*sero_nai.mean_pred[163]
((sero_nai.mean_pred[163] - sero_nai.lb_pred[163])*100,
        (sero_nai.mean_pred[163] + sero_nai.ub_pred[163])*100)



inc_mom = KenyaSerology.incidence_across_samples(mombasa_model_Peff_shortsi,325);
eff_trans_mom = KenyaSerology.create_credible_intervals(inc_mom.eff_transmission)

true_inf_mom = KenyaSerology.create_credible_intervals(inc_mom.true_infecteds)

eff_sus_mom = KenyaSerology.create_credible_intervals(inc_mom.eff_susceptible)
sero_mom = KenyaSerology.create_credible_intervals(inc_mom.sero_converted_samples)


dayaug1 = Date(2020,8,1) - Date(2020,2,20)
100*true_inf_mom.mean_pred[163]./mombasa_model_Peff_shortsi.N
((true_inf_mom.mean_pred[163] - true_inf_mom.lb_pred[163])*100/mombasa_model_Peff_shortsi.N,
        (true_inf_mom.mean_pred[163] + true_inf_mom.ub_pred[163])*100/mombasa_model_Peff_shortsi.N)

100*sero_mom.mean_pred[163]
((sero_mom.mean_pred[163] - sero_mom.lb_pred[163])*100,
        (sero_mom.mean_pred[163] + sero_mom.ub_pred[163])*100)

## R value plot
dayaug21 = Date(2020,8,21) - Date(2020,2,20)

@load("data/projected_contact_data_25082020.jld2")

xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:10 ]
xticklabs = [(monthname(k)[1:3]) for k = 3:12]

plt_Rt = plot(2.01*projected_contactrate_nairobi.contactrate,
        yticks = 0.5:0.25:3.5,
        ylims = (0.5,3.5),
        ls = :dash,
        xticks = (xticktimes,xticklabs),
        lab = "",
        xlims = (1,305),
        color = :red,
        guidefont = 20,
        tickfont = 14,
        titlefont = 16,
        legendfont = 10,
        size = (700,500),
        dpi = 250,
        title = "Inferred Kenyan reproductive ratios",
        ylabel = "R(t)")
plot!(plt_Rt,2.23*projected_contactrate_mombasa.contactrate,
        lab = "",
        ls = :dash,
        color = :green)

mean_R_restofkenya = 1.95
plot!(plt_Rt,mean_R_restofkenya*projected_contactrate_kenya.contactrate,
        lab = "",
        ls = :dash,
        color = :blue)
plot!(plt_Rt, 2.01*projected_contactrate_nairobi.contactrate[1:183],
        lw = 2,
        lab = "Nairobi: basic R(t)",
        color = :red)
plot!(plt_Rt, projected_contactrate_nairobi.contactrate.*eff_trans_nai.mean_pred[1:length(projected_contactrate_nairobi.contactrate)],
        lw = 1,
        ls = :dot,
        lab = "Nairobi: eff. R(t)",
        color = :red)
plot!(plt_Rt, 2.23*projected_contactrate_mombasa.contactrate[1:183],
        lw = 2,
        lab = "Mombasa: basic R(t)",
        color = :green)
plot!(plt_Rt,
        projected_contactrate_mombasa.contactrate.*eff_trans_mom.mean_pred[1:length(projected_contactrate_mombasa.contactrate)],
        lw = 1,
        ls = :dot,
        lab = "Mombasa: eff. R(t)",
        color = :green)
plot!(plt_Rt, 1.95*projected_contactrate_kenya.contactrate[1:183],
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

savefig(plt_Rt,"plotsforpaper/Rt_plot.pdf")

## Effective transmission rate

plt_effRt = plot(projected_contactrate_nairobi.contactrate.*eff_trans_nai.mean_pred[1:length(projected_contactrate_nairobi.contactrate)])

eff_trans_nai.mean_pred
