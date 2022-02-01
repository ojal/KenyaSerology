function groupserobyweek(areamodel::KenyaSerology.CoVAreaModel)
    testingmonths = week.(areamodel.dates[1:size(areamodel.sero_cases,1)])
    unique_weeks = unique(testingweeks)
    posbyweek = [sum(areamodel.sero_cases[testingweeks.==w,1]) for w in unique_weeks]
    numbyweek = [sum(areamodel.sero_cases[testingweeks.==w,2]) for w in unique_weeks] .+ posbyweek
    midpointsofweek = [((Date(2019,12,31) + Week(w) - Day(3)) - Date(2020,2,21)).value for w in unique_weeks]
    upred = zeros(length(posbyweek))
    lpred = similar(upred)
    for i = 1:length(upred)
        if numbyweek[i] > 0
            d = Beta(posbyweek[i] +0.5,numbyweek[i] - posbyweek[i] + 0.5) #Jeffrey's interval
            upred[i] = invlogcdf(d,log(0.995))- ( posbyweek[i]/numbyweek[i] )
            lpred[i] = ( posbyweek[i]/numbyweek[i] ) - invlogcdf(d,log(0.005))
        else
            upred[i] = 0.
            lpred[i] = 0.
        end
    end

    return (posbyweek=posbyweek,
            numbyweek=numbyweek,
            midpointsofweek=midpointsofweek,
            upred = upred,
            lpred=lpred)
end

function groupserobymonth(areamodel::KenyaSerology.CoVAreaModel)
    testingmonths = month.(areamodel.dates[1:size(areamodel.sero_cases,1)])
    unique_months = unique(testingmonths)
    posbymonth = [sum(areamodel.sero_cases[testingmonths.==m,1]) for m in unique_months]
    numbymonth = [sum(areamodel.sero_cases[testingmonths.==m,2]) for m in unique_months] .+ posbymonth
    midpointsofmonth= [((Date(2019,12,31) + Month(m) - Day(15)) - Date(2020,2,21)).value for m in unique_months]
    upred = zeros(length(posbymonth))
    lpred = similar(upred)
    for i = 1:length(upred)
        if numbymonth[i] > 0
            d = Beta(posbymonth[i] +0.5,numbymonth[i] - posbymonth[i] + 0.5) #Jeffrey's interval
            upred[i] = invlogcdf(d,log(0.995))- ( posbymonth[i]/numbymonth[i] )
            lpred[i] = ( posbymonth[i]/numbymonth[i] ) - invlogcdf(d,log(0.005))
        else
            upred[i] = 0.
            lpred[i] = 0.
        end
    end

    return (posbymonth=posbymonth,
            numbymonth=numbymonth,
            midpointsofmonth=midpointsofmonth,
            upred = upred,
            lpred=lpred)
end

function seroweekaverage(areamodel::KenyaSerology.CoVAreaModel)
        dayswithtests = vec(sum(areamodel.sero_cases,dims = 2)) .> 0
        f₁ = findfirst(dayswithtests)
        f₂ = findlast(dayswithtests)
        days = (f₁+3):(f₂-3)
        posbyweek = [sum(areamodel.sero_cases[(d-3):(d+3),1]) for d in days]
        numbyweek = [sum(areamodel.sero_cases[(d-3):(d+3),2]) for d in days] .+ posbyweek
        midpointsofweek = [ d-1 for d in days]
        upred = zeros(length(posbyweek))
        lpred = similar(upred)
    for i = 1:length(upred)
        if numbyweek[i] > 0
            d = Beta(posbyweek[i] +0.5,numbyweek[i] - posbyweek[i] + 0.5) #Jeffrey's interval
            upred[i] = invlogcdf(d,log(0.975))- ( posbyweek[i]/numbyweek[i] )
            lpred[i] = ( posbyweek[i]/numbyweek[i] ) - invlogcdf(d,log(0.025))
        else
            upred[i] = 0.
            lpred[i] = 0.
        end
    end

    return (posbyweek=posbyweek,
            numbyweek=numbyweek,
            midpointsofweek=midpointsofweek,
            upred = upred,
            lpred=lpred)
end

function death_pred(incidence,IFR,p_ID)
    deaths = IFR*simple_conv(incidence,p_ID)
end

function death_pred_err(obs_deaths,incidence,IFR,p_ID)
    err = 0.
    deaths = death_pred(incidence,IFR,p_ID)
    for (i,true_deaths) in enumerate(obs_deaths)
        err -= logpdf(Poisson(deaths[i]),true_deaths)
    end
    return err
end




function create_death_prediction_intervals(death_incidence)
    time_length = size(death_incidence,1)
    mean_pred = zeros(time_length)
    lb_pred = zeros(time_length)
    ub_pred = zeros(time_length)
    for t = 1:time_length
        mean_pred[t] =death_incidence[t]
        lb_pred[t] = quantile(Poisson(death_incidence[t]),0.025)
        ub_pred[t] = quantile(Poisson(death_incidence[t]),0.975)
    end
    return (mean_pred=mean_pred,lb_pred=mean_pred.-lb_pred,ub_pred=ub_pred.-mean_pred)
end

"""
    cornerplot(areamodel::CoVAreaModel)

Output a corner plot for the area model (post-inference)
"""
function cornerplot(areamodel::CoVAreaModel)
    plt = corner(areamodel.MCMC_results.chain,size = (800,600),dpi = 200,plot_title = "MCMC")
end

"""
    incidence_across_samples(areamodel::CoVAreaModel,l)

Solve the underlying ODE model for each of the sampled realisations of the
parameters. This provides a posterior distribution for incidences (sero, PCR, etc)
"""
function incidence_across_samples(areamodel::CoVAreaModel,l)
        n = size(areamodel.MCMC_results.chain,1)
        m = size(areamodel.MCMC_results.chain,2)
        paramnames = areamodel.MCMC_results.chain.name_map.parameters
        rel_test = vcat(areamodel.relative_testing_rate,areamodel.relative_testing_rate[end]*ones(l))

        PCR_incidence_NB = zeros(Int64(l),n)
        PCR_P_BB = zeros(Int64(l),n)
        sero_converted = zeros(Int64(l),n)
        true_infecteds = zeros(Int64(l),n)
        true_incidence = zeros(Int64(l),n)
        true_susceptible = zeros(Int64(l),n)
        eff_susceptible = zeros(Int64(l),n)
        true_infectious = zeros(Int64(l),n)
        eff_transmission = zeros(Int64(l),n)


        for i = 1:n
            #assign to correct values
            #Parameters that exist across models
            R = areamodel.MCMC_results.chain[i,findfirst(paramnames .== :R),1]
            E₀ = areamodel.MCMC_results.chain[i,findfirst(paramnames .== :E₀),1]
            I₀ = areamodel.MCMC_results.chain[i,findfirst(paramnames .== :I₀),1]
            α = 0.01
            p_test = 1e-5
            χ = 1.
            P_eff = 1.
            M_PCR = 20.
            if :α ∈ paramnames
                α = areamodel.MCMC_results.chain[i,findfirst(paramnames .== :α),1]
            end
            if :P_eff ∈ paramnames
                P_eff = areamodel.MCMC_results.chain[i,findfirst(paramnames .== :P_eff),1]
            end
            if :χ ∈ paramnames
                χ = areamodel.MCMC_results.chain[i,findfirst(paramnames .== :χ),1]
            end
            if :p_test ∈ paramnames
                p_test = areamodel.MCMC_results.chain[i,findfirst(paramnames .== :p_test),1]
            end
            if :M_PCR ∈ paramnames
                M_PCR = areamodel.MCMC_results.chain[i,findfirst(paramnames .== :M_PCR),1]
            end


            p = [R,areamodel.σ,areamodel.γ,areamodel.N*P_eff]
            # x₀ = [(areamodel.N*P_eff) - E₀ - I₀,E₀,I₀,0.,0.]
            x₀ = [(areamodel.N*P_eff),E₀,I₀,0.,0.]
            sol = solve(areamodel.prob,BS3();u0 = x₀,p = p,saveat=1,tspan = (0.,l) )
            incidence = get_incidence(sol)
            PCR⁺ = simple_conv(incidence,areamodel.PCR_array)

            #Neg. Bin mean
            PCR⁺_pred_NB = PCR⁺.*p_test.*rel_test[1:length(PCR⁺)]
            #Betabinomial P
            PCR⁺_pred_BB = χ.*PCR⁺./((χ-1).*PCR⁺ .+ areamodel.N)

            p_sero = areamodel.sero_sensitivity*simple_conv(incidence,areamodel.sero_array)./areamodel.N
            PCR_incidence_NB[:,i] = PCR⁺_pred_NB
            PCR_P_BB[:,i] = PCR⁺_pred_BB
            sero_converted[:,i] = p_sero .+ (1 .- p_sero).*(1-areamodel.sero_specificity)
            true_incidence[:,i] = incidence

            for t = 1:Int64(l)
                true_infecteds[t,i] = sum(sol.u[t][2:4])
                true_susceptible[t,i] = (sol.u[t][1]) + areamodel.N*(1-P_eff)
                eff_susceptible[t,i] = sol.u[t][1]/(areamodel.N*P_eff)
                eff_transmission[t,i] = R*eff_susceptible[t,i]
                true_infectious[t,i] = sol.u[t][3]
            end
        end
        return (PCR_incidence_NB=PCR_incidence_NB,
                PCR_P_BB = PCR_P_BB,
                sero_converted=sero_converted,
                true_infecteds=true_infecteds,
                true_incidence=true_incidence,
                true_susceptible =true_susceptible,
                true_infectious=true_infectious,
                eff_susceptible=eff_susceptible,
                eff_transmission=eff_transmission)
end



"""
    function draw_neg_bin(μ,α)

Draws a negative binomial r.v. with mean μ and clustering coefficient α
"""
function draw_neg_bin(μ,α)
    if μ > 0
        σ² = μ + α*μ^2
        p_negbin = 1 - (α*μ^2/σ²)
        r_negbin = 1/α
        return rand(NegativeBinomial(r_negbin,p_negbin))
    else
        return 0
    end

end

function draw_beta_bin(P,M,n)
        return rand(BetaBinomial(n,M*P,M*(1-P)))
end


"""
    function neg_bin_sampling(mean_samples,α_samples,num_draws)

Draws num_draws NegativeBinomial r.v.s using parameters drawn uniformly from the distribution
of means and clustering coefficients.
"""
function neg_bin_sampling(mean_samples,α_samples,num_draws)
    draws = zeros(Int64,num_draws)
    n = length(mean_samples)
    for i = 1:num_draws
        k = rand(1:n)
        draws[i] = draw_neg_bin(mean_samples[k],α_samples[k])
    end
    return draws
end
"""
    draw_week_PCR_NB(pred_mean_PCR_for_week,αs)

Draw a weeks worth of PCR samples from the fit.
"""
function draw_week_PCR_NB(pred_mean_PCR_for_week,αs)
    k = rand(1:size(pred_mean_PCR_for_week,2))
    return sum([KenyaSerology.draw_neg_bin(pred_mean_PCR_for_week[t,k],αs[k]) for t in 1:size(pred_mean_PCR_for_week,1)])
end

function draw_week_PCR_BB(pred_P_PCR_for_week,Ms,ns)
    k = rand(1:size(pred_P_PCR_for_week,2))
    return sum([KenyaSerology.draw_beta_bin(pred_P_PCR_for_week[t,k].+ 0.0001,Ms[k],ns[t]) for t in 1:size(pred_P_PCR_for_week,1)])
end

function create_weekly_prediction_intervals_NB(pred_mean_PCR_for_week,αs,num_samples)
    samples = [KenyaSerology.draw_week_PCR_NB(pred_mean_PCR_for_week,αs) for i = 1:num_samples]
    return (mean = mean(samples),lb = mean(samples) - quantile(samples,0.025),ub = quantile(samples,0.975) - mean(samples))
end

function create_weekly_prediction_intervals_BB(pred_P_PCR_for_week,Ms,ns,num_samples)
    samples = [KenyaSerology.draw_week_PCR_BB(pred_P_PCR_for_week,Ms,ns) for i = 1:num_samples]
    return (mean = mean(samples),lb = mean(samples) - quantile(samples,0.025),ub = quantile(samples,0.975) - mean(samples))
end

function create_all_weekly_prediction_intervals_NB(pred_mean_PCR,αs,num_samples)
    week_on_each_day = [Week(Date(2020,2,20) + Day(k)).value for k = 1:size(pred_mean_PCR,1)]
    week_range = minimum(week_on_each_day):maximum(week_on_each_day)
    weekly_pred = zeros(length(week_range))
    weekly_lb = zeros(length(week_range))
    weekly_ub = zeros(length(week_range))
    for (i,week_num) in enumerate(week_range)
        z = @view pred_mean_PCR[week_on_each_day .== week_num,:]
        weekly_sample = KenyaSerology.create_weekly_prediction_intervals_NB(z,αs,num_samples)
        weekly_pred[i] = weekly_sample.mean
        weekly_lb[i] = weekly_sample.lb
        weekly_ub[i] = weekly_sample.ub
    end
    return (mean_pred=weekly_pred,weekly_lb=weekly_lb,weekly_ub=weekly_ub)
end

function create_all_weekly_prediction_intervals_BB(pred_P_PCR,Ms,n_alldates,num_samples)
    week_on_each_day = [Week(Date(2020,2,20) + Day(k)).value for k = 1:size(pred_P_PCR,1)]
    week_range = minimum(week_on_each_day):maximum(week_on_each_day)
    weekly_pred = zeros(length(week_range))
    weekly_lb = zeros(length(week_range))
    weekly_ub = zeros(length(week_range))
    for (i,week_num) in enumerate(week_range)
        z = @view pred_P_PCR[week_on_each_day .== week_num,:]
        ns = @view n_alldates[week_on_each_day .== week_num]
        weekly_sample = KenyaSerology.create_weekly_prediction_intervals_BB(z,Ms,ns,num_samples)
        weekly_pred[i] = weekly_sample.mean
        weekly_lb[i] = weekly_sample.lb
        weekly_ub[i] = weekly_sample.ub
    end
    return (mean_pred=weekly_pred,weekly_lb=weekly_lb,weekly_ub=weekly_ub)
end

function group_PCR_pos_by_week(PCR_pos)
    week_on_each_day = [Week(Date(2020,2,20) + Day(k)).value for k = 1:size(PCR_pos,1)]
    week_range = minimum(week_on_each_day):maximum(week_on_each_day)
    weekly_PCR_grouping = zeros(length(week_range))
    for (i,week_num) in enumerate(week_range)
        z = @view PCR_pos[week_on_each_day .== week_num]
        weekly_PCR_grouping[i] = sum(z)
    end
    return weekly_PCR_grouping
end

"""
    function create_PCR_prediction_intervals(PCR_incidence,α_samples)

Draws a prediction interval for observed PCR detections (this includes the uncertainty around sampling).
"""
function create_PCR_prediction_intervals(PCR_incidence,α_samples)
    time_length = size(PCR_incidence,1)
    mean_pred = zeros(time_length)
    lb_pred = zeros(time_length)
    ub_pred = zeros(time_length)
    for t = 1:time_length
        samples = neg_bin_sampling(PCR_incidence[t,:],α_samples,100000)
        mean_pred[t] = mean(samples)
        lb_pred[t] = quantile(samples,0.025)
        ub_pred[t] = quantile(samples,0.975)
    end
    return (mean_pred=mean_pred,lb_pred=mean_pred.-lb_pred,ub_pred=ub_pred.-mean_pred)
end

function create_PCR_credible_intervals(PCR_incidence)
    time_length = size(PCR_incidence,1)
    mean_pred = zeros(time_length)
    lb_pred = zeros(time_length)
    ub_pred = zeros(time_length)
    for t = 1:time_length
        PCRonthatday = PCR_incidence[t,:]
        mean_pred[t] = mean(PCRonthatday)
        lb_pred[t] = quantile(PCRonthatday,0.025)
        ub_pred[t] = quantile(PCRonthatday,0.975)
    end
    return (mean_pred=mean_pred,lb_pred=mean_pred.-lb_pred,ub_pred=ub_pred.-mean_pred)
end

"""
    function create_credible_intervals(inc::AbstractMatrix)

Create credible intervals for the `AbstractMatrix` `inc`, where dim 1 is a predicted value on each day and dim 2 is over
    posterior parameter draws.
"""
function create_credible_intervals(inc::AbstractMatrix)
    time_length = size(inc,1)
    mean_pred = zeros(time_length)
    lb_pred = zeros(time_length)
    ub_pred = zeros(time_length)
    for t = 1:time_length
        inconthatday = inc[t,:]
        mean_pred[t] = mean(inconthatday)
        lb_pred[t] = quantile(inconthatday,0.025)
        ub_pred[t] = quantile(inconthatday,0.975)
    end
    return (mean_pred=mean_pred,lb_pred=mean_pred.-lb_pred,ub_pred=ub_pred.-mean_pred)
end

"""
function R_t(areamodel::CoVAreaModel)

This method calculates median R_t and 95% PI width around the median.
"""
function R_t(areamodel::CoVAreaModel)
        R_t_array = zeros(length(areamodel.contactrate_data.contactrate),size(areamodel.MCMC_results.chain,1))
        for t = 1:size(R_t_array,1), n = 1:size(R_t_array,2)
                R_t_array[t,n] = areamodel.contactrate_data.contactrate[t]*areamodel.MCMC_results.chain[n,1,1]
        end
        R_t_med = median(R_t_array,dims = 2)
        R_t_ub = [quantile(R_t_array[t,:],0.975) for t = 1:size(R_t_array,1)] .- R_t_med
        R_t_lb = R_t_med .-[quantile(R_t_array[t,:],0.025) for t = 1:size(R_t_array,1)]

        return (R_t_med = R_t_med,R_t_ub=R_t_ub,R_t_lb=R_t_lb)
end

function weeklyplotfittoPCRdata(model::KenyaSerology.CoVAreaModel)
    xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value  for k = 1:8 ]
    xticklabs = [(monthname(k)[1:3]) for k = 3:10]

    sept_30 = (Date(2020,9,30) - Date(2020,2,20)).value
    inc = KenyaSerology.incidence_across_samples(model,sept_30)
    pred_mean_PCR = inc.PCR_incidence_NB
    pred_P_PCR = inc.PCR_P_BB

    week_on_each_day = [Week(Date(2020,2,20) + Day(k)).value for k = 1:size(pred_mean_PCR,1)]
    firstweek_with_neg_tests = week_on_each_day[findfirst(model.PCR_cases[:,2] .> 0)] - minimum(week_on_each_day)


    Ms = model.MCMC_results.chain[:,findfirst(keys(model.MCMC_results.chain).== :M_PCR),1]
    αs = model.MCMC_results.chain[:,findfirst(keys(model.MCMC_results.chain).== :α),1]
    total_ns = vec(sum(model.PCR_cases,dims=2))
    first_7days_with_neg_tests = collect(1:7)
    if !isempty(findall(model.PCR_cases[:,2] .< 0))
        first_7days_with_neg_tests = findall(model.PCR_cases[:,2] .< 0)[1:7]
    end
    guess_for_total_tests_before_reporting = max(round(Int64,sum(model.PCR_cases[first_7days_with_neg_tests,:])/7),0)
    total_ns[model.PCR_cases[:,2] .< 0] .= guess_for_total_tests_before_reporting# round(Int64,mean(total_ns[model.PCR_cases[:,2] .> 0]))

    weekly_preds_NB = KenyaSerology.create_all_weekly_prediction_intervals_NB(pred_mean_PCR,αs,10000)
    weekly_preds_BB = KenyaSerology.create_all_weekly_prediction_intervals_BB(pred_P_PCR,Ms,total_ns,10000)
    weekly_grouped_pos = KenyaSerology.group_PCR_pos_by_week(model.PCR_cases[:,1])

    #First week is week 8 of year which starts 17th Feb i.e. -4 on the tick scale
    weekdates_on_scale = collect(0:(length(weekly_preds_NB.mean_pred))-1).*7 .- 4

    plt = plot(weekdates_on_scale[1:(firstweek_with_neg_tests+1)] ,weekly_preds_NB.mean_pred[1:(firstweek_with_neg_tests+1)],
                ribbon = (weekly_preds_NB.weekly_lb[1:(firstweek_with_neg_tests+1)],weekly_preds_NB.weekly_ub[1:(firstweek_with_neg_tests+1)]),
                xticks = (xticktimes,xticklabs),
                lab = "Only +ve swab tests available",
                legend = :topleft,
                size = (700,500),dpi = 250,
                guidefontsize = 20,
                tickfontsize = 14,
                legendfontsize = 12,
                titlefontsize = 19,
                left_margin = 5mm,right_margin = 10mm,
                title = "$(model.areaname): PCR positive test prediction and data",
                ylabel = "Weekly positive swab tests ",
                lw = 2.5)
    plot!(plt,weekdates_on_scale[(firstweek_with_neg_tests+1):(end-1)] ,weekly_preds_BB.mean_pred[(firstweek_with_neg_tests+1):(end-1)],
            ribbon = (weekly_preds_BB.weekly_lb[(firstweek_with_neg_tests+1):end],weekly_preds_BB.weekly_ub[(firstweek_with_neg_tests+1):end]),
            lab = "All tests available",
            lw = 2.5)

    scatter!(plt,weekdates_on_scale[1:(end-1)],weekly_grouped_pos[1:(end-1)],lab = "Data",ms = 5)
    return plt
end


function plotfittoPCRdata(areamodel::CoVAreaModel)
        # xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:11 ]
        # xticklabs = vcat([(monthname(k)[1:3])*" 20" for k = 3:12],["Jan 21"])
        xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:10 ]
        xticklabs = [(monthname(k)[1:3]) for k = 3:12]
        inc = incidence_across_samples(areamodel,315);
        PCR_cred = create_PCR_credible_intervals(inc.PCR_incidence_samples)
        PCR_pred = create_PCR_prediction_intervals(inc.PCR_incidence_samples,areamodel.MCMC_results.chain[:,4,1])

        plt_PCRfit = plot(PCR_cred.mean_pred,
                ribbon = (PCR_cred.lb_pred,PCR_cred.ub_pred),
                color = :black,
                lw = 2,
                lab = "Mean number of PCR detectable tested each day (95% CI)",
                ylims = (0.,600),
                xticks = (xticktimes,xticklabs),
                ylabel = "Number of positive swab tests",
                title = "$(areamodel.areaname): Positive swab test prediction and data",
                size = (700,500),dpi = 250)

        plot!(plt_PCRfit, PCR_pred.mean_pred,
                ribbon = (PCR_pred.lb_pred,PCR_pred.ub_pred),
                color = :red,
                fillalpha = 0.1,
                lw = 0,
                lab = "Prediction of swab tests positive each day (95% PI)")

        scatter!(plt_PCRfit,areamodel.PCR_cases,
                color = :blue,
                lw = 0,
                lab = "Data: Swab positive tests by collection date")
        return plt_PCRfit
end


function population_plot(areamodel::KenyaSerology.CoVAreaModel)
    inc = KenyaSerology.incidence_across_samples(areamodel,315);
    sero_pos = KenyaSerology.create_credible_intervals(inc.sero_converted);
    seroweekdata = KenyaSerology.groupserobyweek(areamodel);
    trueinfecteds = KenyaSerology.create_credible_intervals(inc.true_infecteds);
    truesusceptible = KenyaSerology.create_credible_intervals(inc.true_susceptible);

    # xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:11 ]
    # xticklabs = vcat([(monthname(k)[1:3])*" 20" for k = 3:12],["Jan 21"])
    xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:8 ]
    xticklabs = [(monthname(k)[1:3]) for k = 3:10]
    plt = plot(truesusceptible.mean_pred*100/areamodel.N,
            ribbon = (truesusceptible.lb_pred*100/areamodel.N,truesusceptible.ub_pred*100/areamodel.N),
            color = :blue,
            xticks = (xticktimes,xticklabs),
            lab = "Percentage uninfected (95% CI)",
            ylabel = "% of population",
            title = "$(areamodel.areaname): population predictions and data",
            size = (700,500),dpi = 250,
            fillalpha = 0.3)

    plot!(plt,trueinfecteds.mean_pred*100/areamodel.N,
            ribbon = (trueinfecteds.lb_pred*100/areamodel.N,trueinfecteds.ub_pred*100/areamodel.N),
            color = :red,
            xticks = (xticktimes,xticklabs),
            lab = "Percentage exposed to SARS-CoV-2 (95% CI)",
            fillalpha = 0.3)

    plot!(plt,sero_pos.mean_pred*100,
            ribbon = (sero_pos.lb_pred*100,sero_pos.ub_pred*100),
            color = :green,
            xticks = (xticktimes,xticklabs),
            lab = "Prediction: percentage obs. seropositive (95% CI)")

    scatter!(plt,seroweekdata.midpointsofweek,(seroweekdata.posbyweek./seroweekdata.numbyweek)*100,
            ms = 12,
            yerror = (seroweekdata.lpred*100,seroweekdata.upred*100),
            color = :green,
            lab = "Data: weekly seropositive")

    return plt
end

function population_plot(areamodel::KenyaSerology.CoVAreaModel,::Val{:monthlyserology})
    inc = KenyaSerology.incidence_across_samples(areamodel,315);
    sero_pos = KenyaSerology.create_credible_intervals(inc.sero_converted);
    seromonthdata = groupserobymonth(areamodel);

    trueinfecteds = KenyaSerology.create_credible_intervals(inc.true_infecteds);
    truesusceptible = KenyaSerology.create_credible_intervals(inc.true_susceptible);

    # xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:11 ]
    # xticklabs = vcat([(monthname(k)[1:3])*" 20" for k = 3:12],["Jan 21"])
    xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:8 ]
    xticklabs = [(monthname(k)[1:3]) for k = 3:10]
    plt = plot(truesusceptible.mean_pred*100/areamodel.N,
            ribbon = (truesusceptible.lb_pred*100/areamodel.N,truesusceptible.ub_pred*100/areamodel.N),
            color = :blue,
            xticks = (xticktimes,xticklabs),
            lab = "Percentage uninfected (95% CI)",
            ylabel = "% of population",
            title = "$(areamodel.areaname): population predictions and data",
            size = (700,500),dpi = 250,
            fillalpha = 0.3)

    plot!(plt,trueinfecteds.mean_pred*100/areamodel.N,
            ribbon = (trueinfecteds.lb_pred*100/areamodel.N,trueinfecteds.ub_pred*100/areamodel.N),
            color = :red,
            xticks = (xticktimes,xticklabs),
            lab = "Percentage infected (95% CI)",
            fillalpha = 0.3)

    plot!(plt,sero_pos.mean_pred*100,
            ribbon = (sero_pos.lb_pred*100,sero_pos.ub_pred*100),
            color = :green,
            xticks = (xticktimes,xticklabs),
            lab = "Percentage seropositive (95% CI)")

    scatter!(plt,seromonthdata.midpointsofmonth,(seromonthdata.posbymonth./seromonthdata.numbymonth)*100,
            ms = 12,
            yerror = (seromonthdata.lpred*100,seromonthdata.upred*100),
            color = :green,
            lab = "Data: monthly seropositive")

    return plt
end

function population_plot(areamodel::KenyaSerology.CoVAreaModel,::Val{:weeklyaverage})
    inc = KenyaSerology.incidence_across_samples(areamodel,315);
    sero_pos = KenyaSerology.create_credible_intervals(inc.sero_converted_samples);
    seroplotting = seroweekaverage(areamodel)

    trueinfecteds = KenyaSerology.create_credible_intervals(inc.true_infecteds);
    truesusceptible = KenyaSerology.create_credible_intervals(inc.true_susceptible);

    # xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:11 ]
    # xticklabs = vcat([(monthname(k)[1:3])*" 20" for k = 3:12],["Jan 21"])
    xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:10 ]
    xticklabs = [(monthname(k)[1:3]) for k = 3:12]
    plt = plot(truesusceptible.mean_pred*100/areamodel.N,
            ribbon = (truesusceptible.lb_pred*100/areamodel.N,truesusceptible.ub_pred*100/areamodel.N),
            color = :blue,
            xticks = (xticktimes,xticklabs),
            lab = "Percentage uninfected (95% CI)",
            ylabel = "% of population",
            title = "$(areamodel.areaname): population predictions and data",
            size = (700,500),dpi = 250,
            fillalpha = 0.3)

    plot!(plt,trueinfecteds.mean_pred*100/areamodel.N,
            ribbon = (trueinfecteds.lb_pred*100/areamodel.N,trueinfecteds.ub_pred*100/areamodel.N),
            color = :red,
            xticks = (xticktimes,xticklabs),
            lab = "Percentage infected (95% CI)",
            fillalpha = 0.3)

    plot!(plt,sero_pos.mean_pred*100,
            ribbon = (sero_pos.lb_pred*100,sero_pos.ub_pred*100),
            color = :green,
            xticks = (xticktimes,xticklabs),
            lab = "Percentage seropositive (95% CI)")

    scatter!(plt,seroplotting.midpointsofweek,(seroplotting.posbyweek./seroplotting.numbyweek)*100,
            ms = 10,
            yerror = (seroplotting.lpred*100,seroplotting.upred*100),
            color = :green,
            lab = "Data: weekly moving average")

    return plt
end


function plot_incidence(areamodel::CoVAreaModel,areadeaths,p_ID)
        inc = incidence_across_samples(areamodel,315);
        incidence = create_credible_intervals(inc.true_incidence);

        # xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:11 ]
        # xticklabs = vcat([(monthname(k)[1:3])*" 20" for k = 3:12],["Jan 21"])
        xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:10 ]
        xticklabs = [(monthname(k)[1:3]) for k = 3:12]

        f_err = IFR -> death_pred_err(areadeaths,incidence.mean_pred,IFR,p_ID)
        IFR_fit = optimize(f_err,0,0.1)
        IFR_hat = IFR_fit.minimizer

        deaths = create_death_prediction_intervals(death_pred(incidence.mean_pred,IFR_hat,p_ID))

        plt = plot(incidence.mean_pred.+1,
                ribbon = (incidence.lb_pred.+1,incidence.ub_pred.+1),
                yscale = :log10,
                lab = "Infections (95% CI)",
                legend = :topright,
                xticks = (xticktimes,xticklabs),
                color = :red,
                fillalpha = 0.1,
                yticks = ([1,2,11,101,1001,10001],[0,1,10,100,1000,10000]),
                ylims = (0.5,6e4),
                ylabel = "Daily incidence",
                title = "$(areamodel.areaname): incidence of infection and death",
                titlefont = (18,"helvetica"),guidefont = (14,"helvetica"),tickfont = (10,"helvetica"),
                size = (700,500),dpi = 250,
                left_margin = 5mm,
                right_margin = 5mm)

        plot!(cumsum(deaths.mean_pred),
            grid = nothing,
            inset = (1,bbox(-0.1, -0.06, 0.4, 0.2, :center)),
            xticks = (xticktimes[1:end],xticklabs[1:end]),
            subplot = 2,
            lab="",
            color = :black,
            titlefont = (8,"helvetica"),
            title = "Cum. observed deaths",
            bg_inside = nothing)

        scatter!(cumsum(areadeaths),
                color = :black,
                lab = "",
                subplot = 2)

        plot!(deaths.mean_pred.+1,
                ribbon = (deaths.lb_pred,deaths.ub_pred),
                fillalpha = 0.2,
                lab = "Pred. deaths (95% CI; crude IFR = $(round(IFR_hat*100,digits =4))%)",
                color = :black,
                subplot = 1)


        scatter!(areadeaths.+1,
                color = :black,
                lab = "Data: daily deaths",
                subplot = 1)

        return plt
end

function plot_incidenceonly(areamodel::CoVAreaModel)
        sept_30 = (Date(2020,9,30) - Date(2020,2,20)).value

        inc = incidence_across_samples(areamodel,sept_30+1);
        incidence = create_credible_intervals(inc.true_incidence);

        # xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:11 ]
        # xticklabs = vcat([(monthname(k)[1:3])*" 20" for k = 3:12],["Jan 21"])
        xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:8 ]
        xticklabs = [(monthname(k)[1:3]) for k = 3:10]

        plt = plot(incidence.mean_pred,
                ribbon = (incidence.lb_pred,incidence.ub_pred),
                lab = "Infections (95% CI)",
                legend = :topright,
                xticks = (xticktimes,xticklabs),
                color = :red,
                fillalpha = 0.1,
                lw=2,
                # yticks = ([1,2,11,101,1001,10001],[0,1,10,100,1000,10000]),
                # ylims = (0.5,6e4),
                ylabel = "Daily incidence",
                title = "$(areamodel.areaname): incidence of infection",
                titlefontsize = 18,guidefontsize = 20,tickfontsize = 14,
                legendfontsize = 12,
                size = (700,500),dpi = 250,
                left_margin = 5mm,
                right_margin = 5mm)


        return plt
end

function plot_deaths(areamodel::CoVAreaModel,areadeaths,p_ID)
        sept_30 = (Date(2020,9,30) - Date(2020,2,20)).value

        inc = incidence_across_samples(areamodel,sept_30+1);
        incidence = create_credible_intervals(inc.true_incidence);

        # xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:11 ]
        # xticklabs = vcat([(monthname(k)[1:3])*" 20" for k = 3:12],["Jan 21"])
        xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:8 ]
        xticklabs = [(monthname(k)[1:3]) for k = 3:10]

        D = sum(areadeaths[1:165])
        μ_ests = [Erlang(D+1,1/(0.0015 + sum(KenyaSerology.simple_conv(inc.true_incidence[:,n],p_ID)[1:165]))) for n = 1: 10000]
        IFR_hat = mean(mean.(μ_ests))
        deaths = create_death_prediction_intervals(death_pred(incidence.mean_pred,IFR_hat,p_ID))

        plt = plot(deaths.mean_pred,
                ribbon = (deaths.lb_pred,deaths.ub_pred),
                fillalpha = 0.2,
                lw = 3,
                lab = "Pred. deaths (95% CI)",
                color = :black,
                subplot = 1,
                legend = :topright,
                ylabel = "Daily deaths",
                xticks = (xticktimes,xticklabs),
                title = "$(areamodel.areaname): observed and predicted deaths",
                titlefontsize = 18,guidefontsize = 20,tickfontsize = 14,
                legendfontsize = 12,
                size = (700,500),dpi = 250,
                left_margin = 5mm,
                right_margin = 5mm)

        plot!(cumsum(deaths.mean_pred),
            grid = nothing,
            inset = (1,bbox(-0.30, -0.3, 0.3, 0.2, :center)),
            xticks = (xticktimes[1:2:end],xticklabs[1:2:end]),
            subplot = 2,
            lab="",
            color = :black,
            titlefont = 10,
            title = "Cumulative deaths",
            bg_inside = nothing)

        scatter!(cumsum(areadeaths),
                color = :black,
                lab = "",
                subplot = 2)

        scatter!(areadeaths,
                color = :black,
                ms = 10,
                lab = "Data: daily deaths",
                subplot = 1)

        return plt
end

function plot_deaths_weekly(areamodel::CoVAreaModel,areadeaths,p_ID,mainylim)
        sept_30 = (Date(2020,9,30) - Date(2020,2,20)).value
        IFR_priormean = 0.00264
        inc = incidence_across_samples(areamodel,sept_30+1);
        incidence = create_credible_intervals(inc.true_incidence);

        # xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value for k = 1:11 ]
        # xticklabs = vcat([(monthname(k)[1:3])*" 20" for k = 3:12],["Jan 21"])
        xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,21)).value  for k = 1:8 ]
        xticklabs = [(monthname(k)[1:3]) for k = 3:10]

        D = sum(areadeaths)
        μ_ests = [Erlang(D+1,1/((1/IFR_priormean)+ sum(KenyaSerology.simple_conv(inc.true_incidence[:,n],p_ID)[1:length(areadeaths)]))) for n = 1: 10000]
        IFR_hat = mean(mean.(μ_ests))

        weekly_grouped_deaths = KenyaSerology.group_PCR_pos_by_week(death_pred(incidence.mean_pred,IFR_hat,p_ID))
        weekly_grouped_death_data = KenyaSerology.group_PCR_pos_by_week(areadeaths)
        deaths = create_death_prediction_intervals(weekly_grouped_deaths)

        #First week is week 8 of year which starts 17th Feb i.e. -4 on the tick scale
        weekdates_on_scale = collect(0:(length(deaths.mean_pred)-1)).*7 .- 4

        plt = plot(weekdates_on_scale[1:(end-1)],deaths.mean_pred[1:(end-1)],
                ribbon = (deaths.lb_pred[1:(end-1)],deaths.ub_pred[1:(end-1)]),
                fillalpha = 0.2,
                lw = 3,
                lab = "Pred. deaths (95% CI)",
                color = :black,
                subplot = 1,
                legend = :topright,
                ylims = (-2,mainylim),
                ylabel = "Weekly deaths",
                xticks = (xticktimes,xticklabs),
                title = "$(areamodel.areaname): observed and predicted deaths",
                titlefontsize = 18,guidefontsize = 20,tickfontsize = 14,
                legendfontsize = 12,
                size = (700,500),dpi = 250,
                left_margin = 5mm,
                right_margin = 5mm)

        plot!(weekdates_on_scale[1:(end-1)],cumsum(deaths.mean_pred[1:(end-1)]),
            grid = nothing,
            inset = (1,bbox(-0.30, -0.3, 0.3, 0.2, :center)),
            xticks = (xticktimes[1:2:end],xticklabs[1:2:end]),
            subplot = 2,
            lab="",
            color = :black,
            titlefont = 10,
            title = "Cumulative deaths",
            bg_inside = nothing)

        scatter!(weekdates_on_scale[1:(end-1)],cumsum(weekly_grouped_death_data[1:(end-1)]),
                color = :black,
                lab = "",
                subplot = 2)

        scatter!(weekdates_on_scale[1:(end-1)],weekly_grouped_death_data[1:(end-1)],
                color = :black,
                ms = 10,
                lab = "Data: weekly deaths",
                subplot = 1)

        return plt
end
