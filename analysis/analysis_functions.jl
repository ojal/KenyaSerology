function generate_model_for_fit(areaname,popsize,logprior,::Val{:var_testing})
        lwr_name = String([lowercase(s) for s in areaname])
        upr_name = String([uppercase(s) for s in areaname])
        cases = vec(case_data.cases[:,case_data.areas .== upr_name])
        sero = zeros(size(matchedserodata[:,1,:]))
        if sum(sero_data.areas.==lwr_name) > 0
           sero =  Matrix(matchedserodata[:,sero_data.areas.==lwr_name,:][:,1,:])
        end

        return KenyaSerology.CoVAreaModel(areaname = areaname,
                                        PCR_cases = cases,
                                        sero_cases = sero,
                                        dates = case_data.dates,
                                        N = popsize,
                                        γ = 1/(5.5 - 3.1),
                                        contactrate_data = projected_contactrate_kenya,
                                        relative_testing_rate = relative_testing_rate_nairobi.relative_testing_rate, #This is where relative testing rate based on relative to mean testing in Kenya goes in
                                        prob = KenyaSerology.make_odeproblemforinference(projected_contactrate_nairobi,
                                                                                                startdate = Date(2020,2,21),
                                                                                                enddate = Date(2020,8,6)),
                                        sero_array = rel_sero_array_26days,
                                        log_priors = logprior,
                                        log_likelihood = KenyaSerology.loglikelihood_contactratemodelBB_Peff),vec(death_data.deaths[:,death_data.areas .== upr_name])


end


function generate_model_for_fit(areaname,popsize,logprior)
        lwr_name = String([lowercase(s) for s in areaname])
        upr_name = String([uppercase(s) for s in areaname])
        cases = vec(case_data.cases[:,case_data.areas .== upr_name])
        sero = zeros(size(matchedserodata[:,1,:]))
        if sum(sero_data.areas.==lwr_name) > 0
           sero =  Matrix(matchedserodata[:,sero_data.areas.==lwr_name,:][:,1,:])
        end

        return KenyaSerology.CoVAreaModel(areaname = areaname,
                                        PCR_cases = cases,
                                        sero_cases = sero,
                                        dates = case_data.dates,
                                        N = popsize,
                                        γ = 1/(5.5 - 3.1),
                                        contactrate_data = projected_contactrate_kenya,
                                        prob = KenyaSerology.make_odeproblemforinference(projected_contactrate_kenya,
                                                                                                startdate = Date(2020,2,21),
                                                                                                enddate = Date(2020,8,6)),
                                        sero_array = rel_sero_array_26days,
                                        log_priors = logprior,
                                        log_likelihood = KenyaSerology.loglikelihood_contactratemodelBB_Peff),vec(death_data.deaths[:,death_data.areas .== upr_name])

end

function recordpeakinfectionsanddeaths(areamodel::KenyaSerology.CoVAreaModel,areadeaths,IFR_est,p_ID)
   inc = KenyaSerology.incidence_across_samples(areamodel,315);
   incidence = KenyaSerology.create_credible_intervals(inc.true_incidence);
   time_length = size(inc.PCR_incidence_samples,1)
   mean_pred_PCR = zeros(time_length)
   std_pred_PCR = zeros(time_length)
   mean_pred_incidence = zeros(time_length)
   std_pred_incidence = zeros(time_length)
   for t = 1:time_length
        PCRinconthatday = inc.PCR_incidence_samples[t,:]
        inconthatday = inc.true_incidence[t,:]
        mean_pred_PCR[t] = mean(PCRinconthatday)
        std_pred_PCR[t] = std(PCRinconthatday)
        mean_pred_incidence[t] = mean(inconthatday)
        std_pred_incidence[t] = std(inconthatday)
    end

   f_err = IFR -> KenyaSerology.death_pred_err(areadeaths,incidence.mean_pred,IFR,p_ID)
   IFR_fit = optimize(f_err,0,0.1)
   IFR_hat = IFR_fit.minimizer
   deaths = KenyaSerology.create_death_prediction_intervals(KenyaSerology.death_pred(incidence.mean_pred,IFR_hat,p_ID))
   deaths_estIFR = KenyaSerology.create_death_prediction_intervals(KenyaSerology.death_pred(incidence.mean_pred,IFR_est,p_ID))

   post_mean_peak = mean([findmax(inc.true_incidence[:,k])[2] for k = 1:size(inc.true_incidence,2)])

   mean_R = mean(areamodel.MCMC_results.chain[:,1,1])

   return (area = areamodel.areaname,
                peak = post_mean_peak,
                pred_deaths_ftc = deaths.mean_pred,
                pred_deaths_est = deaths_estIFR.mean_pred,
                mean_pred_PCR = mean_pred_PCR,
                std_pred_PCR = std_pred_PCR,
                mean_pred_incidence = mean_pred_incidence,
                std_pred_incidence = std_pred_incidence,
                mean_R = mean_R)
end

function pipeline_for_fit(areaname,popsize,logprior,::Val{:var_testing})
   model,deaths = generate_model_for_fit(areaname,popsize,logprior,Val(:var_testing))
   KenyaSerology.inferparameters!(model,10000,trans_Peff)
   return model,deaths
end

function pipeline_for_fit(areaname,popsize,logprior)
   model,deaths = generate_model_for_fit(areaname,popsize,logprior)
   KenyaSerology.inferparameters!(model,10000,trans_Peff)
   return model,deaths
end


function gatherdatafrommodels(dirname,death_data,IFR_est,p_ID)
    modelnames = readdir(dirname)
    collecteddata = Any[]
    for name in modelnames
        modelpath = joinpath(dirname,name)
        D = FileIO.load(modelpath)
        f = first(keys(D))
        model = D[f]
        println("Gathering data for $(model.areaname)")
        name_upr = String([uppercase(s) for s in model.areaname])
        deaths = death_data.deaths[:,death_data.areas .== name_upr]
        datafrommodel = recordpeakinfectionsanddeaths(model,deaths,IFR_est,p_ID)
        push!(collecteddata,datafrommodel)
    end
    return collecteddata
end

function createpeaktimes(collecteddata,countynames)
    peaktimesbycounty = (-1)*ones(length(countynames))
    for data in collecteddata
        peaktimesbycounty[uppercase.(data.area) .== uppercase.(countynames)] .= data.peak
    end
    return peaktimesbycounty
end

function create_univariate_distributions(post)
    d_R = fit_mle(Gamma,post[:,1,1])
    d_E₀ = fit_mle(Gamma,post[:,2,1])
    d_I₀ = fit_mle(Gamma,post[:,3,1])
    d_α = fit_mle(Gamma,post[:,4,1])
    d_ptest = fit_mle(Gamma,post[:,5,1])
    #MLE fit for the beta distribution
    l(x) = sum(map( p -> -logpdf(Beta(x[1],x[2]),p),post[:,6,1]))
    beta_mle = optimize(l,[2.,1.])
    d_Peff = Beta(beta_mle.minimizer[1],beta_mle.minimizer[2])
    return (d_R = d_R,
            d_E₀= d_E₀,
            d_I₀ = d_I₀,
            d_α = d_α,
            d_ptest = d_ptest,
            d_Peff=d_Peff)
end

function gather_posteriors(modelslist)
    D = load(modelslist[1])
    model = D[first(keys(D))]
    post = deepcopy(model.MCMC_results.chain)
    for modelnames in modelslist[2:end]
        D = load(modelnames)
        model = D[first(keys(D))]
        post = vcat(post,model.MCMC_results.chain[:,:,1])
    end
    return post
end

function gather_and_create_posterior_fits(models_semiurban)
        post = gather_posteriors(models_semiurban)
        return create_univariate_distributions(post)
end


function combinedata(combindeddata)
        kenyaincidence = zeros(size(combindeddata[1].mean_pred_incidence))
        kenyaincidence_std = zeros(size(combindeddata[1].mean_pred_incidence))
        kenyaPCR = zeros(size(combindeddata[1].mean_pred_PCR))
        kenyaPCR_std = zeros(size(combindeddata[1].mean_pred_PCR))
        kenyadeaths = zeros(size(combindeddata[1].pred_deaths_ftc))

        for data in combindeddata
                kenyaincidence .+= data.mean_pred_incidence
                kenyaincidence_std .+= data.std_pred_incidence.^2
                kenyaPCR .+= data.mean_pred_PCR
                kenyaPCR_std .+= data.std_pred_PCR.^2
                kenyadeaths .+= data.pred_deaths_ftc
        end

        kenyaincidence_std = sqrt.(kenyaincidence_std)
        kenyaPCR_std = sqrt.(kenyaPCR_std)
        kenyadeaths_upred = [invlogcdf(Poisson(d),log(0.975)) for d in kenyadeaths]
        kenyadeaths_lpred = [invlogcdf(Poisson(d),log(0.025)) for d in kenyadeaths]

        return (kenyaincidence=kenyaincidence,
                kenyaincidence_std=kenyaincidence_std,
                kenyaPCR = kenyaPCR,
                kenyaPCR_std = kenyaPCR_std,
                kenyadeaths = kenyadeaths,
                kenyadeaths_upred = kenyadeaths_upred .- kenyadeaths,
                kenyadeaths_lpred = kenyadeaths .- kenyadeaths_lpred)
end

function death_mean_pred(obs_deaths,incidence,IFR_array,p_ID)
    unscaleddeaths = [KenyaSerology.death_pred(incidence[:,j],1.,p_ID) for j in 1:1000]
    err = 0.
    for (i,true_deaths) in enumerate(obs_deaths)
        err -= log(mean([pdf(Poisson(max(IFR*unscaleddeaths[j][i],0.)),true_deaths) for (j,IFR) in enumerate(IFR_array[1:1000])]))
    end
    return err
end

function datafortable(areamodel,deaths)
    inc = KenyaSerology.incidence_across_samples(areamodel,315);
    incidence = KenyaSerology.create_credible_intervals(inc.true_incidence);
    total_inf = KenyaSerology.create_credible_intervals(inc.true_infecteds);
    mean_infection_density = total_inf.mean_pred[end]*100/areamodel.N
    infection_density_CI = (mean_infection_density - total_inf.lb_pred[end]*100/areamodel.N,mean_infection_density+total_inf.ub_pred[end]*100/areamodel.N)
    mean_R = mean(areamodel.MCMC_results.chain[:,1,1])
    R_CI = (quantile(areamodel.MCMC_results.chain[:,1,1],0.025),quantile(areamodel.MCMC_results.chain[:,1,1],0.975))
    mean_Peff = mean(areamodel.MCMC_results.chain[:,6,1])
    Peff_CI = (quantile(areamodel.MCMC_results.chain[:,6,1],0.025),quantile(areamodel.MCMC_results.chain[:,6,1],0.975))

    #Fit deaths over first 165 days using simple Bayesian method.
    D = sum(deaths[1:165])
    μ_ests = [Erlang(D+1,1/((1/0.00264) + sum(KenyaSerology.simple_conv(inc.true_incidence[:,n],p_ID)[1:165]))) for n = 1: 10000]
    IFR_hat = mean(mean.(μ_ests))

    log_pred_den = death_mean_pred(deaths,inc.true_incidence,mean.(μ_ests),p_ID)

    return (lpd = log_pred_den,
            IFR_hat = IFR_hat,
            IFR_CI = (quantile(mean.(μ_ests),0.025), quantile(mean.(μ_ests),0.975) ),
            mean_infection_density =mean_infection_density,
            infection_density_CI = infection_density_CI,
            mean_R = mean_R,
            R_CI = R_CI,
            mean_Peff = mean_Peff,
            Peff_CI=Peff_CI
            )
end

function createmeanandCIforeachparameter(model)
	posteriormeans = [mean(model.MCMC_results.chain[:,k,1]) for k = 1:size(model.MCMC_results.chain,2)]
	CIs = [(quantile(model.MCMC_results.chain[:,k,1],0.025), quantile(model.MCMC_results.chain[:,k,1],0.975)) for k = 1:size(model.MCMC_results.chain,2)]
	return (posteriormeans=posteriormeans,CIs=CIs)
end
function getpeakmeanandCI(inc)
	peaktimes = [findmax(inc.true_incidence[:,k])[2] for k = 1:size(inc.true_incidence,2)]
	return (peakmean = mean(peaktimes),
			peakCI = (quantile(peaktimes,0.025), quantile(peaktimes,0.975)))
end
function getIFRmeanandCI(inc,deaths,IFR_priormean,p_ID)
	D = sum(deaths[1:165])
	IFR_ests = mean.([Erlang(D+1,1/((1/IFR_priormean) + sum(KenyaSerology.simple_conv(inc.true_incidence[:,n],p_ID)[1:165]))) for n = 1: 10000])
	return (IFRmean = mean(IFR_ests),
			IFRCI = (quantile(IFR_ests,0.025), quantile(IFR_ests,0.975)))
end

function generate_simulated_death_lpds(model,IFR_priormean,deaths,p_ID)
    inc = KenyaSerology.incidence_across_samples(model,315)
    D= sum(deaths[1:165])
    IFR_array = mean.([Erlang(D+1,1/((1/IFR_priormean) + sum(KenyaSerology.simple_conv(inc.true_incidence[:,n],p_ID)[1:165]))) for n = 1:1000])
    lpd_array = zeros(1000)
    rand_deaths = zeros(Int64,165)
    lpd_actual = death_mean_pred(deaths,inc.true_incidence,IFR_array,p_ID)
    for k = 1:1000
        rand_deaths = [rand(Poisson(max(IFR_array[j]*λ,0.))) for (j,λ) in enumerate(KenyaSerology.simple_conv(inc.true_incidence[:,k],p_ID)[1:165]) ]
        lpd_array[k] = death_mean_pred(rand_deaths,inc.true_incidence,IFR_array,p_ID)
    end
    return lpd_array,lpd_actual
end


function getdataforCSVfile(model,deaths,IFR_priormean,p_ID)
	paramfits = createmeanandCIforeachparameter(model)
	inc = KenyaSerology.incidence_across_samples(model,315)
	peakfits = getpeakmeanandCI(inc)
	IFRfits = getIFRmeanandCI(inc,deaths,IFR_priormean,p_ID)
    posterior_predictive_p= -1.
    try
        lpd_array,lpd_actual = generate_simulated_death_lpds(model,IFR_priormean,deaths,p_ID)
        posterior_predictive_p = sum(lpd_actual .< lpd_array)/1000
    catch
        posterior_predictive_p = -1.
    end
	return (paramfits=paramfits,
			peakfits=peakfits,
			IFRfits=IFRfits,
            posterior_predictive_p=posterior_predictive_p)
end

function createdatarow(model,deaths,IFR_prior,p_ID)
	name = model.areaname
	data = getdataforCSVfile(model,deaths,IFR_prior,p_ID)
	paramfit_str = [string(round(data.paramfits.posteriormeans[k],sigdigits = 3))*" "*string(round.(data.paramfits.CIs[k],sigdigits = 3)) for k = 1:6]
	meanpeakdate = Date(2020,2,20) + Day(round(Int64,data.peakfits.peakmean))
	peakCIdates = (Date(2020,2,20) + Day(round(Int64,data.peakfits.peakCI[1])), Date(2020,2,20) + Day(round(Int64,data.peakfits.peakCI[2]))  )
	peak_str = string(meanpeakdate)*" "*string(peakCIdates)
	IFRfit_str = string(round(data.IFRfits.IFRmean*100,sigdigits = 3))*"% ("*string(round(data.IFRfits.IFRCI[1]*100,sigdigits = 3))*"%, "*string(round(data.IFRfits.IFRCI[2]*100,sigdigits = 3))*"%)"
    posterior_predictive_p = string(data.posterior_predictive_p)
	return vcat(name,paramfit_str,peak_str,IFRfit_str,posterior_predictive_p)
end

function parameterinferenceovercollection(dirname,death_data,IFR_prior,p_ID)
    modelnames = readdir(dirname)
    df = DataFrame(Countyname = String[],
					R0 = String[],
					E0 = String[],
					I0 = String[],
					clusteringfactor = String[],
					p_test = String[],
					EffPopSize = String[],
					dateinfectionpeak = String[],
					IFR = String[],
                    PosteriorpredicitveP = String[])
    for name in modelnames
        modelpath = joinpath(dirname,name)
        D = FileIO.load(modelpath)
        f = first(keys(D))
        model = D[f]
        println("Gathering data for $(model.areaname)")
        name_upr = String([uppercase(s) for s in model.areaname])
        deaths = death_data.deaths[:,death_data.areas .== name_upr]
		row = createdatarow(model,deaths,0.00264,p_ID)
        push!(df,row)
    end
    return df
end
