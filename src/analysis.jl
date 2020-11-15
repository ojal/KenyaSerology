
"""
    generate_model_for_fit(areaname,popsize,logprior,case_data,sero_data,contactrate_data)

Create a CoVAreaModel model from testing, case, and serological data along with population level information e.g. population size.
"""
function generate_model_for_fit(areaname,popsize,logprior,case_data,sero_data,contactrate_data,relative_testing_rate,enddate::Date;
									loglikelihood_func = loglikelihood_contactratemodelBB_Peff)
        lwr_name = String([lowercase(s) for s in areaname])
        upr_name = String([uppercase(s) for s in areaname])
		cases = Array{Float64,ndims(case_data.cases)-1}[]
		if ndims(case_data.cases) == 2
	        cases = vec(case_data.cases[:,case_data.areas .== upr_name])
		end
		if ndims(case_data.cases) == 3
			cases = Matrix(case_data.cases[:,case_data.areas .== upr_name,:][:,1,:])
		end

        sero = zeros(size(case_data.cases,1),2)
        if sum(sero_data.areas.==lwr_name) > 0
           sero =  Matrix(sero_data.serodata[:,sero_data.areas.==lwr_name,:][:,1,:])
        end

        return CoVAreaModel(areaname = areaname,
                                        PCR_cases = cases,
                                        sero_cases = sero,
                                        dates = case_data.dates,
                                        N = popsize,
                                        γ = 1/(5.5 - 3.1),
                                        contactrate_data = contactrate_data,
                                        relative_testing_rate = relative_testing_rate,
                                        prob = make_odeproblemforinference(contactrate_data,
                                                                            startdate = Date(2020,2,21),
                                                                            enddate = enddate),
                                        log_priors = logprior,
                                        log_likelihood = loglikelihood_func)


end

"""
    pipeline_for_fit(areaname,popsize,logprior,case_data,sero_data,contactrate_data)

Create a CoVAreaModel model, populated with relevant data, and fit the model by drawing n MCMC samples.
"""
function pipeline_for_fit(areaname,popsize,logprior,case_data,sero_data,contactrate_data,relative_testing_rate,enddate;
                        n,parameter_transformation,loglikelihood_func = loglikelihood_contactratemodelBB_Peff)
   model = generate_model_for_fit(areaname,popsize,logprior,case_data,sero_data,contactrate_data,relative_testing_rate,enddate;
   									loglikelihood_func = loglikelihood_func)
   inferparameters!(model,n,parameter_transformation)
   #Add a save fit bit
   return model
end

"""
    fit_univariate_distributions(post)

Fit a collection of univariate distributions to the samples from the posterior distribution post using MLE.
"""
function fit_univariate_distributions(post)
    d_R = fit_mle(Gamma,post[:,1,1])
    d_E₀ = fit_mle(Gamma,post[:,2,1])
    d_I₀ = fit_mle(Gamma,post[:,3,1])
    d_α = fit_mle(Gamma,post[:,4,1])
    d_ptest = fit_mle(Gamma,post[:,5,1])
    #MLE fit for the beta distribution
    l(x) = try
            sum(map( p -> -logpdf(Beta(x[1],x[2]),p),post[:,6,1]))
            catch
            Inf
        end
    beta_mle = optimize(l,[2.,1.])
    d_Peff = Beta(beta_mle.minimizer[1],beta_mle.minimizer[2])
    return (d_R = d_R,
            d_E₀= d_E₀,
            d_I₀ = d_I₀,
            d_α = d_α,
            d_ptest = d_ptest,
            d_Peff=d_Peff)
end

"""
    gather_posteriors(modelslist)

Gather groups of samples from posterior distributions.
"""
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

"""
    gather_and_create_posterior_fits(models_semiurban)

Gather groups of samples from posterior distributions and then fit univariate distributions to the grouped samples.
"""
function gather_and_create_posterior_fits(models_semiurban)
        post = gather_posteriors(models_semiurban)
        return fit_univariate_distributions(post)
end

"""
    combinedata(combindeddata)

Combine an array of NamedTuples containing predictions about the incidence, swab testing, deaths, and hospitalisations in each area into a combined prediction over all the areas.
"""
function combinedata(combindeddata)
        kenyaincidence = zeros(size(combindeddata[1].mean_pred_incidence))
        kenyaincidence_std = zeros(size(combindeddata[1].mean_pred_incidence))
        kenyaPCR = zeros(size(combindeddata[1].mean_pred_PCR))
        kenyaPCR_std = zeros(size(combindeddata[1].mean_pred_PCR))
        kenyadeaths = zeros(size(combindeddata[1].pred_deaths_ftc))
        kenyahosps = zeros(size(combindeddata[1].pred_hosps_ftc))

        for data in combindeddata
                kenyaincidence .+= data.mean_pred_incidence
                kenyaincidence_std .+= data.std_pred_incidence.^2
                kenyaPCR .+= data.mean_pred_PCR
                kenyaPCR_std .+= data.std_pred_PCR.^2
                kenyadeaths .+= data.pred_deaths_ftc
                kenyahosps .+= data.pred_hosps_ftc
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

"""
    logdensity_deaths(obs_deaths,incidence,IFR_array,p_ID)

Estimate the negative log-density for each day of observed deaths, and sum over all days of death data.
"""
function logdensity_deaths(obs_deaths,IFR_array,unscaleddeaths::Matrix{Float64})
    # unscaleddeaths = [death_pred(incidence[:,j],1.,p_ID) for j in 1:size(incidence,2)]
	# rate = simple_conv(incidence[:,1],p_ID)
	# unscaleddeaths = zeros(length(rate),size(incidence,2))
	# for j = 1:size(unscaleddeaths,2)
	# 	unscaleddeaths[:,j] .= simple_conv(incidence[:,j],p_ID)
	# end
    err = 0.
    for (i,true_deaths) in enumerate(obs_deaths)
        err -= log(mean([pdf(Poisson(max(IFR*unscaleddeaths[i,j],0.)),true_deaths) for (j,IFR) in enumerate(IFR_array)]))
    end
    return err
end

"""
    datafortable(areamodel,deaths)

Create a NamedTuple containing information extracted from the posterior distribution samples for the model fit.
"""
function datafortable(areamodel,deaths)
    inc = incidence_across_samples(areamodel,315);
    incidence = create_credible_intervals(inc.true_incidence);
    total_inf = create_credible_intervals(inc.true_infecteds);
    mean_infection_density = total_inf.mean_pred[end]*100/areamodel.N
    infection_density_CI = (mean_infection_density - total_inf.lb_pred[end]*100/areamodel.N,mean_infection_density+total_inf.ub_pred[end]*100/areamodel.N)
    mean_R = mean(areamodel.MCMC_results.chain[:,1,1])
    R_CI = (quantile(areamodel.MCMC_results.chain[:,1,1],0.025),quantile(areamodel.MCMC_results.chain[:,1,1],0.975))
    mean_Peff = mean(areamodel.MCMC_results.chain[:,6,1])
    Peff_CI = (quantile(areamodel.MCMC_results.chain[:,6,1],0.025),quantile(areamodel.MCMC_results.chain[:,6,1],0.975))

	x = simple_conv(inc.true_incidence[:,1],p_ID)
	unscaleddeaths = zeros(length(x),size(inc.true_incidence,2))
	for j = 1:size(unscaleddeaths,2)
		unscaleddeaths[:,j] .= simple_conv(inc.true_incidence[:,j],p_ID)
	end

    D = sum(deaths)
    μ_ests = [Erlang(D+1,1/((1/0.00264) + sum(simple_conv(inc.true_incidence[:,n],p_ID)[1:(end-3)]))) for n = 1:size(inc.true_incidence,2)]
    IFR_hat = mean(mean.(μ_ests))

    log_pred_den = logdensity_deaths(deaths,mean.(μ_ests),unscaleddeaths)

    return (lpd = log_pred_den,
            IFR_hat = IFR_hat,
            IFR_CI = (quantile(mean.(μ_ests),0.025), quantile(mean.(μ_ests),0.975) ),
            mean_infection_density =mean_infection_density,
            infection_density_CI = infection_density_CI,
            mean_R = mean_R,
            R_CI = R_CI,
            mean_Peff = mean_Peff,
            Peff_CI=Peff_CI)
end

function createmeanandCIforeachparameter(model)
	paramnames = string.(keys(model.MCMC_results.chain))
	posteriormeans = [mean(model.MCMC_results.chain[:,k,1]) for k = 1:size(model.MCMC_results.chain,2)]
	CIs = [(quantile(model.MCMC_results.chain[:,k,1],0.025), quantile(model.MCMC_results.chain[:,k,1],0.975)) for k = 1:size(model.MCMC_results.chain,2)]
	return (posteriormeans=posteriormeans,CIs=CIs,paramnames=paramnames)
end
function getpeakmeanandCI(inc)
	peaktimes = [findmax(inc.true_incidence[:,k])[2] for k = 1:size(inc.true_incidence,2)]
	return (peakmean = mean(peaktimes),
			peakCI = (quantile(peaktimes,0.025), quantile(peaktimes,0.975)))
end
function getIFRmeanandCI(inc,deaths,IFR_priormean,p_ID)
	D = sum(deaths)
	IFR_ests = mean.([Erlang(D+1,1/((1/IFR_priormean) + sum(simple_conv(inc.true_incidence[:,n],p_ID)[1:length(deaths)]))) for n = 1:size(inc.true_incidence,2)])
	return (IFRmean = mean(IFR_ests),
			IFRCI = (quantile(IFR_ests,0.025), quantile(IFR_ests,0.975)))
end

function getIHRdirectmeanandCI(inc,hosps,IHR_priormean,p_IH)
	H = sum(hosps[1:143]) #Cut off at 12th July (day 143) because of poor data
	IHR_ests = mean.([Erlang(H+1,1/((1/IHR_priormean) + sum(simple_conv(inc.true_incidence[:,n],p_IH)[1:143]))) for n = 1: 10000])
	return (IHRmean = mean(IHR_ests),
			IHRCI = (quantile(IHR_ests,0.025), quantile(IHR_ests,0.975)))
end

"""
    generate_simulated_death_lpds(inc,IFR_priormean,deaths,p_ID;num_sims = 10000)

Generate log-predictive densities (LPDs) for 1000 simulated observed death times series, generated from sampled posterior distributions.

The idea is to compare the distribution of LPDs, conditional on the model be correct, of data simulated using the posterior distribution after data fitting to the LPD of the actual observed deaths.
"""
function generate_simulated_death_lpds(inc,IFR_priormean,deaths,p_ID,num_sims)
    # inc = incidence_across_samples(model,315)
    D= sum(deaths)
    IFR_array = mean.([Erlang(D+1,1/((1/IFR_priormean) + sum(simple_conv(inc.true_incidence[:,n],p_ID)[1:length(deaths)]))) for n = 1:size(inc.true_incidence,2)])
    lpd_array = zeros(num_sims)
    generated_deaths = zeros(Int64,length(deaths))
	x = simple_conv(inc.true_incidence[:,1],p_ID)
	unscaleddeaths = zeros(length(x),size(inc.true_incidence,2))
	for j = 1:size(unscaleddeaths,2)
		unscaleddeaths[:,j] .= simple_conv(inc.true_incidence[:,j],p_ID)
	end

    lpd_actual = logdensity_deaths(deaths,IFR_array,unscaleddeaths)
    for k = 1:num_sims
        generated_deaths .= [rand(Poisson(max(IFR_array[j]*λ,0.))) for (j,λ) in enumerate(simple_conv(inc.true_incidence[:,k],p_ID)[1:length(deaths)]) ]
        lpd_array[k] = logdensity_deaths(generated_deaths,IFR_array,unscaleddeaths)
    end
    return lpd_array,lpd_actual
end

"""
    getdataforCSVfile(model,deaths,IFR_priormean,p_ID)

Extract information on posterior mean of parameters, timing of peak infection rate, IHR, IFR and the posterior predictive P-value from a fitted CoVAreaModel model.
"""
function getdataforCSVfile(model,deaths,IFR_priormean,p_ID,num_sims)
	paramfits = createmeanandCIforeachparameter(model)
	inc = incidence_across_samples(model,315)
	# peakfits = getpeakmeanandCI(inc)
    IFRfits = getIFRmeanandCI(inc,deaths,IFR_priormean,p_ID)
    # IHRfits = getIHRdirectmeanandCI(inc,hosps,IFR_priormean*2,p_IH)
    posterior_predictive_p= -1.
    try
        lpd_array,lpd_actual = generate_simulated_death_lpds(inc,IFR_priormean,deaths,p_ID,num_sims)
        posterior_predictive_p = sum(lpd_actual .< lpd_array)/length(lpd_array)
    catch
        posterior_predictive_p = -1.
    end
	return (paramfits=paramfits,
			IFRfits=IFRfits,
            posterior_predictive_p=posterior_predictive_p)
end

function createdatarow(model,deaths,IFR_prior,p_ID,num_sims)
	name = model.areaname
	data = getdataforCSVfile(model,deaths,IFR_prior,p_ID,num_sims)
	paramfit_str = [string(round(data.paramfits.posteriormeans[k],sigdigits = 3))*" "*string(round.(data.paramfits.CIs[k],sigdigits = 3)) for k = 1:length(data.paramfits.posteriormeans)]
	# meanpeakdate = Date(2020,2,20) + Day(round(Int64,data.peakfits.peakmean))
	# peakCIdates = (Date(2020,2,20) + Day(round(Int64,data.peakfits.peakCI[1])), Date(2020,2,20) + Day(round(Int64,data.peakfits.peakCI[2]))  )
    # peak_str = string(meanpeakdate)*" "*string(peakCIdates)
    # IHRfit_str = string(round(data.IHRfits.IHRmean*100,sigdigits = 3))*"% ("*string(round(data.IHRfits.IHRCI[1]*100,sigdigits = 3))*"%, "*string(round(data.IHRfits.IHRCI[2]*100,sigdigits = 3))*"%)"
	IFRfit_str = string(round(data.IFRfits.IFRmean*100,sigdigits = 3))*"% ("*string(round(data.IFRfits.IFRCI[1]*100,sigdigits = 3))*"%, "*string(round(data.IFRfits.IFRCI[2]*100,sigdigits = 3))*"%)"
    posterior_predictive_p = string(data.posterior_predictive_p)
	return vcat(name,paramfit_str,IFRfit_str,posterior_predictive_p)
end

#
# """
#     parameterinferenceovercollection(dirname,death_data,IFR_prior,p_ID)
#
# Apply getdataforCSVfile to each model fit saved in directory at dirname. Output extracted information as a string and pass to a DataFrame
# """
# function parameterinferenceovercollection(dirname,death_data,IFR_prior,p_ID)
#     modelnames = readdir(dirname)
#     df = DataFrame(Countyname = String[],
# 					R0 = String[],
# 					E0 = String[],
# 					I0 = String[],
# 					clusteringfactor = String[],
# 					p_test = String[],
# 					EffPopSize = String[],
#                     dateinfectionpeak = String[],
#                     IHR = String[],
# 					IFR = String[],
#                     PosteriorpredicitveP = String[])
#     for name in modelnames
#         modelpath = joinpath(dirname,name)
#         D = FileIO.load(modelpath)
#         f = first(keys(D))
#         model = D[f]
#         println("Gathering data for $(model.areaname)")
#         name_upr = String([uppercase(s) for s in model.areaname])
#         deaths = death_data.deaths[:,death_data.areas .== name_upr]
# 		row = createdatarow(model,deaths,IFR_prior,p_ID)
#         push!(df,row)
#     end
#     return df
# end


"""
    parameterinferenceovercollection(dirname,death_data,IFR_prior,p_ID)

Apply getdataforCSVfile to each model fit saved in directory at dirname. Output extracted information as a string and pass to a DataFrame
"""
function parameterinferenceovercollection(modelpaths,death_data,IFR_prior,p_ID;num_sims = 1000,enddate::Date = Date(2020,9,30))
	enddate_int = (enddate - Date(2020,2,20)).value
    df = DataFrame(Countyname = String[],
					R0 = String[],
					E0 = String[],
					I0 = String[],
					PCR_pos_samplingbias = String[],
					PCR_pos_sample_clustering = String[],
					PCR_pos_reporting_clustering = String[],
					PCR_pos_reporting_rate = String[],
					EffPopSize = String[],
					IFR = String[],
                    PosteriorpredicitveP = String[])
    for modelpath in modelpaths
        # modelpath = joinpath(dirname,name)
        D = FileIO.load(modelpath)
        f = first(keys(D))
        model = D[f]
        println("Gathering data for $(model.areaname)")
        name_upr = String([uppercase(s) for s in model.areaname])
        deaths = vec(death_data.deaths[:,death_data.areas .== name_upr][1:enddate_int])
		row = createdatarow(model,deaths,IFR_prior,p_ID,num_sims)
        push!(df,row)
    end
    return df
end


"""
    getincidence_severe_and_deaths(model,hosps,deaths,IFR_priormean,p_IH,p_ID)

Fit a crude IHR and IFR to the data for the area of the fitted model and return posterior mean predictions as a NamedTuple of time series.
"""
function getincidence_severe_and_deaths(model,hosps,deaths,IFR_priormean,p_IH,p_ID)
	inc = incidence_across_samples(model,588)
	IHR_priormean = IFR_priormean*3
	H = sum(hosps[1:143]) #Cut off at 12th July because of poor data
	D = sum(deaths[1:(end-3)])
	IFR_ests = mean.([Erlang(D+1,1/((1/IFR_priormean) + sum(simple_conv(inc.true_incidence[:,n],p_ID)[1:length(deaths[1:(end-3)])]))) for n = 1: 10000])
    IHR_ests = 3*IFR_ests
    if H >= 0 #This catches if we have hosp data for the county
    	IHR_ests = mean.([Erlang(H+1,1/((1/IHR_priormean) + sum(simple_conv(inc.true_incidence[:,n],p_IH)[1:143]))) for n = 1: size(inc.true_incidence,2)])
    end
	inc_predictions = vec(mean(inc.true_incidence,dims = 2))
	death_predictions = mean([death_pred(inc.true_incidence[:,j],IFR_ests[j],p_ID) for j in 1:size(inc.true_incidence,2)])
	hosp_predictions = mean([death_pred(inc.true_incidence[:,j],IHR_ests[j],p_IH) for j in 1:size(inc.true_incidence,2)])
	return (inc_predictions=inc_predictions,hosp_predictions=hosp_predictions,death_predictions=death_predictions)
end

"""
    incidence_severe_death_projectionovercollection(dirname,hosp_data,death_data,IFR_prior,p_IH,p_ID)

Apply getincidence_severe_and_deaths to every fitted model saved in the directory at dirname, and sum their outputs.
"""
function incidence_severe_death_projectionovercollection(dirname,hosp_data,death_data,IFR_prior,p_IH,p_ID)
    modelnames = readdir(dirname)
    totalincidence = zeros(588)
    totalhosp = zeros(588)
    totaldeaths = zeros(588)

    for name in modelnames
        modelpath = joinpath(dirname,name)
        D = FileIO.load(modelpath)
        f = first(keys(D))
        model = D[f]
        println("Gathering data for $(model.areaname)")
        name_upr = String([uppercase(s) for s in model.areaname])
        deaths = death_data.deaths[:,death_data.areas .== name_upr]
        hosps = -1*ones(length(deaths))
        if any(hosp_data.areas .== name_upr)
            hosps = hosp_data.hosps[:,hosp_data.areas .== name_upr]
        end
        data = getincidence_severe_and_deaths(model,hosps,deaths,IFR_prior,p_IH,p_ID)
        totalincidence .+= data.inc_predictions
        totalhosp .+= data.hosp_predictions
        totaldeaths .+= data.death_predictions
    end
    return (totalincidence=totalincidence,
            totalhosp=totalhosp,
            totaldeaths=totaldeaths)
end

"""
    createnoninterventionmodel(model)

Create a deepcopy of the fitted model, and then return a new model where the contact rate stays at a baseline of unity.
"""
function createnoninterventionmodel(model)
    nonintervention_model = deepcopy(model)
    baseline_contactrate = (contactrate = ones(length(model.contactrate_data.contactrate)),
                            date = model.contactrate_data.date,
                            area = model.contactrate_data.area)
    nonintervention_model.contactrate_data = baseline_contactrate
    _prob = make_odeproblemforinference(baseline_contactrate,
                                                            startdate = Date(2020,2,21),
                                                            enddate = Date(2020,8,21))
    nonintervention_model.prob = _prob
    return nonintervention_model
end

"""
    calculatekeyobservables(areamodel::CoVAreaModel,areadeaths,areahosps,IFR_est,p_IH,p_ID)

Calculate the posterior mean and posterior stds of key epidemic observables for model areamodel: peak timing of infection, predicted deaths, predicted hospitalisations, PCR swab testing data and true incidence.
"""
function calculatekeyobservables(areamodel::CoVAreaModel,areadeaths,areahosps,IFR_est,p_IH,p_ID)
   inc = incidence_across_samples(areamodel,588);
   incidence = create_credible_intervals(inc.true_incidence);
   time_length = size(inc.true_incidence,1)

   mean_pred_PCR_NB = zeros(time_length)
   std_pred_PCR_NB = zeros(time_length)
   mean_pred_P_PCR_BB = zeros(time_length)
   std_pred_P_PCR_BB = zeros(time_length)
   mean_pred_incidence = zeros(time_length)
   std_pred_incidence = zeros(time_length)
   mean_true_infected = zeros(time_length)
   std_true_infected = zeros(time_length)

   for t = 1:time_length
        PCRinconthatday_NB = inc.PCR_incidence_NB[t,:]
        PCR_P_onthatday_BB = inc.PCR_P_BB[t,:]
		inconthatday = inc.true_incidence[t,:]
		trueinfectedonday = inc.true_infecteds[t,:]

        mean_pred_PCR_NB[t] = mean(PCRinconthatday_NB)
        std_pred_PCR_NB[t] = std(PCRinconthatday_NB)
		mean_pred_P_PCR_BB[t] = mean(PCR_P_onthatday_BB)
        std_pred_P_PCR_BB[t] = std(PCR_P_onthatday_BB)
        mean_pred_incidence[t] = mean(inconthatday)
        std_pred_incidence[t] = std(inconthatday)
		mean_true_infected[t] = mean(trueinfectedonday)
		std_true_infected[t] = std(trueinfectedonday)
    end

   D = sum(areadeaths[1:(end-3)])
   H = sum(areahosps[1:143]) #Cut off at 12th July because of poor data

   IFR_ests = mean.([Erlang(D+1,1/((1/IFR_est) + sum(simple_conv(inc.true_incidence[:,n],p_ID)[1:length(areadeaths[1:(end-3)])]))) for n = 1:size(inc.true_incidence,2)])
   IHR_ests = 3*IFR_ests
   if H >= 0 #This catches if we have hosp data for the county
       IHR_ests = mean.([Erlang(H+1,1/((1/(3*IFR_est)) + sum(simple_conv(inc.true_incidence[:,n],p_IH)[1:143]))) for n = 1:size(inc.true_incidence,2)])
   end
   death_predictions = mean([death_pred(inc.true_incidence[:,j],IFR_ests[j],p_ID) for j in 1:size(inc.true_incidence,2)])
   hosp_predictions = mean([death_pred(inc.true_incidence[:,j],IHR_ests[j],p_IH) for j in 1:size(inc.true_incidence,2)])

   # deaths_estIFR = create_death_prediction_intervals(death_pred(incidence.mean_pred,IFR_est,p_ID))

   post_mean_peak = mean([findmax(inc.true_incidence[:,k])[2] for k = 1:size(inc.true_incidence,2)])

   mean_R = mean(areamodel.MCMC_results.chain[:,1,1])

   return (area = areamodel.areaname,
            peak = post_mean_peak,
            pred_deaths = death_predictions,
            pred_hosps = hosp_predictions,
            mean_pred_PCR_NB = mean_pred_PCR_NB,
            std_pred_PCR_NB = std_pred_PCR_NB,
			mean_pred_P_PCR_BB = mean_pred_P_PCR_BB,
			std_pred_P_PCR_BB = std_pred_P_PCR_BB,
            mean_pred_incidence = mean_pred_incidence,
            std_pred_incidence = std_pred_incidence,
			mean_true_infected = mean_true_infected,
			std_true_infected = std_true_infected,
            mean_R = mean_R)
end


"""
    calculatekeyobservables(areamodel::CoVAreaModel,areamodel_nointervention::CoVAreaModel,areadeaths,areahosps,IFR_est,p_IH,p_ID)

Calculate the posterior mean and posterior stds of key epidemic observables for model areamodel_nointervention: peak timing of infection, predicted deaths, predicted hospitalisations, PCR swab testing data and true incidence.
NB: IHR and IFR are fitted using the best fit model to the epidemic in that area, areamodel, then used to create counter-factual projections of what might have happened if contact rates had stayed at baseline.
"""
function calculatekeyobservables(areamodel::CoVAreaModel,areamodel_nointervention::CoVAreaModel,areadeaths,areahosps,IFR_est,p_IH,p_ID)
   inc = incidence_across_samples(areamodel,588);
   inc_nointervention = incidence_across_samples(areamodel_nointervention,588);
   time_length = size(inc.true_incidence,1)

   mean_pred_PCR_NB = zeros(time_length)
   std_pred_PCR_NB = zeros(time_length)
   mean_pred_P_PCR_BB = zeros(time_length)
   std_pred_P_PCR_BB = zeros(time_length)
   mean_pred_incidence = zeros(time_length)
   std_pred_incidence = zeros(time_length)
   mean_true_infected = zeros(time_length)
   std_true_infected = zeros(time_length)

   for t = 1:time_length
        PCRinconthatday_NB = inc.PCR_incidence_NB[t,:]
        PCR_P_onthatday_BB = inc.PCR_P_BB[t,:]
		inconthatday = inc.true_incidence[t,:]
		trueinfectedonday = inc.true_infecteds[t,:]

        mean_pred_PCR_NB[t] = mean(PCRinconthatday_NB)
        std_pred_PCR_NB[t] = std(PCRinconthatday_NB)
		mean_pred_P_PCR_BB[t] = mean(PCR_P_onthatday_BB)
        std_pred_P_PCR_BB[t] = std(PCR_P_onthatday_BB)
        mean_pred_incidence[t] = mean(inconthatday)
        std_pred_incidence[t] = std(inconthatday)
		mean_true_infected[t] = mean(trueinfectedonday)
		std_true_infected[t] = std(trueinfectedonday)
    end

   D = sum(areadeaths[1:(end-3)])
   H = sum(areahosps[1:143]) #Cut off at 12th July because of poor data

   IFR_ests = mean.([Erlang(D+1,1/((1/IFR_est) + sum(simple_conv(inc.true_incidence[:,n],p_ID)[1:length(areadeaths[1:(end-3)])]))) for n = 1:size(inc.true_incidence,2)])
   IHR_ests = 3*IFR_ests
   if H >= 0 #This catches if we have hosp data for the county
       IHR_ests = mean.([Erlang(H+1,1/((1/(3*IFR_est)) + sum(simple_conv(inc.true_incidence[:,n],p_IH)[1:143]))) for n = 1:size(inc.true_incidence,2)])
   end
   death_predictions = mean([death_pred(inc.true_incidence[:,j],IFR_ests[j],p_ID) for j in 1:size(inc.true_incidence,2)])
   hosp_predictions = mean([death_pred(inc.true_incidence[:,j],IHR_ests[j],p_IH) for j in 1:size(inc.true_incidence,2)])

   # deaths_estIFR = create_death_prediction_intervals(death_pred(incidence.mean_pred,IFR_est,p_ID))

   post_mean_peak = mean([findmax(inc.true_incidence[:,k])[2] for k = 1:size(inc.true_incidence,2)])

   mean_R = mean(areamodel.MCMC_results.chain[:,1,1])

   return (area = areamodel.areaname,
            peak = post_mean_peak,
            pred_deaths = death_predictions,
            pred_hosps = hosp_predictions,
            mean_pred_PCR_NB = mean_pred_PCR_NB,
            std_pred_PCR_NB = std_pred_PCR_NB,
			mean_pred_P_PCR_BB = mean_pred_P_PCR_BB,
			std_pred_P_PCR_BB = std_pred_P_PCR_BB,
            mean_pred_incidence = mean_pred_incidence,
            std_pred_incidence = std_pred_incidence,
			mean_true_infected = mean_true_infected,
			std_true_infected = std_true_infected,
            mean_R = mean_R)
end


"""
    gatherdatafrommodelsfornonintervention(dirname,death_data,hosp_data,IFR_est,p_IH,p_ID)

Calculate key observables for a counter-factual model (with baseline contact rates set at one) for each fitted model saved in the directory at dirname.
"""
function gatherdatafrommodelsfornonintervention(dirname,death_data,hosp_data,IFR_est,p_IH,p_ID)
    modelnames = readdir(dirname)
    collecteddata = Any[]
    for name in modelnames
        modelpath = joinpath(dirname,name)
        D = FileIO.load(modelpath)
        f = first(keys(D))
        model = D[f]
        noninterventionmodel = createnoninterventionmodel(model)
        println("Gathering data for $(noninterventionmodel.areaname)")
        name_upr = String([uppercase(s) for s in noninterventionmodel.areaname])
        deaths = death_data.deaths[:,death_data.areas .== name_upr]
        hosps = -1*ones(length(deaths))
        if any(hosp_data.areas .== name_upr)
            hosps = hosp_data.hosps[:,hosp_data.areas .== name_upr]
        end
        datafrommodel = calculatekeyobservables(model,noninterventionmodel,deaths,hosps,IFR_est,p_IH,p_ID)
        push!(collecteddata,datafrommodel)
    end
    return collecteddata
end

"""
    gatherdatafrommodels(dirname,death_data,hosp_data,IFR_est,p_IH,p_ID)

Calculate key observables for each fitted model saved in the directory at dirname.
"""
function gatherdatafrommodels(dirname,death_data,hosp_data,IFR_est,p_IH,p_ID)
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
        hosps = -1*ones(length(deaths))
        if any(hosp_data.areas .== name_upr)
            hosps = hosp_data.hosps[:,hosp_data.areas .== name_upr]
        end
        datafrommodel = calculatekeyobservables(model,deaths,hosps,IFR_est,p_IH,p_ID)
        push!(collecteddata,datafrommodel)
    end
    return collecteddata
end
