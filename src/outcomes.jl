"""
    predict_incidence_and_prev(;pred_incidence::Vector{Float64},
                                        pred_hospitalisations::Vector{Float64},
                                        pred_deaths::Vector{Float64},
										symptomaticrate::Float64,
										prop_critical::Float64,
										p_IS::Vector{Float64},
										Q_R::Vector{Float64},
                                        Q_ICUH::Vector{Float64},
                                        Q_HR::Vector{Float64},
                                        Q_ICUR::Vector{Float64},
										F_R::Vector{Float64})

Augment a prediction of incidence of infection, hospitalisation and death in an area with predictions of symptomatic/asymptomatic incidence, critical care (ICU) incidence                                       

                                        # Arguments
- `pred_incidence::Vector{Float64}`: Daily predictions of new infections
- `pred_hospitalisations::Vector{Float64}`: Daily predictions of new hospitalisations.
- `pred_deaths::Vector{Float64}`: Daily predictions of new deaths.
- `symptomaticrate::Float64`: Broad population symptomatic rate.
- `prop_critical::Float64`: Proportion of hospitalised cases that require ICU treatment
- `p_IS::Vector{Float64}`: Distribution of days between infection and symptoms (for symptomatic cases).
- `Q_R::Vector{Float64}`: Upper tail function for days between infection and recovery
- `Q_ICUH::Vector{Float64}`: Upper tail function for days between arriving in ICU, surviving, and returning to less specialist wards
- `Q_HR::Vector{Float64}`: Upper tail function for days between arriving at hospital, remaining severe but not critical and then getting discharged
- `Q_ICUR::Vector{Float64}`: Upper tail function for days between arriving in ICU, surviving, returning to less specialist wards and then discharge
- `F_R::Vector{Float64}`: Distribution function for days until recovery from infection 
"""
function predict_incidence_and_prev(;pred_incidence::Vector{Float64},
                                        pred_hospitalisations::Vector{Float64},
                                        pred_deaths::Vector{Float64},
										symptomaticrate::Float64,
										prop_critical::Float64,
										p_IS::Vector{Float64},
										Q_R::Vector{Float64},
                                        Q_ICUH::Vector{Float64},
                                        Q_HR::Vector{Float64},
                                        Q_ICUR::Vector{Float64},
										F_R::Vector{Float64})
	asymp_incidence = (1 - symptomaticrate) * simple_conv(pred_incidence, p_IS)
	symp_incidence = symptomaticrate * simple_conv(pred_incidence, p_IS)
	severe_incidence = (1 - prop_critical) * pred_hospitalisations
	crit_incidence = prop_critical * pred_hospitalisations
	death_incidence = pred_deaths
	asymp_prev = simple_conv(asymp_incidence, Q_R)
	critical_prev = simple_conv(crit_incidence, Q_ICUH)
	severe_prev = simple_conv(severe_incidence, Q_HR) .+ simple_conv(crit_incidence, Q_ICUR) .- critical_prev
	mild_prev = simple_conv(symp_incidence, Q_R) .- severe_prev .- critical_prev
	recovered_prev = simple_conv(asymp_incidence .+ symp_incidence, F_R)										
	return (asymp_incidence = asymp_incidence,
            symp_incidence = symp_incidence,
            severe_incidence = severe_incidence,
            crit_incidence = crit_incidence,
            death_incidence = death_incidence,
            asymp_prev = asymp_prev,
            critical_prev = critical_prev,
            severe_prev = severe_prev,
            mild_prev = mild_prev,
            recovered_prev = recovered_prev)
end

"""
    predict_incidence_and_prev(modelfit;
                                    pred_incidence::Vector{Float64},
                                    pred_hospitalisations::Vector{Float64},
                                    pred_deaths::Vector{Float64},
                                    symptomaticrate::Float64,
                                    prop_critical::Float64,
                                    p_IS::Vector{Float64},
                                    Q_R::Vector{Float64},
                                    Q_ICUH::Vector{Float64},
                                    Q_HR::Vector{Float64},
                                    Q_ICUR::Vector{Float64},
                                    F_R::Vector{Float64})

Convert a collected modelfit element into a hospitalisation prediction
"""
function predict_incidence_and_prev(modelfit;
                                    symptomaticrate::Float64,
                                    prop_critical::Float64,
                                    p_IS::Vector{Float64},
                                    Q_R::Vector{Float64},
                                    Q_ICUH::Vector{Float64},
                                    Q_HR::Vector{Float64},
                                    Q_ICUR::Vector{Float64},
                                    F_R::Vector{Float64})
	incidence = modelfit.mean_pred_incidence
	hospitalisations = modelfit.pred_hosps_ftc
	deaths = modelfit.pred_deaths_ftc
    name = modelfit.area
    prediction = predict_incidence_and_prev(pred_incidence=incidence,
                                        pred_hospitalisations=hospitalisations,
                                        pred_deaths=deaths,
                                        symptomaticrate=symptomaticrate,
                                        prop_critical=prop_critical,
                                        p_IS=p_IS,
                                        Q_R=Q_R,
                                        Q_ICUH=Q_ICUH,
                                        Q_HR=Q_HR,
                                        Q_ICUR=Q_ICUR,
                                        F_R=F_R)
	return prediction, name
end

"""
    predict_incidence_and_prev(incidence_array,p_IH,p_ID,days;
                                    symptomaticrate::Float64,
                                    prop_critical::Float64,
                                    p_IS::Vector{Float64},
                                    Q_R::Vector{Float64},
                                    Q_ICUH::Vector{Float64},
                                    Q_HR::Vector{Float64},
                                    Q_ICUR::Vector{Float64},
                                    F_R::Vector{Float64})

Apply outcome predictions to each of a set of posterior draws for the incidence curves.                                    
"""
function predict_incidence_and_prev(incidence_array,p_IH,p_ID,days;
                                    symptomaticrate::Float64,
                                    prop_critical::Float64,
                                    p_IS::Vector{Float64},
                                    Q_R::Vector{Float64},
                                    Q_ICUH::Vector{Float64},
                                    Q_HR::Vector{Float64},
                                    Q_ICUR::Vector{Float64},
                                    F_R::Vector{Float64})
        asymp_incidence = zeros(days,size(incidence_array,2))
        symp_incidence = similar(asymp_incidence)
        severe_incidence =similar(asymp_incidence)
        crit_incidence = similar(asymp_incidence)
        death_incidence = similar(asymp_incidence)
        asymp_prev = similar(asymp_incidence)
        critical_prev = similar(asymp_incidence)
        severe_prev = similar(asymp_incidence)
        mild_prev = similar(asymp_incidence)
        recovered_prev = similar(asymp_incidence)
        for k = 1:size(incidence_array,2)
                incidence = incidence_array[:,k]
                hospitalisations = KenyaSerology.simple_conv(incidence,0.001*p_IH)
                deaths = KenyaSerology.simple_conv(incidence,0.001*p_ID)
                prediction = KenyaSerology.predict_incidence_and_prev(pred_incidence=incidence,
                                        pred_hospitalisations=hospitalisations,
                                        pred_deaths=deaths,
                                        symptomaticrate=symptomaticrate,
                                        prop_critical=prop_critical,
                                        p_IS=p_IS,
                                        Q_R=Q_R,
                                        Q_ICUH=Q_ICUH,
                                        Q_HR=Q_HR,
                                        Q_ICUR=Q_ICUR,
                                        F_R=F_R)
                asymp_incidence[:,k] = prediction.asymp_incidence
                symp_incidence[:,k]  = prediction.symp_incidence
                severe_incidence[:,k]  = prediction.severe_incidence
                crit_incidence[:,k]  = prediction.crit_incidence
                death_incidence[:,k]  = prediction.death_incidence
                asymp_prev[:,k]  = prediction.asymp_prev
                critical_prev[:,k]  = prediction.critical_prev
                severe_prev[:,k]  = prediction.severe_prev
                mild_prev[:,k]  = prediction.mild_prev
                recovered_prev[:,k]  = prediction.recovered_prev                
        end
        return (asymp_incidence = asymp_incidence,
            symp_incidence = symp_incidence,
            severe_incidence = severe_incidence,
            crit_incidence = crit_incidence,
            death_incidence = death_incidence,
            asymp_prev = asymp_prev,
            critical_prev = critical_prev,
            severe_prev = severe_prev,
            mild_prev = mild_prev,
            recovered_prev = recovered_prev)
end