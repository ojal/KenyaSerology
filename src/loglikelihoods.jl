
"""
loglikelihood_contactratemodel(θ,model::CoVAreaModel)

This version of the log-likelihood assumes that negative PCR tests on each day are unavailable, and therefore,
uses a NegativeBinomial model (with detection probability and clustering coefficient) to link simulation to observation.
First 21 days didn't have testing, and don't contribute to log-likelihood.
R(t) is assumed to be proportional to the transformed Google mobility data.
"""
function loglikelihood_contactratemodel(θ,model::CoVAreaModel)
    @unpack E₀,I₀,R,α,p_test = θ #Parameters to fit 1) initial conditions and R₀ 2) parameters for Neg. bin., 3) parameters for beta binomial model
    @unpack PCR_cases,sero_cases,N,σ,γ,prob,PCR_array,sero_array,log_priors,PCR_sensitivity,PCR_specificity,sero_specificity,sero_sensitivity = model
    T = eltype(R)
    u0 = convert.(T, [N-E₀-I₀,E₀,I₀,0.,0.])
    p = convert.(T,[R,σ,γ,N])
    try
        sol = solve(prob, BS3(); u0=u0, p=p, saveat = 1)
        incidence = get_incidence(sol)
        num_PCR = simple_conv(incidence,PCR_array)
        p_sero = sero_sensitivity*simple_conv(incidence,sero_array)/N
        LL = T(0.)
        for t in 59:min(length(incidence),size(PCR_cases,1),size(sero_cases,1))
            d_PCR = PCR_cases[t]#Detected number of PCR pos that day
            d_sero_pos = sero_cases[t,1] #Detected number of sero pos that day
            n_sero = sero_cases[t,1] + sero_cases[t,2] #Number sero tested on that day
            pred_num_PCR = num_PCR[t] #Predicted number of PCR pos
            pred_p_sero = p_sero[t] #Prediction proportion sero pos
            #Convert from μ,α parameterisation to p,r parameterisation
            μ = p_test*pred_num_PCR +0.001
            σ² = μ + α*μ^2
            p_negbin = 1 - (α*μ^2/σ²)
            r_negbin = 1/α
            p_hat =  pred_p_sero + (1-pred_p_sero)*(1-sero_specificity)

            LL += logpdf(NegativeBinomial(r_negbin,p_negbin),d_PCR)#likelihood contribution from PCR testing
            if n_sero > 0
                LL += logpdf(Binomial(n_sero,p_hat ),d_sero_pos)#Likelihood contribution from sero testing
                    # println(logpdf(Binomial(n_sero,pred_p_sero),d_sero_pos))
            end
            # println(LL)
        end
        return LL + log_priors(θ) #+ logpdf(Binomial(sum(sero_cases),p_sero[end]),sum(sero_cases[:,1]))
    catch
        return T(-Inf)
    end
end


"""
loglikelihood_contactratemodelBB(θ,model::CoVAreaModel)

This version of the log-likelihood assumes that negative PCR tests on each day are unavailable, and therefore,
uses a NegativeBinomial model (with detection probability and clustering coefficient) to link simulation to observation.
First 21 days didn't have testing, and don't contribute to log-likelihood.

Uncertainty in the serology testing means that we assume that it is drawn from a Betabinomial distribution, with an inferred shrinkage factor M_BB

R(t) is assumed to be proportional to the transformed Google mobility data.
"""
function loglikelihood_contactratemodelBB(θ,model::CoVAreaModel)
    @unpack E₀,I₀,R,α,p_test = θ #Parameters to fit 1) initial conditions and R₀ 2) parameters for Neg. bin., 3) parameters for beta binomial model
    @unpack PCR_cases,sero_cases,N,σ,γ,prob,PCR_array,sero_array,log_priors,PCR_sensitivity,PCR_specificity,sero_specificity,sero_sensitivity,M_BB = model
    T = eltype(R)
    u0 = convert.(T, [N-E₀-I₀,E₀,I₀,0.,0.])
    p = convert.(T,[R,σ,γ,N])
    try
        sol = solve(prob, BS3(); u0=u0, p=p, saveat = 1)
        incidence = get_incidence(sol)
        num_PCR = simple_conv(incidence,PCR_array)
        p_sero = sero_sensitivity*simple_conv(incidence,sero_array)/N
        LL = T(0.)
        for t in 59:min(length(incidence),size(PCR_cases,1))
            d_PCR = PCR_cases[t]#Detected number of PCR pos that day
            d_sero_pos = sero_cases[t,1] #Detected number of sero pos that day
            n_sero = sero_cases[t,1] + sero_cases[t,2] #Number sero tested on that day
            pred_num_PCR = num_PCR[t] #Predicted number of PCR pos
            pred_p_sero = p_sero[t] #Prediction proportion sero pos
            #Convert from μ,α parameterisation to p,r parameterisation
            μ = p_test*pred_num_PCR + 0.001
            σ² = μ + α*μ^2
            p_negbin = 1 - (α*μ^2/σ²)
            r_negbin = 1/α
            #Convert from the μ_BB,M_BB parameterisation of the Betabinomial distribution to the α_BB,β_BB parameterisation
            p_hat =  pred_p_sero + (1-pred_p_sero)*(1-sero_specificity)
            α_BB = M_BB*p_hat
            β_BB = M_BB*(1-p_hat)
            LL += logpdf(NegativeBinomial(r_negbin,p_negbin),d_PCR)#likelihood contribution from PCR testing
            if n_sero > 0
                LL += logpdf(BetaBinomial(n_sero,α_BB,β_BB),d_sero_pos)#Likelihood contribution from sero testing
                    # println(logpdf(Binomial(n_sero,pred_p_sero),d_sero_pos))
            end
            # println(LL)
        end
        return LL + log_priors(θ) #+ logpdf(Binomial(sum(sero_cases),p_sero[end]),sum(sero_cases[:,1]))
    catch
        return T(-Inf)
    end
end

"""
loglikelihood_contactratemodelBB_OR(θ,model::CoVAreaModel)

This version of the log-likelihood assumes that negative PCR tests on each day are unavailable, and therefore,
uses a NegativeBinomial model (with detection probability and clustering coefficient) to link simulation to observation.
First 21 days after 21st Feb didn't have testing, and don't contribute to log-likelihood.

Uncertainty in the serology testing means that we assume that it is drawn from a Betabinomial distribution, with an inferred shrinkage factor M_BB,

There is an odds-ratio parameter (χ) for the relative likelihood of a true seropositive being in the blood donation sample

R(t) is assumed to be proportional to the transformed Google mobility data.
"""
function loglikelihood_contactratemodelBB_OR(θ,model::CoVAreaModel)
    @unpack E₀,I₀,R,α,χ,p_test = θ #Parameters to fit 1) initial conditions and R₀ 2) parameters for Neg. bin., 3) parameters for beta binomial model
    @unpack PCR_cases,sero_cases,N,σ,γ,prob,PCR_array,sero_array,log_priors,PCR_sensitivity,PCR_specificity,sero_specificity,sero_sensitivity,M_BB = model
    T = eltype(R)
    u0 = convert.(T, [N-E₀-I₀,E₀,I₀,0.,0.])
    p = convert.(T,[R,σ,γ,N])
    try
        sol = solve(prob, BS3(); u0=u0, p=p, saveat = 1)
        incidence = get_incidence(sol)
        num_PCR = simple_conv(incidence,PCR_array)
        N_sero = simple_conv(incidence,sero_array)
        LL = T(0.)
        for t in 59:min(length(incidence),size(PCR_cases,1))
            d_PCR = PCR_cases[t]#Detected number of PCR pos that day
            d_sero_pos = sero_cases[t,1] #Detected number of sero pos that day
            n_sero = sero_cases[t,1] + sero_cases[t,2] #Number sero tested on that day
            pred_num_PCR = num_PCR[t] #Predicted number of PCR pos
            pred_N_sero = N_sero[t] #Prediction number sero pos
            #Convert from μ,α parameterisation to p,r parameterisation
            μ = p_test*pred_num_PCR + 0.001
            σ² = μ + α*μ^2
            p_negbin = 1 - (α*μ^2/σ²)
            r_negbin = 1/α
            #Convert from the μ_BB,M_BB parameterisation of the Betabinomial distribution to the α_BB,β_BB parameterisation
            #account for OR of sero-pos
            p_hat = PCR_sensitivity*(χ*pred_N_sero)/(N + (χ-1)*pred_N_sero)
            p_hat += (1-p_hat)*(1-sero_specificity)
            α_BB = M_BB*p_hat
            β_BB = M_BB*(1-p_hat)
            LL += logpdf(NegativeBinomial(r_negbin,p_negbin),d_PCR)#likelihood contribution from PCR testing
            if n_sero > 0
                LL += logpdf(BetaBinomial(n_sero,α_BB,β_BB),d_sero_pos)#Likelihood contribution from sero testing
                    # println(logpdf(Binomial(n_sero,pred_p_sero),d_sero_pos))
            end
            # println(LL)
        end
        return LL + log_priors(θ) #+ logpdf(Binomial(sum(sero_cases),p_sero[end]),sum(sero_cases[:,1]))
    catch
        return T(-Inf)
    end
end

"""
    loglikelihood_contactratemodelBB_Peff(θ,model::CoVAreaModel)

This version of the log-likelihood assumes that negative PCR tests on each day are unavailable, and therefore,
uses a NegativeBinomial model (with detection probability and clustering coefficient) to link simulation to observation.
First 59 days after 21st Feb didn't have testing or had noticeably lower testing rates, and don't contribute to log-likelihood.

Uncertainty in the serology testing means that we assume that it is drawn from a Betabinomial distribution, with an inferred shrinkage factor M_BB,

We assume that the extra heterogeneity in the model is captured by having an effective population size P_eff

R(t) is assumed to be proportional to the transformed Google mobility data.
"""
function loglikelihood_contactratemodelBB_Peff(θ,model::CoVAreaModel)
    @unpack E₀,I₀,R,α,p_test,P_eff = θ #Parameters to fit 1) initial conditions and R₀ 2) parameters for Neg. bin., 3) parameters for beta binomial model
    @unpack PCR_cases,sero_cases,N,σ,γ,prob,PCR_array,sero_array,log_priors,PCR_sensitivity,PCR_specificity,sero_specificity,sero_sensitivity,M_BB,relative_testing_rate = model
    T = eltype(R)
    N_eff = N*P_eff
    u0 = convert.(T, [N_eff-E₀-I₀,E₀,I₀,0.,0.])
    p = convert.(T,[R,σ,γ,N_eff])
    try
        sol = solve(prob, BS3(); u0=u0, p=p, saveat = 1)
        incidence = get_incidence(sol)
        num_PCR = simple_conv(incidence,PCR_array)
        p_sero = sero_sensitivity*simple_conv(incidence,sero_array)/N
        LL = T(0.)
        for t in 59:min(length(incidence),size(PCR_cases,1))
            d_PCR = PCR_cases[t]#Detected number of PCR pos that day
            d_sero_pos = 0
            n_sero = 0
            if t<= size(sero_cases,1)
                d_sero_pos = sero_cases[t,1] #Detected number of sero pos that day
                n_sero = sero_cases[t,1] + sero_cases[t,2] #Number sero tested on that day
            end
            pred_num_PCR = num_PCR[t] #Predicted number of PCR pos
            pred_p_sero = p_sero[t] #Prediction proportion sero pos
            #Convert from μ,α parameterisation to p,r parameterisation
            μ = relative_testing_rate[t]*p_test*pred_num_PCR + 0.001
            σ² = μ + α*μ^2
            p_negbin = 1 - (α*μ^2/σ²)
            r_negbin = 1/α
            #Convert from the μ_BB,M_BB parameterisation of the Betabinomial distribution to the α_BB,β_BB parameterisation
            p_hat =  pred_p_sero + (1-pred_p_sero)*(1-sero_specificity)
            α_BB = M_BB*p_hat
            β_BB = M_BB*(1-p_hat)
            LL += logpdf(NegativeBinomial(r_negbin,p_negbin),d_PCR)#likelihood contribution from PCR testing
            if n_sero > 0
                LL += logpdf(BetaBinomial(n_sero,α_BB,β_BB),d_sero_pos)#Likelihood contribution from sero testing
                    # println(logpdf(Binomial(n_sero,pred_p_sero),d_sero_pos))
            end
            # println(LL)
        end
        return LL + log_priors(θ) #+ logpdf(Binomial(sum(sero_cases),p_sero[end]),sum(sero_cases[:,1]))
    catch
        return T(-Inf)
    end
end




"""
    loglikelihood_without_negPCR(θ,model::CoVAreaModel)

This version of the log-likelihood assumes that negative PCR tests on each day are unavailable, and therefore,
uses a NegativeBinomial model (with detection probability and clustering coefficient) to link simulation to observation
"""
function loglikelihood_without_negPCR(θ,model::CoVAreaModel)
    @unpack E₀,I₀,R_eff1,R_eff2,α,p_test = θ #Parameters to fit 1) initial conditions and R₀ 2) parameters for Neg. bin.
    @unpack PCR_cases,sero_cases,N,σ,γ,prob,PCR_array,sero_array,log_priors = model
    T = eltype(R_eff1)
    u0 = convert.(T, [N-E₀-I₀,E₀,I₀,0.,0.])
    p = convert.(T,[R_eff1,R_eff2,σ,γ,N])
    try
        sol = solve(prob, BS3(); u0=u0, p=p, saveat = 1)
        incidence = get_incidence(sol)
        num_PCR = simple_conv(incidence,PCR_array)
        p_sero = simple_conv(incidence,sero_array)/N
        log_likelihood = T(0.)
        for t in 1:min(length(incidence),size(PCR_cases,1))
            d_PCR = PCR_cases[t]#Detected number of PCR pos that day
            d_sero_pos = sero_cases[t,1] #Detected number of sero pos that day
            n_sero = sero_cases[t,1] + sero_cases[t,2] #Number sero tested on that day
            pred_num_PCR = num_PCR[t] #Predicted number of PCR pos
            pred_p_sero = p_sero[t] #Prediction proportion sero pos
            #Convert from μ,α parameterisation to p,r parameterisation
            μ = p_test*pred_num_PCR + 0.001
            σ² = μ + α*μ^2
            p_negbin = 1 - (α*μ^2/σ²)
            r_negbin = 1/α
            log_likelihood += logpdf(NegativeBinomial(r_negbin,p_negbin),d_PCR)#likelihood contribution from PCR testing
            if n_sero > 0
                log_likelihood += logpdf(Binomial(n_sero,pred_p_sero),d_sero_pos)#Likelihood contribution from sero testing
                    # println(logpdf(Binomial(n_sero,pred_p_sero),d_sero_pos))
            end
        end
        return log_likelihood + log_priors(θ) #+ logpdf(Binomial(sum(sero_cases),p_sero[end]),sum(sero_cases[:,1]))
    catch
        return T(-Inf)
    end
end

"""
    loglikelihood_with_definitenegPCR(θ,model::CoVAreaModel)

Calculate the log-likelihood, or posterior log-density if priors are non-flat, for parameters `θ`. A key assumption is that negative test rates are always available.
This version assumes that the PCR tests each day are a biased sample of the PCR⁺ individuals in the population, and moreover, there is intra-cluster correlation which inflates the day-to-day variance of the PCR testing data.

Parameters:

* `R`: Baseline reproductive ratio.
* `E₀`: Number of exposed infecteds on 21st Feburary 2020.
* `I₀`: Number of infectious infecteds on 21st Feburary 2020.
* `χ`: Odds-ratio of a PCR⁺ invidual being accurately tested relative to a PCR⁻ individual.
* `M_PCR`: The 'effective sample size' parameter for the PCR testing. The intra-cluster correlation is ρ = 1/(M_PCR + 1).
* `P_eff`: Effective population size parameter.
"""
function loglikelihood_with_definitenegPCR(θ,model::KenyaSerology.CoVAreaModel)
    @unpack R,E₀,I₀,χ,M_PCR,P_eff = θ #Parameters to fit 1) initial conditions and R₀ 2) parameters for Neg. bin.
    @unpack PCR_cases,sero_cases,N,σ,γ,prob,PCR_array,sero_array,log_priors,PCR_sensitivity,PCR_specificity,sero_specificity,sero_sensitivity,M_BB,relative_testing_rate = model
    T = eltype(R)
    N_eff = N*P_eff
    u0 = convert.(T, [N_eff,E₀,I₀,0.,0.])
    p = convert.(T,[R,σ,γ,N_eff])
    try
        sol = solve(prob, BS3(); u0=u0, p=p, saveat = 1)
        incidence = KenyaSerology.get_incidence(sol)
        num_PCR = KenyaSerology.simple_conv(incidence,PCR_array)
        p_sero = sero_sensitivity*KenyaSerology.simple_conv(incidence,sero_array)/N
        LL = T(0.)
        for t in 59:min(length(incidence),size(PCR_cases,1))
            d_PCR_pos = PCR_cases[t,1]#Detected number of PCR pos that day
            n_PCR = PCR_cases[t,1] + PCR_cases[t,2]#Number PCR tested on that day
            d_sero_pos = 0
            n_sero = 0
            if t<= size(sero_cases,1)
                d_sero_pos = sero_cases[t,1] #Detected number of sero pos that day
                n_sero = sero_cases[t,1] + sero_cases[t,2] #Number sero tested on that day
            end
            #Convert from p_PCR,inv_M_PCR parameterisation of Beta-Binomial to standard α,β parameterisation
            if n_PCR > 0
                p_PCR = χ*num_PCR[t]/((χ-1)*num_PCR[t] + N)
                LL += logpdf(BetaBinomial(n_PCR,p_PCR*M_PCR,(1-p_PCR)*M_PCR),d_PCR_pos)#likelihood contribution from PCR testing
            end
            #Convert from the p_hat,M_BB parameterisation of the Betabinomial distribution to the α_BB,β_BB parameterisation
            if n_sero > 0
                p_hat =   p_sero[t] + (1- p_sero[t])*(1-sero_specificity)
                α_BB = M_BB*p_hat
                β_BB = M_BB*(1-p_hat)
                LL += logpdf(BetaBinomial(n_sero,α_BB,β_BB),d_sero_pos)#Likelihood contribution from sero testing
            end
        end
        return LL + log_priors(θ)
    catch
        return T(-Inf)
    end
end


"""
    loglikelihood_with_negPCR(θ,model::CoVAreaModel)

Calculate the log-likelihood, or posterior log-density if priors are non-flat, for parameters `θ`. Negative test rates are used when available, otherwise there is a default to using a Neg. Binomial model.
This version assumes that the PCR tests each day are a biased sample of the PCR⁺ individuals in the population, and moreover, there is intra-cluster correlation which inflates the day-to-day variance of the PCR testing data.

Parameters:

* `R`: Baseline reproductive ratio.
* `E₀`: Number of exposed infecteds on 21st Feburary 2020.
* `I₀`: Number of infectious infecteds on 21st Feburary 2020.
* `χ`: Odds-ratio of a PCR⁺ invidual being accurately tested relative to a PCR⁻ individual.
* `M_PCR`: The 'effective sample size' parameter for the PCR testing. The intra-cluster correlation is ρ = 1/(M_PCR + 1).
* `p_test`: The effective sampling rate of PCR⁺ individuals when negative tests are not available.
* `α`: The clustering factor of the Neg. binomial model used when negative tests are not available.
* `P_eff`: Effective population size parameter.
"""
function loglikelihood_with_negPCR(θ,model::KenyaSerology.CoVAreaModel)
    @unpack R,E₀,I₀,χ,M_PCR,α,p_test,P_eff = θ #Parameters to fit 1) initial conditions and R₀ 2) parameters for Neg. bin.
    @unpack PCR_cases,sero_cases,N,σ,γ,prob,PCR_array,sero_array,log_priors,PCR_sensitivity,PCR_specificity,sero_specificity,sero_sensitivity,M_BB,relative_testing_rate = model
    T = eltype(R)
    N_eff = N*P_eff
    u0 = convert.(T, [N_eff,E₀,I₀,0.,0.])
    p = convert.(T,[R,σ,γ,N_eff])
    try
        sol = solve(prob, BS3(); u0=u0, p=p, saveat = 1)
        incidence = KenyaSerology.get_incidence(sol)
        num_PCR = KenyaSerology.simple_conv(incidence,PCR_array)
        p_sero = sero_sensitivity*KenyaSerology.simple_conv(incidence,sero_array)/N
        LL = T(0.)
        for t in 59:min(length(incidence),size(PCR_cases,1))
            d_PCR_pos = PCR_cases[t,1]#Detected number of PCR pos that day
            d_PCR_neg = PCR_cases[t,2]#Detected number of PCR neg that day, -1 === no negative tests available
            d_sero_pos = 0
            n_sero = 0
            if t<= size(sero_cases,1)
                d_sero_pos = sero_cases[t,1] #Detected number of sero pos that day
                n_sero = sero_cases[t,1] + sero_cases[t,2] #Number sero tested on that day
            end
            if d_PCR_neg >= 0
                #Convert from p_PCR,inv_M_PCR parameterisation of Beta-Binomial to standard α,β parameterisation
                p_PCR = χ*num_PCR[t]/((χ-1)*num_PCR[t] + N)
                LL += logpdf(BetaBinomial(d_PCR_pos+d_PCR_neg,p_PCR*M_PCR,(1-p_PCR)*M_PCR),d_PCR_pos)#likelihood contribution from PCR testing --- given
            end
            if d_PCR_neg == -1
                #Convert from μ,α parameterisation to p,r parameterisation
                μ = relative_testing_rate[t]*p_test*num_PCR[t] + 0.001
                σ² = μ + α*μ^2
                p_negbin = 1 - (α*μ^2/σ²)
                r_negbin = 1/α
                LL += logpdf(NegativeBinomial(r_negbin,p_negbin),d_PCR_pos)#likelihood contribution from PCR testing
            end
            #Convert from the p_hat,M_BB parameterisation of the Betabinomial distribution to the α_BB,β_BB parameterisation
            if n_sero > 0
                p_hat =   p_sero[t] + (1- p_sero[t])*(1-sero_specificity)
                α_BB = M_BB*p_hat
                β_BB = M_BB*(1-p_hat)
                LL += logpdf(BetaBinomial(n_sero,α_BB,β_BB),d_sero_pos)#Likelihood contribution from sero testing
            end
        end
        return LL + log_priors(θ)
    catch
        return T(-Inf)
    end
end
