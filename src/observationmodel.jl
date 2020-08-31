
"""
Callable on CoVAreaModel object. This uses the log_likelihood of the object.
"""
function (model::CoVAreaModel)(θ)
        return model.log_likelihood(θ,model)
end



function inferparameters!(areamodel::CoVAreaModel,samples,trans,q₀)
        println("Starting parameter inference")
        n = length(q₀)
        ℓ = TransformedLogDensity(trans, areamodel)#transformed log-likelihood
        ∇ℓ = LogDensityProblems.ADgradient(:ForwardDiff, ℓ)#transformed log-likelihood gradient wrt the parameters
        D = Diagonal([0.1,1.,1.,0.1,0.1,0.1].^2)
        results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇ℓ, samples,
                                initialization = (q = q₀,κ = GaussianKineticEnergy(D) ))
                                # warmup_stages = default_warmup_stages(; M =  DynamicHMC.Symmetric))
        transformed_results = TransformVariables.transform.(trans,results.chain)
        val = zeros(length(transformed_results),length(transformed_results[1]),1)
        for i = 1:size(val,1),j = 1:size(val,2)
                val[i,j,1] = transformed_results[i][j]
        end
        chn = Chains(val,[String(k) for k in keys(transformed_results[1])])
        areamodel.MCMC_results = MCMCResults(chn,
                                                [areamodel(transformed_results[i]) - areamodel.log_priors(transformed_results[i]) for i = 1:length(transformed_results)],
                                                results.tree_statistics)
        return nothing
end

function inferparameters!(areamodel::CoVAreaModel,samples,trans)
        println("Searching for a good initial condition")
        f(x) = -areamodel((R=x[1],E₀=x[2],I₀ = x[3],α=x[4],p_test=x[5],P_eff=x[6]))
        searchrange = [(1.,4.),(0.,10.),(0.,10.),(0.,3.),(0.00000,0.001),(0.25,1.)]
        res = bboptimize(f; SearchRange = searchrange,PopulationSize=100)
        q₀ = best_candidate(res)
        inferparameters!(areamodel,samples,trans,q₀)
        return nothing
end


function modeldic(areamodel::CoVAreaModel)
        D = -2*areamodel.MCMC_results.logL
        D̄ = mean(D)
        p_D = 0.5*var(D)
        return D̄ + p_D
end
