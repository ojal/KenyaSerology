
"""
    function (model::CoVAreaModel)(θ)

Callable on CoVAreaModel object. This uses the log_likelihood (or log density if non-flat priors are used) of the object.
"""
function (model::CoVAreaModel)(θ)
        return model.log_likelihood(θ,model)
end


"""
    function inferparameters!(areamodel::CoVAreaModel,samples,trans,q₀)

Infer the unknown transmission and observation parameters for the Kenyan county described by `areamodel`. Inference is performed by drawing `samples` number of replicates from the posterior distribution with log-density
    given in the model. The parameters and their transformation are encoded in `trans`. Initial parameter search begins at `q₀`. Inference output is stored in-place.
"""
function inferparameters!(areamodel::CoVAreaModel,samples,trans,q₀)
        println("Starting parameter inference")
        n = length(q₀)
        ℓ = TransformedLogDensity(trans, areamodel)#transformed log-likelihood
        ∇ℓ = LogDensityProblems.ADgradient(:ForwardDiff, ℓ)#transformed log-likelihood gradient wrt the parameters
        # D = Diagonal([0.1,1.,1.,0.1,0.1,0.1].^2)
        results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇ℓ, samples,
                                initialization = (q = q₀,))
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

"""
    function inferparameters!(areamodel::CoVAreaModel,samples,trans,q₀,D)

Infer the unknown transmission and observation parameters for the Kenyan county described by `areamodel`. Inference is performed by drawing `samples` number of replicates from the posterior distribution with log-density
    given in the model. The parameters and their transformation are encoded in `trans`. Initial parameter search begins at `q₀`, and the initial kinetic energy search begins with the `Diagonal` matrix `D`. Inference output is stored in-place.
"""
function inferparameters!(areamodel::CoVAreaModel,samples,trans,q₀,D::Diagonal)
        println("Starting parameter inference")
        n = length(q₀)
        ℓ = TransformedLogDensity(trans, areamodel)#transformed log-likelihood
        ∇ℓ = LogDensityProblems.ADgradient(:ForwardDiff, ℓ)#transformed log-likelihood gradient wrt the parameters
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


"""
    function inferparameters!(areamodel::CoVAreaModel,samples,trans)

Infer the unknown transmission and observation parameters for the Kenyan county described by `areamodel`. Inference is performed by drawing `samples` number of replicates from the posterior distribution with log-density
    given in the model. The parameters and their transformation are encoded in `trans`. Initial parameter search has a pre-optimisation using the default adaptive differential evolution algorithm implemented by the `BlackBoxOptim` package. Inference output is stored in-place.
"""
function inferparameters!(areamodel::CoVAreaModel,samples,trans)
        println("Searching for a good initial condition")
        f(x) = -transform_logdensity(trans, areamodel,x)
        searchrange = fill((-5.,5.),TransformVariables.dimension(trans))
        res = bboptimize(f; SearchRange = searchrange,PopulationSize=100)
        q₀ = best_candidate(res)
        inferparameters!(areamodel,samples,trans,q₀)
        return nothing
end

"""
    function inferparameters!(areamodel::CoVAreaModel,samples,trans,D::Diagonal)

Infer the unknown transmission and observation parameters for the Kenyan county described by `areamodel`. Inference is performed by drawing `samples` number of replicates from the posterior distribution with log-density
    given in the model. The parameters and their transformation are encoded in `trans`. Initial parameter search has a pre-optimisation using the default adaptive differential evolution algorithm implemented by the `BlackBoxOptim` package, whilst the initial kinetic energy search begins at `D`.
    Inference output is stored in-place.
"""
function inferparameters!(areamodel::CoVAreaModel,samples,trans,D::Diagonal)
        println("Searching for a good initial condition")
        f(x) = -transform_logdensity(trans, areamodel,x)
        searchrange = fill((-5.,5.),TransformVariables.dimension(trans))
        res = bboptimize(f; SearchRange = searchrange,PopulationSize=100)
        q₀ = best_candidate(res)
        inferparameters!(areamodel,samples,trans,q₀,D)
        return nothing
end

"""
    function modeldic(areamodel::CoVAreaModel)

Calculate the deviance information criterion (DIC) for the posterior parameter draws in `areamodel`.
"""
function modeldic(areamodel::CoVAreaModel)
        D = -2*areamodel.MCMC_results.logL
        D̄ = mean(D)
        p_D = 0.5*var(D)
        return D̄ + p_D
end
