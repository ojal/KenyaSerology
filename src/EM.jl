function stochastic_Ct_negmeanLL(ct,model::KenyaSerology.CoVAreaModel)
    model.prob = KenyaSerology.make_odeproblemforinference_simple(ct,#Method for defining the ODE problem underlying the inference
                                                            startdate = Date(2020,2,21),#Don't change from Feb 21st as start date!
                                                            enddate = Date(2020,10,10))
        transformed_results = [(R=model.MCMC_results.chain[k,1,1],
                    E₀ =model.MCMC_results.chain[k,2,1],
                    I₀ =model.MCMC_results.chain[k,3,1],
                    χ = model.MCMC_results.chain[k,4,1],
                    M_PCR = model.MCMC_results.chain[k,5,1],
                    α = model.MCMC_results.chain[k,6,1],
                    p_test = model.MCMC_results.chain[k,7,1],
                    P_eff = model.MCMC_results.chain[k,8,1]) for k = rand(1:size(model.MCMC_results.chain,1))]

                    return -mean([model(transformed_results[i]) for i = 1:length(transformed_results)])
end


function get_transformed_results(model::KenyaSerology.CoVAreaModel)
    transformed_results = [(R=model.MCMC_results.chain[k,1,1],
                            E₀ =model.MCMC_results.chain[k,2,1],
                            I₀ =model.MCMC_results.chain[k,3,1],
                            χ = model.MCMC_results.chain[k,4,1],
                            M_PCR = model.MCMC_results.chain[k,5,1],
                            α = model.MCMC_results.chain[k,6,1],
                            p_test = model.MCMC_results.chain[k,7,1],
                            P_eff = model.MCMC_results.chain[k,8,1]) for k = 1:500:size(model.MCMC_results.chain,1)]
end


function new_Ct_negmeanLL(ct,model::KenyaSerology.CoVAreaModel,transformed_results)
    model.prob = make_odeproblemforinference_simple(ct,#Method for defining the ODE problem underlying the inference
                                                            startdate = Date(2020,2,21),#Don't change from Feb 21st as start date!
                                                            enddate = Date(2020,10,10))
    # return -mean([model(transformed_results[i]) - model.log_priors(transformed_results[i]) for i = 1:length(transformed_results)])
    return -mean([model(transformed_results[i]) for i = 1:length(transformed_results)])
end

function get_model_and_optimise(filename,projected_contact_rate;varname = "model")
    model = load(filename)[varname]
    transformed_results = get_transformed_results(model)
    function cost(x)
        if all(x .> 0)
            return new_Ct_negmeanLL(vcat(projected_contact_rate.contactrate[1:30],x),model,transformed_results)
        else
            return Inf
        end
    end
    n = (Date(2020,9,30) - Date(2020,2,20)).value
    x₀ = projected_contact_rate.contactrate[31:n]
    opt = optimize(cost,x₀,Optim.Options(iterations=10^4,show_trace = true,show_every = 1000))
    ct = opt.minimizer
    smoothed_ct = vcat(projected_contact_rate.contactrate[1:3],[mean(vcat(projected_contact_rate.contactrate[1:30],ct)[(t-3):(t+3)]) for t = 4:(length(vcat(projected_contact_rate.contactrate[1:30],ct))-3)])

    return smoothed_ct,model
end

function get_model_and_optimise(filename,projected_contact_rate,length_scale::Float64,σ::Float64,priorweight::Float64;varname = "model")
    model = load(filename)[varname]
    n = (Date(2020,9,30) - Date(2020,2,20)).value
    x₀ = projected_contact_rate.contactrate[31:n]
    transformed_results = get_transformed_results(model)
    d_ct = build_mv_lognormalprior(projected_contact_rate.contactrate[1:n],length_scale,σ)

    function cost(x)
        if all(x .> 0)
            y = vcat(projected_contact_rate.contactrate[1:30],x)
            return new_Ct_negmeanLL(y,model,transformed_results) - priorweight*logpdf(d_ct,y)
        else
            return Inf
        end
    end

    function cost_bbopt(x)
        if all(x .> 0)
            y = vcat(projected_contact_rate.contactrate[1:30],x)
            return stochastic_Ct_negmeanLL(y,model) - priorweight*logpdf(d_ct,y)
        else
            return Inf
        end
    end
    #Step 1: BlackBoxOptim
    searchrange = [(x*0.5,x*1.5) for x in projected_contact_rate.contactrate[31:n]]
    opt_bb = bboptimize(cost_bbopt; SearchRange = searchrange,PopulationSize=10*n,MaxSteps=100000)
    y  = best_candidate(opt_bb)
    #Step 2: Nelder-mead

    opt = optimize(cost,y,Optim.Options(iterations=10^4,show_trace = true,show_every = 1000))
    ct = opt.minimizer
    smoothed_ct = vcat(projected_contact_rate.contactrate[1:3],[mean(vcat(projected_contact_rate.contactrate[1:30],ct)[(t-3):(t+3)]) for t = 4:(length(vcat(projected_contact_rate.contactrate[1:30],ct))-3)])

    return smoothed_ct,model
end

function optimise_ct_model(model,projected_contact_rate)
    transformed_results = get_transformed_results(model)
    n = (Date(2020,9,30) - Date(2020,2,20)).value
    x₀ = projected_contact_rate.contactrate[31:n]
    function cost(x)
        if all(x .> 0)
            y = vcat(projected_contact_rate.contactrate[1:30],x)
            return new_Ct_negmeanLL(y,model,transformed_results)
        else
            return Inf
        end
    end

    opt = optimize(cost,x₀,Optim.Options(iterations=10^4,show_trace = true,show_every = 1000))
    ct = opt.minimizer
    smoothed_ct = vcat(projected_contact_rate.contactrate[1:3],[mean(vcat(projected_contact_rate.contactrate[1:30],ct)[(t-3):(t+3)]) for t = 4:(length(vcat(projected_contact_rate.contactrate[1:30],ct))-3)])
    return smoothed_ct
    # return vcat(projected_contact_rate.contactrate[1:30],ct)
end


function optimise_ct_model(model,projected_contact_rate,length_scale,σ,priorweight::Float64)
    transformed_results = get_transformed_results(model)
    n = (Date(2020,9,30) - Date(2020,2,20)).value
    x₀ = projected_contact_rate.contactrate[31:n]
    d_ct = build_mv_lognormalprior(projected_contact_rate.contactrate[1:n],length_scale,σ)

    function cost(x)
        if all(x .> 0)
            y = vcat(projected_contact_rate.contactrate[1:30],x)
            return new_Ct_negmeanLL(y,model,transformed_results) - priorweight*logpdf(d_ct,y)
        else
            return Inf
        end
    end
    function cost_bbopt(x)
        if all(x .> 0)
            y = vcat(projected_contact_rate.contactrate[1:30],x)
            return stochastic_Ct_negmeanLL(y,model) - priorweight*logpdf(d_ct,y)
        else
            return Inf
        end
    end
    #Step 1: BlackBoxOptim
    searchrange = [(x*0.5,x*1.5) for x in projected_contact_rate.contactrate[31:n]]
    opt_bb = bboptimize(cost_bbopt; SearchRange = searchrange,PopulationSize=10*n,MaxSteps=100000)
    y  = best_candidate(opt_bb)
    #Step 2: Nelder-mead

    opt = optimize(cost,y,Optim.Options(iterations=10^4,show_trace = true,show_every = 1000))
    ct = opt.minimizer
    smoothed_ct = vcat(projected_contact_rate.contactrate[1:3],[mean(vcat(projected_contact_rate.contactrate[1:30],ct)[(t-3):(t+3)]) for t = 4:(length(vcat(projected_contact_rate.contactrate[1:30],ct))-3)])

    # return smoothed_ct
    return vcat(projected_contact_rate.contactrate[1:30],ct)

end


function EM_optimise(filename,projected_contact_rate,trans;num_steps=3,varname="model")
    smoothed_ct,model = get_model_and_optimise(filename,projected_contact_rate;varname = varname)
    model.prob = make_odeproblemforinference_simple(smoothed_ct,#Method for defining the ODE problem underlying the inference
                                                            startdate = Date(2020,2,21),#Don't change from Feb 21st as start date!
                                                            enddate = Date(2020,10,1))
    model.contactrate_data = copy(smoothed_ct)
    KenyaSerology.inferparameters!(model,10000,trans)
    for i = 1:(num_steps-1)
        smoothed_ct = optimise_ct_model(model,projected_contact_rate)
        model.prob = make_odeproblemforinference_simple(smoothed_ct,#Method for defining the ODE problem underlying the inference
                                                                startdate = Date(2020,2,21),#Don't change from Feb 21st as start date!
                                                                enddate = Date(2020,10,1))
        model.contactrate_data = copy(smoothed_ct)
        KenyaSerology.inferparameters!(model,10000,trans)
    end
    return model
end

function EM_optimise(filename,projected_contact_rate,trans,length_scale::Float64,σ::Float64,priorweight::Float64;num_steps=3,varname="model")
    smoothed_ct,model = get_model_and_optimise(filename,projected_contact_rate,length_scale,priorweight,σ;varname = varname)
    model.prob = make_odeproblemforinference_simple(smoothed_ct,#Method for defining the ODE problem underlying the inference
                                                            startdate = Date(2020,2,21),#Don't change from Feb 21st as start date!
                                                            enddate = Date(2020,10,1))
    model.contactrate_data = copy(smoothed_ct)
    KenyaSerology.inferparameters!(model,10000,trans)
    for i = 1:(num_steps-1)
        smoothed_ct = optimise_ct_model(model,projected_contact_rate,length_scale,σ,priorweight)
        model.prob = make_odeproblemforinference_simple(smoothed_ct,#Method for defining the ODE problem underlying the inference
                                                        startdate = Date(2020,2,21),#Don't change from Feb 21st as start date!
                                                        enddate = Date(2020,10,1))
        model.contactrate_data = copy(smoothed_ct)
        KenyaSerology.inferparameters!(model,10000,trans)
    end
    return model
end
