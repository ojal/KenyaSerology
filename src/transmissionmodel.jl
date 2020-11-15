"""
make_odeproblemforinfence(;startdate = startdate::Date,changedate = changedate::Date,enddate = enddate::Date)

The default constructor for an ODE Problem to be used in the inference method - two R values
"""
function make_odeproblemforinference(;startdate = startdate::Date,changedate = changedate::Date,enddate = enddate::Date)
    time_change_start = (changedate - startdate).value
    t_end = (enddate - startdate).value
    function two_R(t,R1,R2)
        if t< time_change_start
            return R1
        elseif t>=time_change_start
            return R2
        end
    end
    function f(du,u,p,t)
        S,E,I,R,cum_inf = u #Get current state of model
        R_eff1,R_eff2,σ,γ,N = p #Get parameters from parameter vector
        ι = two_R(t,R_eff1,R_eff2)*γ*S*I/N #Incidence rate
        du[1] = -ι #Loss of susceptibles
        du[2] = ι - σ*E #Flux of exposed
        du[3] = σ*E - γ*I #Flux of infecteds
        du[4] = γ*I #Gain of recovered/immunes
        du[5] = ι #Rate of new infections
    end
    ff = ODEFunction(f,syms = [:S,:E,:I,:R,:C]) #This just augments the vector field function with names for each variable
    #Placeholder parameters and initial conditions
    p = [1.3,2.,1/3.1,1/9.,4.2e6]
    x₀ = [4.2e6,100,100,0,0]
    #Define the ODEProblem for the SEIR model
    prob = ODEProblem(ff,x₀,(0.,t_end),p)
    return prob
end

"""
make_odeproblemforinference(contactrate_data;startdate = startdate::Date,enddate = enddate::Date)

The default constructor for an ODE Problem to be used in the inference method with daily changing contact rates
"""
function make_odeproblemforinference(contactrate_data;startdate = startdate::Date,enddate = enddate::Date)
    contactrate = vec(contactrate_data.contactrate)
    dates = contactrate_data.date
    t_end = (enddate - startdate).value
    function getcontactrate(t)
        return contactrate[max(min(ceil(Int64,t),length(contactrate)),1)]
    end
    function f(du,u,p,t)
        S,E,I,R,cum_inf = u #Get current state of model
        R,σ,γ,N = p #Get parameters from parameter vector
        ι = R*getcontactrate(t)*γ*S*I/N #Incidence rate
        du[1] = -ι #Loss of susceptibles
        du[2] = ι - σ*E #Flux of exposed
        du[3] = σ*E - γ*I #Flux of infecteds
        du[4] = γ*I #Gain of recovered/immunes
        du[5] = ι #Rate of new infections
    end
    ff = ODEFunction(f,syms = [:S,:E,:I,:R,:C]) #This just augments the vector field function with names for each variable
    #Placeholder parameters and initial conditions
    p = [3.,1/3.1,1/9.,4.3e6]
    x₀ = [4.3e6,100,100,0,0]
    #Define the ODEProblem for the SEIR model
    prob = ODEProblem(ff,x₀,(0.,t_end),p)
    return prob
end

"""
    make_odeproblemforinference_parameter_ct(projected_contactrate;startdate::Date,enddate::Date)

Construct an ODE Problem to simulate the county level epidemic, where the effective contact rates are parameters rather than data
"""
function make_odeproblemforinference_parameter_ct(projected_contactrate;startdate::Date,enddate::Date)
    dates = projected_contactrate.date
    t_end = (enddate - startdate).value
    function getcontactrate(t,p)
        contactrate = @view p[5:end]
        return contactrate[max(min(ceil(Int64,t),length(contactrate)),1)]
    end
    function f(du,u,p,t)
        S,E,I,R,cum_inf = u #Get current state of model
        R,σ,γ,N = @view p[1:4] #Get parameters from parameter vector
        ι = R*getcontactrate(t,p)*γ*S*I/N #Incidence rate
        du[1] = -ι #Loss of susceptibles
        du[2] = ι - σ*E #Flux of exposed
        du[3] = σ*E - γ*I #Flux of infecteds
        du[4] = γ*I #Gain of recovered/immunes
        du[5] = ι #Rate of new infections
    end
    ff = ODEFunction(f,syms = [:S,:E,:I,:R,:C]) #This just augments the vector field function with names for each variable
    #Placeholder parameters and initial conditions
    p = [3.,1/3.1,1/9.,4.3e6]
    x₀ = [4.3e6,100,100,0,0]
    #Define the ODEProblem for the SEIR model
    prob = ODEProblem(ff,x₀,(0.,t_end),p)
    return prob
end

function make_odeproblemforinference_simple(contactrate;startdate = startdate::Date,enddate = enddate::Date)
    t_end = (enddate - startdate).value
    function getcontactrate(t)
        return contactrate[max(min(ceil(Int64,t),length(contactrate)),1)]
    end
    function f(du,u,p,t)
        S,E,I,R,cum_inf = u #Get current state of model
        R,σ,γ,N = p #Get parameters from parameter vector
        ι = R*getcontactrate(t)*γ*S*I/N #Incidence rate
        du[1] = -ι #Loss of susceptibles
        du[2] = ι - σ*E #Flux of exposed
        du[3] = σ*E - γ*I #Flux of infecteds
        du[4] = γ*I #Gain of recovered/immunes
        du[5] = ι #Rate of new infections
    end
    ff = ODEFunction(f,syms = [:S,:E,:I,:R,:C]) #This just augments the vector field function with names for each variable
    #Placeholder parameters and initial conditions
    p = [3.,1/3.1,1/9.,4.3e6]
    x₀ = [4.3e6,100,100,0,0]
    #Define the ODEProblem for the SEIR model
    prob = ODEProblem(ff,x₀,(0.,t_end),p)
    return prob
end

"""
function simple_conv(incidence,kern)

Direct implementation of convolution [incidence ⋆ kern](t)
"""
function simple_conv(incidence,kern)
    z = similar(incidence)
    n = length(incidence)
    m = length(kern)
    for t = 1:n
        z[t] = 0.
        for s = 1:(t-1)
            if t-s <= m && t-s >= 1
                z[t] += kern[t-s]*incidence[s]
            end
        end
    end
    return z
end


"""
function get_incidence(sol)
Convert ODE solution into a vector of new infections each day
"""
function get_incidence(sol)
    incidence = [sol.u[t][5] - sol.u[t-1][5]  for t = 2:length(sol.u)]
end
