function basic_prior_contactrate_PCR(θ)
        @unpack E₀,I₀,R,α,p_test = θ
        return  logpdf(Gamma(2,0.5),I₀/3*E₀)
end

function basic_prior_contactrate_PCR_OR(θ)
        @unpack E₀,I₀,R,α,p_test,χ = θ
        return  logpdf(Gamma(2,0.5),I₀/3*E₀) + logpdf(Normal(0.,2.),log(χ))
end

function basic_prior_contactrate_PCR_Peff(θ)
        @unpack E₀,I₀,R,α,p_test,P_eff = θ
        return logpdf(Gamma(2,1000/2),E₀ + I₀) + logpdf(Beta(2,1),P_eff)
end

# plot(Gamma(10,1000/4.3e7))
# plot(Gamma(5,2*250/1.2e7))
# plot(Gamma(2,500/2))
function basic_prior_contactrate_PCR_Peff_mombasa(θ)
        @unpack E₀,I₀,R,α,p_test,P_eff = θ
        return logpdf(Gamma(1,100/1),E₀ + I₀) +
                logpdf(Beta(2,1),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(5,2*250/1.2e7),p_test)+
                logpdf(Gamma(3,0.5/3),α)
end


function basic_prior_contactrate_PCR_Peff_nairobi(θ)
        @unpack E₀,I₀,R,α,p_test,P_eff = θ
        return logpdf(Gamma(1,100/1),E₀ + I₀) +
                logpdf(Beta(2,1),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(10,3000/4.3e7),p_test) +
                logpdf(Gamma(3,0.5/3),α)
end

function basic_prior_contactrate_PCR_Peff_semiurban(θ)
        @unpack E₀,I₀,R,α,p_test,P_eff = θ
        return logpdf(Gamma(1,5/1),E₀) +
                logpdf(Gamma(1,5/1), I₀) +
                logpdf(Beta(2,1),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(1,1e-4),p_test) +
                logpdf(Gamma(3,0.5/3),α)
end

function basic_prior_contactrate_PCR_Peff_rural(θ)
        @unpack E₀,I₀,R,α,p_test,P_eff = θ
        return logpdf(Gamma(1,0.5/1),E₀ ) +
                 logpdf(Gamma(1,0.5/1),I₀) +
                logpdf(Beta(2,1),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(1,1e-4),p_test) +
                logpdf(Gamma(3,0.5/3),α)
end
function basic_prior_nairobi(θ)
        @unpack E₀,I₀,R,α,p_test,P_eff = θ
        return logpdf(Gamma(1,100/1),E₀ + I₀) +
                logpdf(Beta(2,1),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(10,3000/4.3e7),p_test) +
                logpdf(Gamma(3,0.5/3),α)
end

"""
        build_prior_for_nairobi(projected_contactrate_nairobi)

Build a prior for Nairobi which has the log-Normal distribution for the effective contact rate.
"""
function build_prior_for_nairobi(projected_contactrate_nairobi)
        n = length(projected_contactrate_nairobi.contactrate)
        μ = log.(vec(projected_contactrate_nairobi.contactrate))
        Σ = [0.01*exp(-((i-j)^2)/(2*7)) for i = 1:n,j=1:n]
        d_prior_random_ct = MvLogNormal(MvNormal(μ,Σ))
        function basic_prior_nairobi_variable_ct(θ)
                @unpack E₀,I₀,R,α,p_test,P_eff,contactrate = θ
                return logpdf(Gamma(10,50/10),E₀) +
                        logpdf(Gamma(10,50/10),I₀) +
                        logpdf(Beta(2,1),P_eff) +
                        logpdf(Gamma(2,2.5/2),R) +
                        logpdf(Gamma(10,3000/4.3e7),p_test) +
                        logpdf(Gamma(3,0.5/3),α) +
                        logpdf(d_prior_random_ct,contactrate)
                end

        return basic_prior_nairobi_variable_ct
end

function build_prior_for_mombasa_NB_PCR_model(projected_contactrate)
        n = length(projected_contactrate.contactrate)
        μ = log.(vec(projected_contactrate.contactrate))
        Σ = [0.01*exp(-((i-j)^2)/(2*7)) for i = 1:n,j=1:n]
        d_prior_random_ct = MvLogNormal(MvNormal(μ,Σ))
        function basic_prior_nairobi_variable_ct(θ)
                @unpack E₀,I₀,R,α,p_test,P_eff,contactrate = θ
                return logpdf(Gamma(10,50/10),E₀) +
                        logpdf(Gamma(10,50/10),I₀) +
                        logpdf(Beta(2,1),P_eff) +
                        logpdf(Gamma(2,2.5/2),R) +
                        logpdf(Gamma(10,3000/4.3e7),p_test) +
                        logpdf(Gamma(3,0.5/3),α) +
                        logpdf(d_prior_random_ct,contactrate)
                end

        return basic_prior_nairobi_variable_ct
end

function build_prior_for_mombasa_BB_PCR_model(projected_contactrate)
        n = length(projected_contactrate.contactrate)
        μ = log.(vec(projected_contactrate.contactrate))
        Σ = [0.01*exp(-((i-j)^2)/(2*7)) for i = 1:n,j=1:n]
        d_prior_random_ct = MvLogNormal(MvNormal(μ,Σ))
        function basic_prior_nairobi_variable_ct(θ)
                @unpack E₀,I₀,R,α,p_test,P_eff,contactrate = θ
                return logpdf(Gamma(10,50/10),E₀) +
                        logpdf(Gamma(10,50/10),I₀) +
                        logpdf(Beta(2,1),P_eff) +
                        logpdf(Gamma(2,2.5/2),R) +
                        logpdf(Gamma(10,3000/4.3e7),p_test) +
                        logpdf(Gamma(3,0.5/3),α) +
                        logpdf(d_prior_random_ct,contactrate)
                end

        return basic_prior_nairobi_variable_ct
end


function mombasa_denominator_prior_def_negs(θ)
    @unpack E₀,I₀,R,χ,M_PCR,P_eff = θ
        return logpdf(Gamma(1,100. /1.),E₀) +
                logpdf(Gamma(1,100/1),I₀) +
                logpdf(Beta(20,5),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(3,2/3),χ)+
                logpdf(Gamma(2,10/2),M_PCR)
end

function rural_denominator_prior_def_negs(θ)
    @unpack E₀,I₀,R,χ,M_PCR,P_eff = θ
        return logpdf(Gamma(1,0.5/1),E₀ ) +
                 logpdf(Gamma(1,0.5/1),I₀) +
                logpdf(Beta(2,1),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(1,1e-4),p_test) +
                logpdf(Gamma(3,0.5/3),α)
end

function nairobi_denominator_prior(θ)
    @unpack E₀,I₀,R,χ,M_PCR,P_eff,α,p_test = θ
        return logpdf(Gamma(1,100. /1.),E₀) +
                logpdf(Gamma(1,100/1),I₀) +
                logpdf(Beta(20,5),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(3,2/3),χ)+
                logpdf(Gamma(2,10/2),M_PCR)+
                logpdf(Gamma(3,0.5/3),α) +
                logpdf(Gamma(10,3000/4.3e7),p_test)
end

function semiurban_denominator_prior(θ)
        @unpack E₀,I₀,R,χ,M_PCR,P_eff,α,p_test = θ
        return logpdf(Gamma(1,5/1),E₀) +
                logpdf(Gamma(1,5/1), I₀) +
                logpdf(Beta(2,1),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(3,2/3),χ)+
                logpdf(Gamma(2,10/2),M_PCR)+
                logpdf(Gamma(3,0.5/3),α) +
                logpdf(Gamma(1,1e-5),p_test) +
                logpdf(Gamma(3,0.5/3),α)
end

function rural_denominator_prior(θ)
        @unpack E₀,I₀,R,χ,M_PCR,P_eff,α,p_test = θ
        return logpdf(Gamma(1,0.5/1),E₀ ) +
                 logpdf(Gamma(1,0.5/1),I₀) +
                logpdf(Beta(2,1),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(3,2/3),χ)+
                logpdf(Gamma(2,10/2),M_PCR)+
                logpdf(Gamma(1,1e-5),p_test) +
                logpdf(Gamma(3,0.5/3),α)
end

"""
        build_mv_lognormalprior(projected_contactrate,l,σ)

Create a multi-variate log-normal distribution with a daily log-mean `μ = log.(projected_contactrate)`, and a covariance kernel
which is squared-exponential with length scale `l` and std. dev. `σ`.
"""
function build_mv_lognormalprior(projected_contactrate,l,σ)
        n = length(projected_contactrate)
        μ = log.(vec(projected_contactrate)) .- (σ^2/2) #Mean correction
        Σ = [(σ^2)*exp(-((i-j)^2)/(2*l^2)) for i = 1:n,j=1:n]
        d_prior_random_ct = MvLogNormal(MvNormal(μ,Σ .+ Diagonal(0.0001*ones(n)))) #Extra small diagonal term increases numerical stability
        return d_prior_random_ct
end
