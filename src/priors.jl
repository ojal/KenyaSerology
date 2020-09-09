function basic_prior_contactrate_PCR(θ)
        E₀,I₀,R,α,p_test = θ
        return  logpdf(Gamma(2,0.5),I₀/3*E₀)
end

function basic_prior_contactrate_PCR_OR(θ)
        E₀,I₀,R,α,p_test,χ = θ
        return  logpdf(Gamma(2,0.5),I₀/3*E₀) + logpdf(Normal(0.,2.),log(χ))
end

function basic_prior_contactrate_PCR_Peff(θ)
        E₀,I₀,R,α,p_test,P_eff = θ
        return logpdf(Gamma(2,1000/2),E₀ + I₀) + logpdf(Beta(2,1),P_eff)
end

# plot(Gamma(10,1000/4.3e7))
# plot(Gamma(5,2*250/1.2e7))
# plot(Gamma(2,500/2))
function basic_prior_contactrate_PCR_Peff_mombasa(θ)
        E₀,I₀,R,α,p_test,P_eff = θ
        return logpdf(Gamma(1,100/1),E₀ + I₀) +
                logpdf(Beta(2,1),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(5,2*250/1.2e7),p_test)+
                logpdf(Gamma(3,0.5/3),α)
end


function basic_prior_contactrate_PCR_Peff_nairobi(θ)
        E₀,I₀,R,α,p_test,P_eff = θ
        return logpdf(Gamma(1,100/1),E₀ + I₀) +
                logpdf(Beta(2,1),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(10,3000/4.3e7),p_test) +
                logpdf(Gamma(3,0.5/3),α)
end

function basic_prior_contactrate_PCR_Peff_semiurban(θ)
        E₀,I₀,R,α,p_test,P_eff = θ
        return logpdf(Gamma(1,10/1),E₀ + I₀) +
                logpdf(Beta(2,1),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(1,1e-4),p_test) +
                logpdf(Gamma(3,0.5/3),α)
end

function basic_prior_contactrate_PCR_Peff_rural(θ)
        E₀,I₀,R,α,p_test,P_eff = θ
        return logpdf(Gamma(1,1/1),E₀ + I₀) +
                logpdf(Beta(2,1),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(1,1e-4),p_test) +
                logpdf(Gamma(3,0.5/3),α)
end
function basic_prior_nairobi(θ)
        E₀,I₀,R,α,p_test,P_eff = θ
        return logpdf(Gamma(1,100/1),E₀ + I₀) +
                logpdf(Beta(2,1),P_eff) +
                logpdf(Gamma(2,2.5/2),R) +
                logpdf(Gamma(10,3000/4.3e7),p_test) +
                logpdf(Gamma(3,0.5/3),α)
end
