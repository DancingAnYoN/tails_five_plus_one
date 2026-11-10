module ScalarField
using ..Types: Config, IPHI, IPI, IFIELD
export gaussian, scalar_rhs!
"""Gaussian pulse multiplied by a C-infinity cutoff; return (phi, dphi/dr)."""
function gaussian(r::Float64, p::Config)
    q = (r - p.pulse_r0) / p.pulse_sigma
    g = p.pulse_amplitude * exp(-q*q)
    gp = -2*q/p.pulse_sigma * g
    if r <= p.cutoff_start
        return 0.0, 0.0
    elseif r >= p.cutoff_end
        return g, gp
    end
    w = p.cutoff_end-p.cutoff_start
    x = (r-p.cutoff_start)/w
    # Standard smooth bump: chi=exp(-1/x)/(exp(-1/x)+exp(-1/(1-x))).
    e1 = exp(-1/x)
    e2 = exp(-1/(1-x))
    chi = e1 / (e1+e2)
    chip = chi*(1-chi)*(1/x^2 + 1/(1-x)^2)/w
    return chi*g, chi*gp+chip*g
end

"""Scalar-field portion of the semi-discrete wave equations.
Flux arrays and their derivatives are precomputed by the coupled RHS.
"""
function scalar_rhs!(du, u, J, fd, dFFlux, dPFlux, fFlux)
    @inbounds for i in eachindex(J)
        du[i,IPHI] = dFFlux[i]
        du[i,IPI] = dPFlux[i]/J[i]^4 - 4*fd*u[i,IPI]/J[i]
        du[i,IFIELD] = fFlux[i]
    end
    return nothing
end
end # module ScalarField
