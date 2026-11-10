module BoundaryConditions
using ..Types: Config, IA,IB,IPHI,IPI,IFIELD
using ..Geometry: fdot_horizon
export outer_projection!, outer_rhs!

"""Approximate far-field outgoing scalar mode in 5 spatial dimensions.
It uses (d_t+c_out d_r)(J^2 phi) ~= 0; at a finite radius this is
not exact and may cause reflections. Increase rmax and check convergence.
"""
function outer_projection!(u::Matrix{Float64}, f::Float64, p::Config, r::Vector{Float64})
    i=length(r)
    fd=fdot_horizon(u,f,p)
    J=r[i]+f
    q=u[i,IA]*J*u[i,IB]
    1+q>0 || error("Outer gauge denominator <= 0")
    beta=(fd+q)/(1+q)
    omb=1-beta
    omb>0 || error("Outer lapse factor <= 0")
    c=1-2beta
    u[i,IPI]=-u[i,IPHI]-2*(fd+c)*u[i,IFIELD]/(J*omb)
    return nothing
end

"""Heuristic boundary differential conditions for a and B; scalar Phi flux
comes from differentiating the 5+1 asymptotic outgoing condition.
These are NOT constraint-preserving continuum boundary conditions.
"""
function outer_rhs!(du,u,J,beta,fd,dAmetric,dB,dPhi,dBeta)
    i=size(u,1)
    Ji=J[i]; b=beta[i]; c=1-2b
    du[i,IA]=-c*dAmetric[i]-(fd+c)*(u[i,IA]-1)/Ji
    du[i,IB]=-c*dB[i]-4*(fd+c)*u[i,IB]/Ji
    du[i,IPHI]=-c*dPhi[i] - (2*(fd+c)-2Ji*dBeta[i])*u[i,IPHI]/Ji +
                  2*(fd+c+2Ji*dBeta[i])*u[i,IFIELD]/Ji^2
    du[i,IPI]=-du[i,IPHI] # Pi is projected after every RK4 stage
    return nothing
end
end
