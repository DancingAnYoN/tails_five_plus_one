module Geometry
using ..Types: Config, IA, IPHI, IPI
export schwarzschild_fields, schwarzschild_Bprime, fdot_horizon,
       lapse_shift, misner_sharp_mass, physical_mass, horizon_radius, initial_mu

"""5+1 Schwarzschild--Tangherlini in the ingoing t=v-J gauge.
q=(rh/J)^3, a=sqrt(1+q), B=q/(J*a), A=-3q(1+q/2)/(J*a^3).
"""
function schwarzschild_fields(J::Real,rh::Real)
    q=(rh/J)^3
    a=sqrt(1+q)
    B=q/(J*a)
    A=-3q*(1+q/2)/(J*a^3)
    beta=q/(1+q)
    return (a=a,B=B,A=A,beta=beta)
end
@inline function schwarzschild_Bprime(J::Real,rh::Real)
    q=(rh/J)^3
    return -q*(4+5q/2)/(J^2*(1+q)^(3/2))
end
@inline initial_mu(p::Config)=p.rh0^3/2
@inline horizon_radius(p::Config,f::Real)=p.rh0+f

function lapse_shift(a::Real,J::Real,B::Real,fd::Real)
    q=a*J*B
    1+q>0 || error("MMIEF denominator 1+aJB <= 0")
    beta=(fd+q)/(1+q)
    alpha=a*(1-beta)
    alpha>0 || error("Nonpositive MMIEF lapse")
    return (alpha=alpha,beta=beta)
end

"""Geometrical mass aspect mu (dimensions length^3): F=1-2mu/J^3."""
@inline misner_sharp_mass(J,a,B)=J^3/2*(1-inv(a*a)+(J*B)^2)
"""Physical quasi-local mass M = 4pi mu/(3 G6), not equal to mu."""
@inline physical_mass(mu,G6)=4π*mu/(3G6)

"""Six-dimensional horizon-locking law evaluated at grid index one."""
@inline function fdot_horizon(u::Matrix{Float64},f::Float64,p::Config)
    if p.fixed_background
        return 0.0
    end
    J=horizon_radius(p,f)
    return (2π*p.G6/3)*(J/u[1,IA])^2*(u[1,IPHI]+u[1,IPI])^2
end
end
