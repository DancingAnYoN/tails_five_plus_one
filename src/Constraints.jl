module Constraints
using ..Types: Config, IA,IB,IPHI,IPI
export reconstruct_A!, hamiltonian_residual, mass_rhs, project_hamiltonian!

"""Radial momentum constraint: A=B+JB'-(2piG6 J/a)Phi Pi."""
function reconstruct_A!(A,u,J,dB,G6::Float64)
    @inbounds for i in eachindex(J)
        A[i]=u[i,IB]+J[i]*dB[i]-2π*G6*J[i]*u[i,IPHI]*u[i,IPI]/u[i,IA]
    end
    return A
end

"""The radial form of the six-dimensional Hamiltonian constraint."""
@inline function hamiltonian_residual(J,a,ap,B,A,Phi,Pi,G6)
    return ap+3*(a^3-a)/(2J)+a^3*J*(A*B+3B^2/2)-π*G6*a*J*(Phi^2+Pi^2)
end

"""Mass constraint mu' in 5+1; J'=1, F=1-2mu/J^3+(JB)^2."""
@inline function mass_rhs(J::Float64,mu::Float64,B::Float64,
                          Phi::Float64,Pi::Float64,G6::Float64)
    ainvsq=1-2mu/J^3+(J*B)^2
    (isfinite(ainvsq) && ainvsq>0) ||
        error("Nonpositive projected 1/a² at J=$J (mu=$mu)")
    ainv=sqrt(ainvsq)
    return 2π*G6*J^4*(0.5*(Phi^2+Pi^2)*ainvsq+J*B*ainv*Phi*Pi)
end

"""Constrained evolution: radial RK4 integration of mu' with fixed horizon mass.
Reconstruct a from the Misner--Sharp relation. No inward integration or
clamping of potentially unphysical inverse-metric values.
"""
function project_hamiltonian!(u::Matrix{Float64},f::Float64,p::Config,r::Vector{Float64})
    n=length(r)
    mu=(p.rh0+f)^3/2
    h=p.dr
    @inbounds for i in 1:n
        J=r[i]+f
        if i>1
            j=i-1
            J0=r[j]+f
            B0,P0,Pi0=u[j,IB],u[j,IPHI],u[j,IPI]
            B1,P1,Pi1=u[i,IB],u[i,IPHI],u[i,IPI]
            if 2<=j<=n-2
                # Centered cubic midpoint interpolation.
                Bmid=(-u[j-1,IB]+9B0+9B1-u[j+2,IB])/16
                Pmid=(-u[j-1,IPHI]+9P0+9P1-u[j+2,IPHI])/16
                Pimid=(-u[j-1,IPI]+9Pi0+9Pi1-u[j+2,IPI])/16
            else
                Bmid=(B0+B1)/2
                Pmid=(P0+P1)/2
                Pimid=(Pi0+Pi1)/2
            end
            k1=mass_rhs(J0,mu,B0,P0,Pi0,p.G6)
            k2=mass_rhs(J0+h/2,mu+h*k1/2,Bmid,Pmid,Pimid,p.G6)
            k3=mass_rhs(J0+h/2,mu+h*k2/2,Bmid,Pmid,Pimid,p.G6)
            k4=mass_rhs(J,mu+h*k3,B1,P1,Pi1,p.G6)
            mu+=h*(k1+2k2+2k3+k4)/6
        end
        inva2=1-2mu/J^3+(J*u[i,IB])^2
        (isfinite(inva2) && inva2>0) || error("Hamiltonian projection failed at J=$J")
        u[i,IA]=inv(sqrt(inva2))
    end
    return nothing
end
end
