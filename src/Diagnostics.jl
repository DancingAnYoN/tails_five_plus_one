module Diagnostics
using ..Types: Config, IA, IB, IPHI, IPI, IFIELD
using ..FiniteDifferences: deriv4!
using ..Geometry: misner_sharp_mass,physical_mass
using ..Constraints: reconstruct_A!, hamiltonian_residual
export diagnostics
"""Six-dimensional mass aspect mu (not ADM mass) and residual diagnostics.
For fixed-background runs the Hamiltonian residual includes the *test-field*
stress tensor by design; it is not an evolution constraint of that mode.
"""
function diagnostics(u::Matrix{Float64},f::Float64,r::Vector{Float64},p::Config)
    n=length(r)
    da=zeros(n);dB=zeros(n);dphi=zeros(n);A=zeros(n)
    @views begin
        deriv4!(da,u[:,IA],p.dr)
        deriv4!(dB,u[:,IB],p.dr)
        deriv4!(dphi,u[:,IFIELD],p.dr)
    end
    J=r .+ f
    reconstruct_A!(A,u,J,dB,p.G6)
    hsum=0.0; defs=0.0; count=0
    @inbounds for i in 3:n-2
        h=hamiltonian_residual(J[i],u[i,IA],da[i],u[i,IB],A[i],
                               u[i,IPHI],u[i,IPI],p.G6)
        hsum+=h*h
        defs+=(u[i,IPHI]-dphi[i])^2
        count+=1
    end
    jh=J[1];ah=u[1,IA]; bh=u[1,IB]
    muout=misner_sharp_mass(J[end],u[end,IA],u[end,IB])
    muh=jh^3/2
    return (ham_l2=sqrt(hsum/max(count,1)),
            scalar_constraint_l2=sqrt(defs/max(count,1)),
            horizon_error=ah*jh*bh-1,
            mu_h=muh,mu_out=muout,
            M_h=physical_mass(muh,p.G6), M_out=physical_mass(muout,p.G6),
            horizon_radius=jh)
end
end
