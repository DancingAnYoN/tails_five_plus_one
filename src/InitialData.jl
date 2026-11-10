module InitialData
using ..Types: Config, IA,IB,IPHI,IPI,IFIELD
using ..Grid: radial_grid
using ..Geometry: schwarzschild_fields, schwarzschild_Bprime, initial_mu
using ..Constraints: mass_rhs
using ..ScalarField: gaussian
export initial_data

"""Approximate ingoing compact initial wave, with 5+1 far-field r^(-2) scaling.
This Sommerfeld condition is asymptotic, not exact for the 6D wave equation.
"""
@inline function ingoing_Pi(Phi,phi,J,a,B)
    return Phi+2*(1+a*J*B)*phi/J
end

"""Integrate the 5+1 radial Hamiltonian and momentum constraints.
The prescribed A(mu,J) is the vacuum ingoing-EF profile with local mu,
not an exact stationary Kerr--Schild ansatz for perturbed geometry.
"""
function initial_ode(J::Float64,dmu::Float64,dB::Float64,p::Config)
    mu=initial_mu(p)+dmu
    B=schwarzschild_fields(J,p.rh0).B+dB
    inva2=1-2mu/J^3+(J*B)^2
    (isfinite(inva2) && inva2>0) || error("Invalid initial inverse a² at J=$J")
    a=inv(sqrt(inva2))
    q=2mu/J^3
    A=-3q*(1+q/2)/(J*(1+q)^(3/2))
    phi,Phi=gaussian(J,p)
    Pi=ingoing_Pi(Phi,phi,J,a,B)
    mup=mass_rhs(J,mu,B,Phi,Pi,p.G6)
    Bp=(A-B)/J+2π*p.G6*Phi*Pi/a-schwarzschild_Bprime(J,p.rh0)
    return mup,Bp
end

function initial_data(p::Config)
    r=radial_grid(p)
    n=length(r)
    u=Matrix{Float64}(undef,n,5)
    dmu=zeros(Float64,n)
    dB=zeros(Float64,n)
    if !p.fixed_background
        @inbounds for i in 1:n-1
            J=r[i]; h=p.dr
            m0,b0=dmu[i],dB[i]
            km1,kb1=initial_ode(J,m0,b0,p)
            km2,kb2=initial_ode(J+h/2,m0+h*km1/2,b0+h*kb1/2,p)
            km3,kb3=initial_ode(J+h/2,m0+h*km2/2,b0+h*kb2/2,p)
            km4,kb4=initial_ode(J+h,m0+h*km3,b0+h*kb3,p)
            dmu[i+1]=m0+h*(km1+2km2+2km3+km4)/6
            dB[i+1]=b0+h*(kb1+2kb2+2kb3+kb4)/6
        end
    end
    @inbounds for i in 1:n
        J=r[i]
        mu=initial_mu(p)+dmu[i]
        B=schwarzschild_fields(J,p.rh0).B+dB[i]
        inva2=1-2mu/J^3+(J*B)^2
        inva2>0 || error("Invalid metric initial data at J=$J")
        a=inv(sqrt(inva2))
        phi,Phi=gaussian(J,p)
        Pi=ingoing_Pi(Phi,phi,J,a,B)
        u[i,IA]=a; u[i,IB]=B; u[i,IPHI]=Phi
        u[i,IPI]=Pi; u[i,IFIELD]=phi
    end
    return r,u
end
end
