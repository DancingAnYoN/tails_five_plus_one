module Evolution
using Printf
using ..Types: Config,IA,IB,IPHI,IPI,IFIELD
using ..Grid: radial_grid
using ..FiniteDifferences: deriv4!,ko6!
using ..Geometry: schwarzschild_fields,fdot_horizon,initial_mu
using ..Constraints: project_hamiltonian!,reconstruct_A!
using ..ScalarField: scalar_rhs!
using ..EinsteinEquations: gravity_rhs!
using ..InitialData: initial_data
using ..BoundaryConditions: outer_projection!,outer_rhs!
using ..Integrators: RK4Workspace
using ..MethodOfLines: MOLSystem,mol_step!,cfl_dt
using ..Observers: sample_fixed_J
using ..Diagnostics: diagnostics
using ..DataOutput: RunData,store_sample!,store_frame!
export run_evolution,vacuum_test,rhs!,Workspace,setup_vacuum_balance!,
       evolve_one_step!,EvolutionWorkspace

mutable struct Workspace
    r::Vector{Float64}
    J::Vector{Float64}
    A::Vector{Float64}
    beta::Vector{Float64}
    dAmetric::Vector{Float64}
    dB::Vector{Float64}
    dPhi::Vector{Float64}
    dBeta::Vector{Float64}
    dABeta::Vector{Float64}
    dFFlux::Vector{Float64}
    dPFlux::Vector{Float64}
    aBeta::Vector{Float64}
    fFlux::Vector{Float64}
    pFlux::Vector{Float64}
    vacuum_a_rhs::Vector{Float64}
    vacuum_B_rhs::Vector{Float64}
end
function Workspace(r::Vector{Float64})
    x=[zeros(length(r)) for _ in 1:15]
    return Workspace(r,x...)
end

"""Balance discrete stationary Schwarzschild--Tangherlini source at f=0.
This is a numerical well-balancing operation, not a proof of consistency.
"""
function setup_vacuum_balance!(p::Config,w::Workspace)
    n=length(w.r)
    uv=zeros(n,5); dv=similar(uv)
    @inbounds for i in 1:n
        s=schwarzschild_fields(w.r[i],p.rh0)
        uv[i,IA]=s.a; uv[i,IB]=s.B
    end
    rhs!(dv,uv,0.0,p,w;subtract_background=false)
    @views w.vacuum_a_rhs .= dv[:,IA]
    @views w.vacuum_B_rhs .= dv[:,IB]
    return nothing
end

"""Five 5+1 coupled MMIEF field right-hand sides and horizon-shift ODE.
If fixed_background=true, a,B,f are held at their exact vacuum values.
"""
function rhs!(du::Matrix{Float64},u::Matrix{Float64},f::Float64,
              p::Config,w::Workspace;subtract_background::Bool=true)
    n=length(w.r);h=p.dr
    @views begin
        a=u[:,IA];B=u[:,IB];Phi=u[:,IPHI];Pi=u[:,IPI]
        deriv4!(w.dB,B,h)
        deriv4!(w.dAmetric,a,h)
        deriv4!(w.dPhi,Phi,h)
        fd=fdot_horizon(u,f,p)
        fd<1 || error(@sprintf("MMIEF lapse collapse: fdot=%.6g",fd))
        @inbounds for i in 1:n
            Ji=w.r[i]+f;ai=a[i]
            (isfinite(ai) && ai>0 && Ji>0) || error("Invalid geometry at index $i")
            w.J[i]=Ji
        end
        reconstruct_A!(w.A,u,w.J,w.dB,p.G6)
        @inbounds for i in 1:n
            Ji=w.J[i]; ai=a[i]; Bi=B[i]
            q=ai*Ji*Bi; den=1+q
            den>0 || error("Invalid MMIEF denominator at index $i")
            b=(fd+q)/den
            w.beta[i]=b
            w.aBeta[i]=ai*b
            w.fFlux[i]=b*Phi[i]+(1-b)*Pi[i]
            w.pFlux[i]=Ji^4*(b*Pi[i]+(1-b)*Phi[i])
        end
        deriv4!(w.dBeta,w.beta,h)
        deriv4!(w.dABeta,w.aBeta,h)
        deriv4!(w.dFFlux,w.fFlux,h)
        deriv4!(w.dPFlux,w.pFlux,h)
        fill!(du,0.0)
        if !p.fixed_background
            gravity_rhs!(du,u,w.J,w.A,w.beta,w.dB,w.dBeta,w.dABeta)
        end
        scalar_rhs!(du,u,w.J,fd,w.dFFlux,w.dPFlux,w.fFlux)
        ko6!(du,u,h,p.ko_strength)
        outer_rhs!(du,u,w.J,w.beta,fd,w.dAmetric,w.dB,w.dPhi,w.dBeta)
        if p.fixed_background
            @views fill!(du[:,IA],0.0)
            @views fill!(du[:,IB],0.0)
        elseif subtract_background
            @inbounds for i in 1:n
                du[i,IA]-=w.vacuum_a_rhs[i]
                du[i,IB]-=w.vacuum_B_rhs[i]
            end
        end
        return fd
    end
end

struct EvolutionWorkspace
    rk4::RK4Workspace
    rhswork::Workspace
end
EvolutionWorkspace(u::Matrix{Float64},r::Vector{Float64})=
    EvolutionWorkspace(RK4Workspace(u),Workspace(r))

function evolve_one_step!(u::Matrix{Float64},f::Float64,dt::Float64,
                          p::Config,w::EvolutionWorkspace)
    rhs_eval! = (du,x,shift) -> rhs!(du,x,shift,p,w.rhswork)
    project! = (x,shift) -> begin
        if !p.fixed_background
            project_hamiltonian!(x,shift,p,w.rhswork.r)
        end
        nothing
    end
    boundary! = (x,shift) -> outer_projection!(x,shift,p,w.rhswork.r)
    preproj=(x,shift)->x[1,IA]*(w.rhswork.r[1]+shift)*x[1,IB]-1
    model=MOLSystem(rhs_eval!,project!,boundary!,preproj)
    return mol_step!(u,f,dt,w.rk4,model)
end

function run_evolution(p::Config;verbose::Bool=true)
    r,u=initial_data(p)
    f=0.0
    outer_projection!(u,f,p,r)
    w=EvolutionWorkspace(u,r)
    if !p.fixed_background
        setup_vacuum_balance!(p,w.rhswork)
        project_hamiltonian!(u,f,p,r)
        outer_projection!(u,f,p,r)
    end
    d0=diagnostics(u,f,r,p)
    initial_outer_mu=d0.mu_out
    mu_budget_tol=max(p.mass_budget_fraction*initial_mu(p),
                      0.1*max(0.0,initial_outer_mu-initial_mu(p)))
    data=RunData(p,r,u)
    store_sample!(data,0.0,u,f)
    store_frame!(data,0.0,u,f)
    if verbose
        @printf("Initial: points=%d, mu_out=%.10g, H_L2=%.4g, horizon error=%.4g\n",
                length(r),d0.mu_out,d0.ham_l2,d0.horizon_error)
    end
    t=0.0;step=0
    next_sample=p.output_dt
    next_frame=p.snapshot_dt
    next_report=25*p.rh0
    while t < p.tmax-1e-13
        fd=fdot_horizon(u,f,p)
        maxspeed=1.0
        @inbounds for i in eachindex(r)
            q=u[i,IA]*(r[i]+f)*u[i,IB]
            b=(fd+q)/(1+q)
            maxspeed=max(maxspeed,abs(1-2b))
        end
        dt=cfl_dt(p.dr,p.cfl,maxspeed,p.tmax-t)
        f,preproj_error=evolve_one_step!(u,f,dt,p,w)
        t+=dt;step+=1
        if !p.fixed_background && p.check_mass_budget
            muh=(r[1]+f)^3/2
            if !(isfinite(muh) && muh<=initial_outer_mu+mu_budget_tol)
                error(@sprintf("Unphysical 5+1 horizon mass-aspect growth at t/rh0=%.6g; mu_h=%.6g, mu_outer_initial=%.6g, f=%.6g. Check CFL, constraints and boundaries.",
                               t/p.rh0,muh,initial_outer_mu,f))
            end
        end
        for Jobs in p.observers
            x=Jobs-f
            if !(isfinite(x) && r[1]<=x<=r[end])
                sample_fixed_J(u,r,f,Jobs,p;time=t)
            end
        end
        if step%20==0
            if !all(isfinite,u) || !isfinite(f) || minimum(@view(u[:,IA]))<=0
                error("Non-finite or non-positive metric detected at t=$t")
            end
            fdot_horizon(u,f,p)<1 || error("MMIEF lapse factor 1-fdot <=0")
        end
        if t+1e-12>=next_sample || t>=p.tmax-1e-12
            store_sample!(data,t,u,f,preproj_error)
            next_sample+=p.output_dt
        end
        if t+1e-12>=next_frame || t>=p.tmax-1e-12
            store_frame!(data,t,u,f)
            next_frame+=p.snapshot_dt
        end
        if verbose && t+1e-12>=next_report
            d=diagnostics(u,f,r,p)
            @printf("t/rh0=%7.2f  mu_h=%.8g  mu_out=%.8g  fdot=%8.2e H_L2=%.3e\n",
                    t/p.rh0,d.mu_h,d.mu_out,fdot_horizon(u,f,p),d.ham_l2)
            next_report+=25*p.rh0
        end
    end
    data.final_state=copy(u)
    return data
end

"""Return max discrete field residual for zero-field Tangherlini initial state."""
function vacuum_test(;dr::Float64=0.25,rmax::Float64=100.0)
    p=Config(dr=dr,rmax=rmax,tmax=1.0,pulse_amplitude=0.0)
    r=radial_grid(p)
    u=zeros(length(r),5)
    for i in eachindex(r)
        s=schwarzschild_fields(r[i],p.rh0)
        u[i,IA]=s.a;u[i,IB]=s.B
    end
    du=similar(u)
    w=Workspace(r)
    setup_vacuum_balance!(p,w)
    rhs!(du,u,0.0,p,w)
    return [maximum(abs,@view(du[:,j])) for j in 1:5]
end
end
