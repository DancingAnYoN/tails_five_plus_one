module MethodOfLines
using ..Integrators: RK4Workspace, rk4_step!
export MOLSystem, mol_step!, cfl_dt

"""The semi-discrete spatial RHS and numerical constraints/boundary hooks.
No dependency on any particular PDE, gravitational geometry, or plotting.
"""
struct MOLSystem{R,P,B,D}
    rhs!::R
    project!::P
    boundary!::B
    preprojection::D
end

"""Evolve one time step of an arbitrary method-of-lines semidiscretization."""
function mol_step!(u::Matrix{Float64},f::Float64,dt::Float64,
                   work::RK4Workspace,model::MOLSystem)
    return rk4_step!(u,f,dt,work;
                     rhs! = model.rhs!, project! = model.project!,
                     boundary! = model.boundary!, preprojection=model.preprojection)
end

"""Conservative CFL selection with final-step truncation."""
function cfl_dt(dr::Real,C::Real,maxspeed::Real,remaining::Real)
    dr>0 && C>0 && maxspeed>0 && remaining>0 ||
        throw(ArgumentError("Invalid CFL inputs"))
    return min(C*dr/maxspeed,remaining)
end
end # module MethodOfLines
