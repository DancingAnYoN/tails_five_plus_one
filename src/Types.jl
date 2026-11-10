module Types
export Config, IA, IB, IPHI, IPI, IFIELD, SimulationState, fieldview

# Six-dimensional spacetime: five spatial dimensions; SO(5) spherical orbits S^4.
const D_SPACIAL = 5
const SPHERE_DIM = D_SPACIAL - 1 # k=4

Base.@kwdef struct Config
    # Horizon radius, NOT ADM mass. Lengths are in arbitrary consistent units.
    rh0::Float64 = 2.0
    G6::Float64 = 1.0
    rmax::Float64 = 600.0
    dr::Float64 = 0.25
    cfl::Float64 = 0.12
    tmax::Float64 = 400.0
    pulse_r0::Float64 = 20.0
    pulse_sigma::Float64 = 2.0
    pulse_amplitude::Float64 = 1e-5
    cutoff_start::Float64 = 6.0
    cutoff_end::Float64 = 9.0
    ko_strength::Float64 = 0.04
    output_dt::Float64 = 0.5
    snapshot_dt::Float64 = 2.0
    profile_stride::Int = 4
    observers::Vector{Float64} = [10.0, 30.0, 50.0]
    # Decouple the scalar from gravity, holding the exact vacuum metric fixed.
    fixed_background::Bool = false
    check_mass_budget::Bool = true
    mass_budget_fraction::Float64 = 0.01
end

const IA, IB, IPHI, IPI, IFIELD = 1,2,3,4,5
mutable struct SimulationState
    fields::Matrix{Float64}
    f::Float64
end
function fieldview(state::SimulationState,name::Symbol)
    index=Dict(:a=>IA,:B=>IB,:Phi=>IPHI,:Pi=>IPI,:phi=>IFIELD)
    haskey(index,name) || throw(ArgumentError("Unknown state field $name"))
    return view(state.fields,:,index[name])
end
end
