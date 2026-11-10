module Grid
using ..Types: Config
export RadialGrid, radial_grid, check_config
struct RadialGrid
    r::Vector{Float64}
    dr::Float64
end
function check_config(p::Config)
    rh=p.rh0
    p.G6>0 && rh>0 || error("G6 and rh0 must be positive")
    p.dr>0 && p.rmax>rh+12*p.dr || error("Invalid radial grid")
    n=(p.rmax-rh)/p.dr
    isapprox(n,round(n);atol=1e-8,rtol=0) || error("(rmax-rh0)/dr must be an integer")
    0<p.cfl<=0.35 || error("Require 0 < cfl <= 0.35 for this prototype")
    p.tmax>0 && p.output_dt>0 && p.snapshot_dt>0 || error("Invalid time/output settings")
    p.pulse_sigma>0 && p.pulse_amplitude>=0 || error("Invalid Gaussian data")
    rh<p.cutoff_start<p.cutoff_end<p.pulse_r0 || error("Cutoff must lie between horizon and pulse")
    all(J -> rh<J<p.rmax,p.observers) || error("Fixed-J observers must be in the domain")
    0<=p.ko_strength<=0.15 || error("KO strength outside experimental interval")
    p.profile_stride>=1 || error("Invalid profile_stride")
    p.mass_budget_fraction>=0 || error("Invalid mass budget fraction")
    return nothing
end
function radial_grid(p::Config)
    check_config(p)
    n=round(Int,(p.rmax-p.rh0)/p.dr)+1
    return [p.rh0+(i-1)*p.dr for i in 1:n]
end
RadialGrid(p::Config)=RadialGrid(radial_grid(p),p.dr)
end
