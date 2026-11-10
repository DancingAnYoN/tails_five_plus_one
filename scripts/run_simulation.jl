# Julia --project=. scripts/run_simulation.jl --smoke --vacuum
# Julia --project=. scripts/run_simulation.jl --standard
# Julia --project=. scripts/run_simulation.jl --standard --fixed
using MMIEFD6
using Printf

function configuration(mode::String,vacuum::Bool,fixed::Bool)
    amp=vacuum ? 0.0 : 1e-5
    if mode=="--smoke"
        return Config(rmax=70.0,dr=0.5,cfl=0.10,tmax=2.0,output_dt=0.25,
                      snapshot_dt=0.25,profile_stride=2,ko_strength=0.03,
                      pulse_amplitude=amp,fixed_background=fixed)
    elseif mode=="--quick"
        return Config(rmax=160.0,dr=0.5,cfl=0.10,tmax=65.0,output_dt=0.25,
                      snapshot_dt=1.0,profile_stride=2,ko_strength=0.03,
                      pulse_amplitude=amp,fixed_background=fixed)
    elseif mode=="--standard"
        return Config(rmax=700.0,dr=0.25,cfl=0.10,tmax=350.0,output_dt=0.5,
                      snapshot_dt=2.0,profile_stride=4,ko_strength=0.03,
                      pulse_amplitude=amp,fixed_background=fixed)
    elseif mode=="--long"
        return Config(rmax=1700.0,dr=0.25,cfl=0.10,tmax=800.0,output_dt=0.5,
                      snapshot_dt=4.0,profile_stride=8,ko_strength=0.03,
                      pulse_amplitude=amp,fixed_background=fixed)
    else
        error("Unknown mode $mode")
    end
end

function main(args::Vector{String}=ARGS)
    if "--help" in args || "-h" in args
        println("""
        5+1 MMIEF Einstein--Klein--Gordon prototype
        julia --project=. scripts/run_simulation.jl [--smoke|--quick|--standard|--long] [--vacuum] [--fixed]
        --fixed : test scalar on exact fixed Tangherlini vacuum (no backreaction)
        --vacuum: zero scalar amplitude (vacuum invariance)
        --standard: full BH+scalar evolution, default --quick
        CSV/TOML output in results/<mode>[_fixed][_vacuum]
        """)
        return nothing
    end
    valid=Set(["--smoke","--quick","--standard","--long","--vacuum","--fixed"])
    all(x->x in valid,args) || error("Unknown flag; see --help")
    modes=[x for x in args if x in ("--smoke","--quick","--standard","--long")]
    length(modes)<=1 || error("Specify exactly one size mode")
    mode=isempty(modes) ? "--quick" : only(modes)
    fixed="--fixed" in args
    vac="--vacuum" in args
    p=configuration(mode,vac,fixed)
    tag=mode[3:end]*(fixed ? "_fixed" : "")*(vac ? "_vacuum" : "")
    target=normpath(joinpath(@__DIR__,"..","results",tag))
    @printf("Running %s | 5+1 | dr=%.4g | tmax=%.4g | rh0=%.4g | fixed=%s\n",
            tag,p.dr,p.tmax,p.rh0,string(fixed))
    data=run_evolution(p)
    output=save_run(data,target)
    @printf("Saved %s | final physical M_h=%.9g | final mu_h=%.9g\n",
            output,data.masses_h[end],data.masses_h[end]*3*p.G6/(4π))
    return data
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main(ARGS)
end
