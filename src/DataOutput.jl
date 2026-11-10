module DataOutput
using Printf, TOML
using ..Types: Config, IA, IB, IPHI, IPI, IFIELD
using ..Diagnostics: diagnostics
using ..Observers: sample_fixed_J
export RunData, store_sample!, store_frame!, save_csv, save_run, load_run

const SCHEMA_VERSION = 2
mutable struct RunData
    settings::Config
    r::Vector{Float64}
    times::Vector{Float64}
    signals::Vector{Vector{Float64}}
    shifts::Vector{Float64}
    masses_h::Vector{Float64}
    masses_out::Vector{Float64}
    ham_l2::Vector{Float64}
    scalar_constraint_l2::Vector{Float64}
    horizon_error::Vector{Float64}
    horizon_preprojection::Vector{Float64}
    frame_times::Vector{Float64}
    frame_shift::Vector{Float64}
    frame_profiles::Vector{Vector{Float32}}
    final_state::Matrix{Float64}
end

function RunData(p::Config,r::Vector{Float64},u::Matrix{Float64})
    return RunData(p,r,Float64[],[Float64[] for _ in p.observers],Float64[],
                   Float64[],Float64[],Float64[],Float64[],Float64[],
                   Float64[],Float64[],Float64[],Vector{Float32}[],u)
end

function store_sample!(data::RunData,t::Float64,u::Matrix{Float64},f::Float64,
                       preprojection_error::Float64=0.0)
    # Sample all observers before modifying RunData, so a failed interpolation
    # cannot leave its arrays with inconsistent lengths.
    wave = [sample_fixed_J(u,data.r,f,J,data.settings;time=t)
            for J in data.settings.observers]
    d=diagnostics(u,f,data.r,data.settings)
    push!(data.times,t)
    push!(data.shifts,f)
    push!(data.masses_h,d.M_h)
    push!(data.masses_out,d.M_out)
    push!(data.ham_l2,d.ham_l2)
    push!(data.scalar_constraint_l2,d.scalar_constraint_l2)
    push!(data.horizon_error,d.horizon_error)
    push!(data.horizon_preprojection,preprojection_error)
    for k in eachindex(data.settings.observers)
        push!(data.signals[k],wave[k])
    end
end

function store_frame!(data::RunData,t::Float64,u::Matrix{Float64},f::Float64)
    push!(data.frame_times,t)
    push!(data.frame_shift,f)
    stride=data.settings.profile_stride
    push!(data.frame_profiles, Float32[
        (data.r[i]+f)^2*u[i,IFIELD] for i in 1:stride:length(data.r)])
end
function save_csv(data::RunData,filename::AbstractString)
    open(filename,"w") do io
        cols=["t","T_asymptotic","f","M_h","M_outer","hamiltonian_L2","scalar_definition_L2","horizon_error",
              "horizon_preprojection_error"]
        for J in data.settings.observers
            push!(cols,@sprintf("phi_J%.2f",J))
        end
        println(io,join(cols,","))
        for i in eachindex(data.times)
            vals=[data.times[i],data.times[i]-data.shifts[i],data.shifts[i],data.masses_h[i],data.masses_out[i],
                  data.ham_l2[i],data.scalar_constraint_l2[i],
                  data.horizon_error[i],data.horizon_preprojection[i]]
            append!(vals,[sig[i] for sig in data.signals])
            println(io,join([@sprintf("%.16e",x) for x in vals],","))
        end
    end
    return filename
end

"""Write portable, versioned CSV+TOML run data. Snapshot rows contain
(t,f,r,J,J2phi); rh0 is the initial horizon radius; M is physical 6D mass.
No GLMakie dependency and no global mutable output path.
"""
function save_run(data::RunData, directory::AbstractString)
    mkpath(directory)
    config=Dict(string(k)=>getfield(data.settings,k) for k in fieldnames(Config))
    open(joinpath(directory,"config.toml"),"w") do io
        TOML.print(io,config)
    end
    open(joinpath(directory,"metadata.toml"),"w") do io
        TOML.print(io,Dict(
            "schema_version"=>SCHEMA_VERSION,
            "code"=>"MMIEFD6 v0.1.0",
            "units"=>"c=1; G6 and rh0 in config.toml; masses M=4*pi*mu/(3G6)",
            "time_coordinate"=>"t; T_asymptotic=t-f",
            "grid_coordinate"=>"r; areal radius J=r+f",
            "evolution"=>(data.settings.fixed_background ? "fixed Schwarzschild--Tangherlini test field" : "partially constrained, RK4+Hamiltonian projection"),
            "spacetime_dimension"=>6,
            "snapshot_quantity"=>"J2phi=J^2*phi",
            "column_order"=>"a,B,Phi,Pi,phi"
        ))
    end
    open(joinpath(directory,"diagnostics.csv"),"w") do io
        println(io,"t,T_asymptotic,f,M_h,M_outer,hamiltonian_L2,scalar_definition_L2,horizon_error,horizon_preprojection_error")
        for i in eachindex(data.times)
            fields=(data.times[i],data.times[i]-data.shifts[i],data.shifts[i],
                    data.masses_h[i],data.masses_out[i],data.ham_l2[i],
                    data.scalar_constraint_l2[i],data.horizon_error[i],
                    data.horizon_preprojection[i])
            println(io,join((@sprintf("%.17g",x) for x in fields),","))
        end
    end
    open(joinpath(directory,"observers.csv"),"w") do io
        println(io,"t,T_asymptotic,J_obs,phi")
        for i in eachindex(data.times), j in eachindex(data.signals)
            fields=(data.times[i], data.times[i]-data.shifts[i],
                    data.settings.observers[j], data.signals[j][i])
            println(io,join((@sprintf("%.17g",x) for x in fields),","))
        end
    end
    open(joinpath(directory,"snapshots.csv"),"w") do io
        println(io,"t,f,r,J,J2phi")
        grid=view(data.r,1:data.settings.profile_stride:length(data.r))
        for k in eachindex(data.frame_times)
            profile=data.frame_profiles[k]
            length(profile)==length(grid) || error("Snapshot/grid length mismatch")
            t=data.frame_times[k]; f=data.frame_shift[k]
            for i in eachindex(profile)
                fields=(t,f,grid[i],grid[i]+f,profile[i])
                println(io,join((@sprintf("%.17g",x) for x in fields),","))
            end
        end
    end
    open(joinpath(directory,"final_state.csv"),"w") do io
        println(io,"r,a,B,Phi,Pi,phi")
        for i in eachindex(data.r)
            fields=(data.r[i],data.final_state[i,1],data.final_state[i,2],
                    data.final_state[i,3],data.final_state[i,4],data.final_state[i,5])
            println(io,join((@sprintf("%.17g",x) for x in fields),","))
        end
    end
    return abspath(directory)
end

"""Read the versioned output independent of the simulation process.
Uses only Julia stdlib; no GLMakie or simulation required to load the data.
"""
function load_run(directory::AbstractString)
    meta=TOML.parsefile(joinpath(directory,"metadata.toml"))
    get(meta,"schema_version",nothing)==SCHEMA_VERSION ||
        error("Unsupported MMIEFD6 output schema version")
    raw=TOML.parsefile(joinpath(directory,"config.toml"))
    opts=(Symbol(k)=>v for (k,v) in raw)
    p=Config(;opts...)
    data_rows=readlines(joinpath(directory,"final_state.csv"))
    length(data_rows)>=2 || error("No final state stored")
    r=Float64[]
    u=Matrix{Float64}(undef,length(data_rows)-1,5)
    for (i,line) in enumerate(data_rows[2:end])
        values=parse.(Float64,split(line,','))
        length(values)==6 || error("Invalid final_state.csv at row $i")
        push!(r,values[1])
        u[i,:].=values[2:6]
    end
    data=RunData(p,r,u)
    for line in Iterators.drop(eachline(joinpath(directory,"diagnostics.csv")),1)
        v=parse.(Float64,split(line,','))
        length(v)==9 || error("Invalid diagnostics.csv")
        push!(data.times,v[1]); push!(data.shifts,v[3]); push!(data.masses_h,v[4])
        push!(data.masses_out,v[5]); push!(data.ham_l2,v[6])
        push!(data.scalar_constraint_l2,v[7]); push!(data.horizon_error,v[8])
        push!(data.horizon_preprojection,v[9])
    end
    rows=0
    for line in Iterators.drop(eachline(joinpath(directory,"observers.csv")),1)
        v=parse.(Float64,split(line,','))
        length(v)==4 || error("Invalid observers.csv")
        j=findfirst(x->isapprox(x,v[3];atol=1e-10,rtol=0),p.observers)
        j===nothing && error("Unknown observer J=$(v[3])")
        push!(data.signals[j],v[4]); rows+=1
    end
    rows==length(data.times)*length(p.observers) || error("Missing observer rows")
    for line in Iterators.drop(eachline(joinpath(directory,"snapshots.csv")),1)
        v=parse.(Float64,split(line,','))
        length(v)==5 || error("Invalid snapshots.csv")
        if isempty(data.frame_times) || v[1]!=data.frame_times[end]
            push!(data.frame_times,v[1]); push!(data.frame_shift,v[2])
            push!(data.frame_profiles,Float32[])
        end
        push!(data.frame_profiles[end],Float32(v[5]))
    end
    expected=length(1:p.profile_stride:length(r))
    all(length(profile)==expected for profile in data.frame_profiles) ||
        error("Snapshot profile length mismatch")
    data.final_state=copy(u)
    return data
end
end # module DataOutput
