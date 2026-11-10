# REPL-only convenience helper. In VS Code Julia REPL, at project root:
# using Pkg; Pkg.activate("viz"); Pkg.develop(path=pwd()); Pkg.instantiate()
# include("scripts/repl_visualizations.jl")
# visualize_run("results/quick_fixed"; fps=24, every=1, radial_max=100.0)
# This does not run the PDE. It reads saved schema-v2 results and saves figures.
import MMIEFD6
include(joinpath(@__DIR__, "..", "viz", "Visualization.jl"))
import .MMIEFD6Visualization

"""Read a saved 5+1 run and export the available GLMakie visualizations.

`movie=false` skips MP4 encoding; `dashboard=false` or `history=false` skips PNGs.
`fps`, `every` and `radial_max` only affect rendering, not the PDE evolution.
Returns a named tuple of paths (`nothing` for skipped products).
"""
function visualize_run(directory::AbstractString="results/quick_fixed";
                       dashboard::Bool=true, history::Bool=true,
                       movie::Bool=true, fps::Int=24, every::Int=1,
                       radial_max::Real=100.0, show_window::Bool=false)
    location=abspath(directory)
    isfile(joinpath(location,"metadata.toml")) ||
        error("No saved run at $location. Run the core solver first.")
    data=MMIEFD6.load_run(location)
    p1=dashboard ? MMIEFD6Visualization.make_dashboard(data;
                      outdir=location,show_window=show_window) : nothing
    p2=history ? MMIEFD6Visualization.make_observer_history(data;
                      outdir=location,show_window=show_window) : nothing
    p3=movie ? MMIEFD6Visualization.make_pulse_movie(data;
                      outdir=location,fps=fps,every=every,radial_max=radial_max,
                      show_window=show_window) : nothing
    return (dashboard=p1, observer_history=p2, movie=p3)
end

# Dedicated late-time tail figures and tables. Keep core analysis independent
# of GLMakie; only these wrappers load viz/TailVisualization.jl.
include(joinpath(@__DIR__, "..", "viz", "TailVisualization.jl"))
import .MMIEFD6TailVisualization

"""From the Julia REPL, generate the power-law overlays, local power indices,
fit-window sensitivity and CSV tables for a saved 5+1 run.

`tail_diagnostics_run("results/standard_fixed"; tlo=180.0, thi=300.0)`

The numeric endpoints are raw observer time T=t-f (NOT T/rh0).
For the optional window scan, starts/widths are in the same units. The default
reference p=10 is drawn only for `fixed_background=true`; no p=6 curve is
automatically claimed for coupled black-hole evolution.
"""
function tail_diagnostics_run(directory::AbstractString;
                              tlo::Real,thi::Real,
                              method::Symbol=:absolute,
                              envelope_halfwidth::Real=0.0,
                              amplitude_floor::Real=0.0,
                              peak_halfwidth::Real=0.0,
                              lpi_halfwidth::Real=12.0,
                              window_scan::Bool=true,
                              scan_starts=nothing,
                              scan_widths=nothing,
                              reference_power=nothing,
                              show_window::Bool=false)
    location=abspath(directory)
    run=MMIEFD6.load_run(location)
    times=run.times .- run.shifts
    !isempty(times) && first(times)<=tlo && last(times)>=thi ||
        throw(ArgumentError("Requested [$tlo,$thi] is outside saved observer times. Choose a longer run or an in-range window."))
    a=MMIEFD6TailVisualization.make_tail_diagnostics(run;
        outdir=location,tlo=tlo,thi=thi,method=method,
        envelope_halfwidth=envelope_halfwidth,
        amplitude_floor=amplitude_floor,peak_halfwidth=peak_halfwidth,
        lpi_halfwidth=lpi_halfwidth,reference_power=reference_power,
        show_window=show_window)
    b=window_scan ? MMIEFD6TailVisualization.make_tail_window_stability(run;
        outdir=location,tlo=tlo,thi=thi,starts=scan_starts,
        widths=scan_widths,method=method,
        envelope_halfwidth=envelope_halfwidth,
        amplitude_floor=amplitude_floor,reference_power=reference_power,
        show_window=show_window) : nothing
    return (diagnostics=a,window_scan=b)
end

"""Compare three matched runs, supplied in coarse/medium/fine dr order.

`tail_convergence_runs(["results/h04", "results/h02", "results/h01"];
                       tlo=180.0, thi=300.0)`

It is essential to use distinct saved runs with the *same physics* and
refinement factor. Merely resampling one run is not a convergence test.
"""
function tail_convergence_runs(directories::AbstractVector{<:AbstractString};
                               tlo::Real,thi::Real,
                               outdir::AbstractString="results/tail_convergence",
                               method::Symbol=:absolute,
                               envelope_halfwidth::Real=0.0,
                               amplitude_floor::Real=0.0,
                               reference_power=nothing,
                               show_window::Bool=false)
    length(directories)==3 || throw(ArgumentError("Give exactly three run directories: coarse, medium, fine"))
    runs=map(path->MMIEFD6.load_run(abspath(path)),directories)
    return MMIEFD6TailVisualization.make_tail_convergence(runs[1],runs[2],runs[3];
        outdir=abspath(outdir),tlo=tlo,thi=thi,method=method,
        envelope_halfwidth=envelope_halfwidth,
        amplitude_floor=amplitude_floor,reference_power=reference_power,
        show_window=show_window)
end
