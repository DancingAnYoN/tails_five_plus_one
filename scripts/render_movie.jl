# julia --project=viz scripts/render_movie.jl results/quick
import MMIEFD6
include(joinpath(@__DIR__,"..","viz","Visualization.jl"))
import .MMIEFD6Visualization
function main(args::Vector{String}=ARGS)
    isempty(args) && error("Usage: render_movie.jl RUN_DIR [--show]")
    all(x->x=="--show",args[2:end]) || error("Unknown option")
    run=MMIEFD6.load_run(abspath(args[1]))
    return MMIEFD6Visualization.make_pulse_movie(run;outdir=abspath(args[1]),
                                                 show_window=("--show" in args))
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main(ARGS)
end
