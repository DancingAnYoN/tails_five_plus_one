# REPL-first launcher for only the tail-analysis visualizations.
# 1. Activate the viz environment as documented in README.md.
# 2. include("scripts/repl_tail_analysis.jl")
# 3. tail_diagnostics_run("results/standard_fixed"; tlo=180.0,thi=300.0)
# 4. tail_convergence_runs(["results/h050", "results/h025", "results/h0125"];
#                          tlo=180.0,thi=300.0)
#
# Do not include scripts/repl_visualizations.jl and this file in the SAME
# REPL session: they both load the visualization modules. Use one or the other.
include(joinpath(@__DIR__, "repl_visualizations.jl"))
