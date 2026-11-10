module MMIEFD6
include("Types.jl")
include("Grid.jl")
include("FiniteDifferences.jl")
include("Geometry.jl")
include("ScalarField.jl")
include("EinsteinEquations.jl")
include("Constraints.jl")
include("InitialData.jl")
include("BoundaryConditions.jl")
include("Integrators.jl")
include("MethodOfLines.jl")
include("Observers.jl")
include("Diagnostics.jl")
include("DataOutput.jl")
include("Evolution.jl")
include("Analysis.jl")
include("TailExtraction.jl")
using .Types: Config, SimulationState
using .Grid: radial_grid, RadialGrid
using .Geometry: schwarzschild_fields, misner_sharp_mass, fdot_horizon, physical_mass
using .InitialData: initial_data
using .Diagnostics: diagnostics
using .DataOutput: RunData, save_run, load_run, save_csv
using .Evolution: run_evolution, vacuum_test
using .Analysis: fit_qnm, local_power_index, fit_tail
using .TailExtraction: fit_power_law, predict_power_law, amplitude_envelope,
                       scan_tail_windows, observer_tail_fits, compare_resolutions,
                       save_tail_fits, save_window_scan, save_convergence
export Config, SimulationState, radial_grid, RadialGrid, schwarzschild_fields,
       misner_sharp_mass, physical_mass, fdot_horizon, initial_data,
       diagnostics, RunData, save_run, load_run, save_csv,
       run_evolution, vacuum_test, fit_qnm, local_power_index, fit_tail,
       fit_power_law, predict_power_law, amplitude_envelope,
       scan_tail_windows, observer_tail_fits, compare_resolutions,
       save_tail_fits, save_window_scan, save_convergence
end
