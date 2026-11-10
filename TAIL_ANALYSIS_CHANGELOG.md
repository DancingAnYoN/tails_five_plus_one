# Dedicated 5+1 late-time tail analysis — changes

The source-based MMIEF5+1 repository's existing numerical evolution equations,
initial data, RK4 integrator, output schema and original visualization functions
are **unchanged**. This update adds an independent analysis and plotting layer.

## New files

- `src/TailExtraction.jl` — dimensionally normalized power-law fits, optional
  running-maximum amplitude envelope, window scans, three-level self-convergence,
  and CSV writers (does not require GLMakie).
- `viz/TailVisualization.jl` — fit-overlay/LPI/residual figure, fitting-window
  sensitivity figure, three-resolution difference/convergence figure.
- `scripts/repl_tail_analysis.jl` — REPL launcher for these functions.
- `docs/LATE_TIME_TAIL_EXTRACTION.md` — detailed VS Code Julia REPL workflow,
  including scripts to generate independent coarse/medium/fine runs.

## Changed files

- `src/MMIEFD6.jl` — includes and exports the new core module.
- `scripts/repl_visualizations.jl` — adds `tail_diagnostics_run` and
  `tail_convergence_runs` wrappers.
- `test/runtests.jl` — known power-law and analytic second-order convergence
  synthetic tests, plus CSV tests.
- `README.md`, `docs/VISUALIZATION_GUIDE.md`, `VALIDATION.md` — usage and
  limitations of the new analysis.

## Generated files (after using the REPL)

- `tail_extraction_5plus1.png`
- `tail_fit_summary.csv`
- `tail_local_power_index.csv`
- `tail_window_stability_5plus1.png`
- `tail_window_scan.csv`
- `tail_convergence_5plus1.png`
- `tail_convergence_summary.csv`
- `tail_convergence_curves.csv`

**Verification limits:** source files were reviewed and the output ZIP integrity
checked, but Julia and GLMakie are not installed in this environment; source
execution and graphical rendering have not been performed here. Even a passing
unit test is not a convergence-validation of Price tails in the EKG solver.
