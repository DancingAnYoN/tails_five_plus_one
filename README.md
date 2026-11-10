# Late-time tails in 5+1 Einstein--Klein--Gordon (MMIEF)

**Research prototype · Julia 1.10+ (developed for VS Code/Julia 1.12) · 5 spatial + 1 time dimensions · MIT**

A separate six-dimensional Julia project adapted from the [original 3+1 MMIEF repository](https://github.com/DancingAnYoN/late-time_tails_in_einstein_klein_gordon_system). The code evolves a massless, spherically symmetric scalar field around a Schwarzschild--Tangherlini black hole and records waveforms at fixed *areal* radii. Everything below is run **inside the Julia REPL**: no PowerShell commands or shell flags are required.

Two physics modes are available:

- **Fixed Schwarzschild--Tangherlini background:** a test scalar on an exact six-dimensional vacuum black hole (`fixed_background=true`, or `main(["--quick", "--fixed"])`). The generic spherically symmetric linear tail has the reference exponent `p=10`.
- **Partially constrained Einstein--Klein--Gordon prototype:** includes gravitational backreaction (`fixed_background=false`). It has not been convergence-validated. The weak-data `p=6` nonlinear tail of Bizoń--Chmaj--Rostworowski was obtained for dispersal around Minkowski spacetime, **not** automatically for an existing black hole.

**Validation status:** the continuum vacuum checks were performed independently, but Julia and GLMakie were unavailable when this prototype was created. The complete executable suite and the visualizations must be tested locally; no measured Price-law or QNM result is certified.

## 1. Start Julia in VS Code

1. Unzip the project, then use **File → Open Folder** in VS Code to open the folder containing `Project.toml`, `src/`, `scripts/`, and `viz/`.
2. Start the **Julia REPL** with the Julia VS Code extension (`Julia: Start REPL` from the Command Palette, or the extension's Start REPL button). Look for a `julia>` prompt, **not** `PS C:\\>`.
3. Check the working directory with `pwd()`. If necessary, change it to your *actual* project folder:

   ```julia
   cd(raw"C:\Git\tails_five_plus_one") # example: replace with your folder
   pwd()
   isfile("Project.toml")   # should be true
   ```

   The rest of this README assumes the REPL's working directory is the project root. The absolute path may differ on your machine.

### 1.1 Install the core environment (once)

```julia
using Pkg
root = pwd()                             # save the project root
Pkg.activate(root)                      # core physics environment
Pkg.instantiate()                       # install any missing dependencies
Pkg.test()                              # run tests (recommended)
using MMIEFD6
```

The core simulation environment does not load GLMakie. You can return to it at any point with `Pkg.activate(root)`.

### 1.2 Generate a simulation from the REPL

```julia
include("scripts/run_simulation.jl")     # defines main(args)
data = main(["--smoke", "--vacuum"])  # short vacuum check
# Or choose one of the following:
data = main(["--quick", "--fixed"])    # linear test scalar, 65 time units
data = main(["--quick"])              # coupled EKG prototype
data = main(["--standard", "--fixed"]) # longer fixed-background run
# data = main(["--standard"])          # coupled, more expensive
# data = main(["--long", "--fixed"])   # expensive exploratory integration
```

The last executed command writes a *separate* directory under `results/`. Examples: `results/smoke_vacuum`, `results/quick_fixed`, `results/quick`, `results/standard_fixed`. The ZIP contains **no precomputed runs**. Saved data include `diagnostics.csv`, `observers.csv`, `snapshots.csv`, `final_state.csv`, `config.toml`, and `metadata.toml`.

**Alternative for custom experiments:** rather than editing the runner, construct a `Config` and save explicitly:

```julia
p = MMIEFD6.Config(rmax=160.0, dr=0.5, cfl=0.10, tmax=65.0,
                   output_dt=0.25, snapshot_dt=1.0, profile_stride=2,
                   observers=[10.0, 30.0, 50.0],
                   pulse_amplitude=1e-5, fixed_background=true)
custom = MMIEFD6.run_evolution(p)
MMIEFD6.save_run(custom, joinpath(root, "results", "my_fixed_run"))
```

`output_dt` controls waveform/diagnostic sampling; `snapshot_dt` controls saved movie frame spacing; `profile_stride` reduces radial points per frame. These require a **new simulation** to change. The Julia Courant number is `cfl`, equivalent to `C=Δt/Δr`.

## 2. Visualize from the same Julia REPL (GLMakie)

### 2.1 Prepare the visualization environment (first time only)

After a run is saved, switch environments **within the existing Julia REPL**:

```julia
using Pkg
Pkg.activate(joinpath(root, "viz"))
Pkg.develop(path=root)                  # register local MMIEFD6 package
Pkg.instantiate()                      # includes GLMakie
```

The first GLMakie setup may take some time and require internet access. Changing the active environment **does not change** the working directory. If `root` was lost after restarting Julia, use `root = pwd()` while in the project folder.

### 2.2 One-time load of visualization functions

```julia
using MMIEFD6
include(joinpath(root, "viz", "Visualization.jl"))
using .MMIEFD6Visualization
run_dir = joinpath(root, "results", "quick_fixed")
run = MMIEFD6.load_run(run_dir)
```

You can inspect the data directly:

```julia
run.settings.observers                 # default: [10.0, 30.0, 50.0]
length(run.times)                      # observer/diagnostic samples
length(run.frame_times)                # movie snapshots
run.times[end]                         # final coordinate time
run.settings.rh0                       # initial horizon radius, not ADM mass
```

**No simulation rerun is needed** to change or regenerate plots from an existing run. If you restart Julia, repeat §1's `cd(...)`, set `root = pwd()`, then repeat §2.2 in the `viz` environment.

### 2.3 Six-panel scientific dashboard

```julia
make_dashboard(run; outdir=run_dir, show_window=true)
```

Saves `results/quick_fixed/dashboard.png`. The panels show: (1) last `J²φ` radial profile and horizon; (2) horizon/outer mass aspects; (3) signed fixed-`J` waveforms; (4) log–log `|φ|`; (5) local power index with `p=10` reference; (6) Hamiltonian and scalar-definition residuals. `show_window=true` requests a GLMakie window on a local graphical desktop. The file is saved even with `show_window=false`.

A quick run has `tmax=65`, but the dashboard LPI begins near `T=50rh0=100` for the default `rh0=2`. The LPI panel may **correctly be empty** in a quick run; use a longer dataset for later-time analysis.

### 2.4 Full scalar-pulse animation

```julia
make_pulse_movie(run;
    outdir=run_dir,
    fps=24,             # frames per second in the MP4
    every=1,            # 1 = use all saved snapshots; 2 = every second
    radial_max=100.0,   # physical maximum J (same length units as rh0)
    show_window=false)
```

Saves `scalar_pulse_5plus1.mp4` and `scalar_pulse_final.png`. Panels show the evolving normalized `J²φ` radial profile and horizon line, growing horizon mass aspect, and accumulated fixed-`J` observer signals. `fps`, `every`, and `radial_max` control *rendering*, not the PDE. At least two stored snapshots are required. To save a shorter movie:

```julia
make_pulse_movie(run; outdir=run_dir, fps=30, every=3, radial_max=70.0)
```

### 2.5 Fixed-areal-radius observer history

```julia
make_observer_history(run; outdir=run_dir, show_window=true)
```

Saves `observer_history_5plus1.png` showing signed `φ(T,J_obs)` and `log₁₀|φ|` on *linear time* axes. Here `T=t-f(t)`, and the observers remain at fixed physical `J_obs`, with moving computational locations `r_obs(t)=J_obs-f(t)`. This is a semilog plot, whereas the dashboard's tail-amplitude panel is log–log.

### 2.6 One REPL command to render all three

The repository also provides the small REPL helper `scripts/repl_visualizations.jl`:

```julia
include(joinpath(root, "scripts", "repl_visualizations.jl"))
visualize_run(run_dir; show_window=false)  # dashboard, history, MP4
visualize_run(run_dir; movie=false)        # two PNGs only
visualize_run(run_dir; fps=30, every=2, radial_max=120.0)
```

Run this helper **instead of** §2.2's `include("viz/Visualization.jl")` if you are starting a fresh session: it imports the visualization module itself. To avoid redefining the module, include each setup file once per REPL session. The helper returns paths to the saved outputs.

## 3. Analyze waveforms directly in the REPL

The built-in tail and QNM fits are **text/number analyses**, not separate figure generators. The core functions remain available from `MMIEFD6` while the `viz` environment is active.

```julia
T = run.times .- run.shifts
for (J, phi) in zip(run.settings.observers, run.signals)
    tail = MMIEFD6.fit_tail(T, phi; tlo=180.0, thi=300.0)
    println("J=$J: ", tail)
end
```

Choose a dataset that actually spans the requested window, for example `results/standard_fixed`. `fit_tail` returns `nothing` if the window lacks enough usable points. An analogous exploratory QNM fit is:

```julia
for (J, phi) in zip(run.settings.observers, run.signals)
    qnm = MMIEFD6.fit_qnm(T, phi;
        tlo=20.0, thi=60.0, length_scale=run.settings.rh0)
    println("J=$J: ", qnm)
end
```

Use consistent windows and check for oscillations and noise; fitted frequencies are not automatically physical QNMs. **Time-window endpoints here are in raw simulation units**, not multiples of `rh0`. For detailed inspection and a custom GLMakie plot of an individual observer, see [the full visualization guide](docs/VISUALIZATION_GUIDE.md).


## Dedicated late-time tail extraction (Julia REPL)

The project now has a **separate tail-extraction module** with fits over chosen
time windows, fit overlays, local power indices, window-stability scans, and
independent three-resolution waveform comparisons. These are post-processing
tasks: existing data can be analyzed without rerunning the PDE.

After installing/activating the GLMakie `viz` environment from §2, **in a fresh
REPL session** (or before including `scripts/repl_visualizations.jl`):

```julia
include(joinpath(root, "scripts", "repl_tail_analysis.jl"))

report = tail_diagnostics_run("results/standard_fixed";
    tlo=180.0, thi=300.0, window_scan=true, show_window=true)

# Requires three genuinely independent runs at dr=0.50, 0.25, 0.125:
convergence = tail_convergence_runs(
    ["results/h050", "results/h025", "results/h0125"];
    tlo=180.0, thi=300.0, show_window=true)
```

Outputs include `tail_extraction_5plus1.png`, `tail_fit_summary.csv`, `tail_local_power_index.csv`,
`tail_window_stability_5plus1.png`, `tail_window_scan.csv`,
`tail_convergence_5plus1.png`, and convergence summary/curve CSVs.

**Full REPL instructions, resolution-run recipes, scientific meaning and
limitations:** [Late-Time Tail Extraction Guide](docs/LATE_TIME_TAIL_EXTRACTION.md).


## 4. Coordinates, mass, and interpretation

The metric and MMIEF gauge are

\[
 ds^2=-N^2dt^2+a^2(dr+\beta\,dt)^2+J^2d\Omega_4^2,
 \qquad J=r+f(t),\quad N=a(1-\beta).
\]

The evolved grid variables are `(a,B,Phi,Pi,phi)` and the global coordinate shift `f`. The radial curvature `A` is reconstructed from the momentum constraint, while `a` is projected using the radial mass constraint during RK4 stages. The fixed-background mode keeps the geometry equal to Schwarzschild--Tangherlini and sets `f=0`.

- `rh0=2.0` is the initial **horizon radius**, not the physical ADM mass. Time and distance should be normalized by `rh0` (not by a four-dimensional `M`).
- The geometric mass aspect is \(\mu=(J^3/2)(1-a^{-2}+J^2B^2)\), and the physical mass is \(M=4\pi\mu/(3G_6)\).
- Observers `J_obs=[10,30,50]` are fixed *areal radius* sampling points; their numerical grid coordinate changes with the shift.
- Saved snapshot profiles are `J²φ`, **not** the `Jφ` used in the original 3+1 project.
- The fixed six-dimensional reference exponent `p=10` is theoretical; the `p=6` dispersive, nonlinear Minkowski result does not apply automatically to black-hole-plus-scalar initial data.

## 5. Repository contents

```text
Project.toml                  Julia physics package, no Makie dependency
src/                          MMIEF equations, RK4, I/O, diagnostics, fits
src/TailExtraction.jl         Power-law fits, window scans and convergence in core Julia
scripts/run_simulation.jl      Simulation runner; include in the REPL
scripts/repl_visualizations.jl One-session dashboard/history/movie + tail helpers
scripts/repl_tail_analysis.jl REPL-first standalone tail-workflow entry point
scripts/render_*.jl           Optional standalone renderer entry points
scripts/analyze_{tails,qnm}.jl Optional standalone fit entry points
viz/Project.toml              GLMakie dependency (separate environment)
viz/Visualization.jl          make_dashboard / make_pulse_movie / make_observer_history
viz/TailVisualization.jl      Tail fit, window-stability and convergence figures
docs/LATE_TIME_TAIL_EXTRACTION.md  Dedicated Julia REPL tail-analysis instructions
docs/VISUALIZATION_GUIDE.md   Full Julia REPL visualization tutorial
VALIDATION.md                 Validation limitations and tests to perform
ARCHITECTURE.md               Migration from 3+1 to 5+1
results/                      Generated CSV, TOML, PNG, and MP4 files
```

The output schema is `schema_version=2`: `config.toml`, `metadata.toml`, `diagnostics.csv`, `observers.csv`, `snapshots.csv`, and `final_state.csv`. Previous 3+1 result directories are not compatible.

## 6. Scientific sources and validation

- Marsa & Choptuik, *Black Hole--Scalar Field Interactions in Spherical Symmetry*, PRD **54** (1996): https://arxiv.org/abs/gr-qc/9607034.
- Cardoso, Yoshida, Dias & Lemos, *Late-Time Tails of Wave Propagation in Higher Dimensional Spacetimes*, PRD **68** (2003): https://arxiv.org/abs/hep-th/0307122.
- Bizoń, Chmaj & Rostworowski, *Late-Time Tails of a Self-Gravitating Massless Scalar Field, Revisited*, CQG **26** (2009): https://arxiv.org/abs/0812.4333.
- Yoshino & Shibata, *Higher-Dimensional Numerical Relativity: Formulation and Code Tests*, PRD **80** (2009): https://arxiv.org/abs/0907.2760.

**Do not report numerical QNM frequencies or tail exponents as verified results without resolution, boundary-position, gauge, and fit-window convergence studies.** To publish this as a separate repository, see [`GITHUB_SETUP.md`](GITHUB_SETUP.md).
