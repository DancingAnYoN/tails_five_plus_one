# Julia REPL visualization guide — six-dimensional MMIEF tails

This guide is intentionally **REPL-first**. All runnable commands below are Julia expressions entered at the `julia>` prompt in the VS Code Julia extension; no PowerShell, Bash, or CLI `--flags` are needed. The physics solver and the GLMakie renderer use separate Julia project environments but can be used sequentially from the **same long-lived Julia REPL**.

## 1. What the project actually visualizes

| Output | REPL function | Main scientific use |
|---|---|---|
| `dashboard.png` | `make_dashboard(data; ...)` | Six-panel overview: scalar profile, masses, observer waveforms, log tails, LPI, constraints |
| `observer_history_5plus1.png` | `make_observer_history(data; ...)` | Signed and logarithmic amplitude histories at fixed physical areal radii |
| `scalar_pulse_5plus1.mp4` | `make_pulse_movie(data; ...)` | Time-resolved radial propagation, horizon, mass and accumulated signals |
| `scalar_pulse_final.png` | Saved by `make_pulse_movie` | Static image of final animation frame |
| Printed numerical fit | `MMIEFD6.fit_tail(T, phi; ...)` | Exploratory late-time exponent and logarithmic fit error |
| Printed numerical fit | `MMIEFD6.fit_qnm(T, phi; ...)` | Exploratory damped-mode frequency and damping estimate |

**Important distinction:** the first three are implemented GLMakie graphics; the last two are numerical fitting functions that do **not** automatically produce a fitted-line figure. The examples below also show you how to build your own plot interactively.

## 2. Starting from an empty Julia REPL

### 2.1 Open the right folder in VS Code

Open the extracted `late_time_tails_5plus1_mmief` repository in VS Code. Open the Command Palette and choose `Julia: Start REPL` (or the Julia extension's REPL button). Julia expressions belong at the `julia>` prompt, not at `PS>`.

If your terminal opens somewhere else, change directory from Julia:

```julia
cd(raw"C:\Git\tails_five_plus_one") # example; replace with the extracted project path
pwd()
isfile("Project.toml")          # must be true
isfile("viz/Visualization.jl") # must be true
root = pwd()                   # save path for the whole session
```

The raw-string syntax is useful for Windows paths because you do not need to double each backslash.

### 2.2 Activate and install the simulation environment

The physics package at the project root uses Julia standard libraries and does not require GLMakie:

```julia
using Pkg
Pkg.activate(root)
Pkg.instantiate()
Pkg.test()               # recommended first check, not a validation proof
using MMIEFD6
```

Run a smoke vacuum evolution and a short fixed-background scalar test:

```julia
include("scripts/run_simulation.jl")
main(["--smoke", "--vacuum"]) # very short, geometry only
main(["--quick", "--fixed"])   # linear scalar pulse, tmax=65
```

The runner's `main` function accepts `Vector{String}`; in the REPL **do not** type `--quick` bare. For a longer linear experiment, run `main(["--standard", "--fixed"])`. For self-gravity, run `main(["--standard"])`. To try the longest predefined configuration, use `main(["--long", "--fixed"])` only once the smaller cases work. The combined numerical and graphical workflow is **not** a validated Price-tail measurement.

The resulting saved directories are:

```text
results/smoke_vacuum/
results/quick_fixed/
results/standard_fixed/     (only after a standard fixed run)
results/standard/           (only after a coupled standard run)
```

You can generate custom sampling resolution without editing scripts:

```julia
p = MMIEFD6.Config(
    rh0=2.0, rmax=160.0, dr=0.5, cfl=0.10, tmax=65.0,
    pulse_r0=20.0, pulse_sigma=2.0, pulse_amplitude=1e-5,
    output_dt=0.25, snapshot_dt=0.5, profile_stride=2,
    observers=[10.0,30.0,50.0], fixed_background=true)
custom = MMIEFD6.run_evolution(p)
MMIEFD6.save_run(custom, joinpath(root,"results","custom_fixed"))
```

Increase `tmax` if you need later-time data; extend `rmax` if you need to delay outer-boundary contamination. Reduce `dr` and `cfl` only after checking resource usage. `snapshot_dt` controls *how many frames are stored*; `output_dt` controls waveform/diagnostic time sampling. `profile_stride=2` means keep every second grid point in stored movie profiles. None of those settings can be repaired retroactively by a renderer.

### 2.3 Activate the visualization environment

Once core results exist, **stay in the same REPL**:

```julia
Pkg.activate(joinpath(root,"viz"))
Pkg.develop(path=root)         # link to this exact local MMIEFD6 code
Pkg.instantiate()             # install GLMakie and other viz requirements
using MMIEFD6
```

The `Pkg.develop` step matters: the GLMakie environment contains a dependency on the local MMIEFD6 package. If `Package MMIEFD6 not found` occurs, this step is usually missing. `Pkg.activate` switches dependencies; it does not execute simulation code and does not move the working directory.

The visualization functions live in `viz/Visualization.jl` as a **separate module**, not inside `MMIEFD6`:

```julia
include(joinpath(root,"viz","Visualization.jl"))
using .MMIEFD6Visualization
```

Do that **once per Julia session**. Re-including a module repeatedly can produce a `replacing module` warning. If you already included `scripts/repl_visualizations.jl` instead, skip this manual module include.

### 2.4 Load one recorded run into memory

```julia
run_dir = joinpath(root,"results","quick_fixed")
isfile(joinpath(run_dir,"metadata.toml"))   # should be true
run = MMIEFD6.load_run(run_dir)
```

An existing run is read from CSV/TOML without reintegrating the field equations. Inspect it:

```julia
run.settings.rh0              # default 2.0 = initial horizon radius
run.settings.observers        # fixed physical J, typically [10,30,50]
run.times[1], run.times[end]  # start/end coordinate time
run.shifts[end]               # f(t) at the last saved diagnostic sample
length(run.times)             # diagnostic/observer time samples
length(run.frame_times)       # saved pulse movie frames
size(run.final_state)         # rows: radial points; columns: fields
run.frame_profiles[end][1:10] # sample of the saved J²φ profile
```

Data from the original four-dimensional repository are **not** compatible with this output schema.

## 3. Task A — six-panel diagnostic dashboard

```julia
make_dashboard(run; outdir=run_dir, show_window=true)
```

This creates `results/quick_fixed/dashboard.png`. `show_window=true` also calls `display(fig)` to display the graphics window when GLMakie is working on the desktop; set it to `false` if you only want to save the PNG.

The panels and their meaning:

| Panel | Plotted quantity | Interpretation / caveat |
|---|---|---|
| 1. Last scalar profile | `J²φ` versus `J/rh0` and horizon line | Where the pulse ended; `J²φ` is not the bare scalar and not a guaranteed conserved amplitude |
| 2. Mass aspects | `μ_h/μ0`, `μ_outer/μ0` versus `t/rh0` | Horizon growth and mass budget; the fixed test-field mode deliberately lacks backreaction |
| 3. Observer histories | signed `φ(T,J_obs)` versus `T/rh0` | Pulse arrival, zero crossings, potential damped ringing |
| 4. Tail amplitude | `log10|φ|` versus `log10(T/rh0)` | Candidate asymptotic power law; zero crossings and noise make dips/spikes |
| 5. Local power index | `p(T)=-d ln(env)/d ln T` versus `T/rh0` | Stabilization toward a constant would support a tail hypothesis; the horizontal `p=10` is theoretical only |
| 6. Constraint diagnostics | `log10` Hamiltonian and scalar-definition L² residuals versus `t/rh0` | Detect numerical inconsistency; fixed-background Hamiltonian residual includes scalar source effects |

The dashboard uses `MMIEFD6.local_power_index(T, φ; peak_halfwidth=12rh0, fit_halfwidth=24rh0, tstart=50rh0)`. With default `rh0=2`, it begins around `T=100`; a quick run ends at `t=65`, so **its LPI panel is expected to be empty**. That is not a GLMakie error. The default line at `p=10` is a reference for the linear six-dimensional Schwarzschild spherical mode, **not** an automatically recovered numerical result.

To regenerate a different dataset, simply load it and call the same function:

```julia
run_dir = joinpath(root,"results","standard_fixed")
run = MMIEFD6.load_run(run_dir)
make_dashboard(run; outdir=run_dir, show_window=true)
```

## 4. Task B — saved scalar pulse MP4

```julia
make_pulse_movie(run;
    outdir=run_dir,
    fps=24,
    every=1,
    radial_max=100.0,
    show_window=false)
```

This creates **two files** in `run_dir`: `scalar_pulse_5plus1.mp4` and `scalar_pulse_final.png`. The function uses the saved snapshots (`run.frame_times`, `run.frame_shift`, `run.frame_profiles`) and does not call the PDE solver.

The movie displays three synchronized views:

1. **Large radial panel**: scaled `J²φ/(rh0² ε)` against `J/rh0` for the selected snapshot, with the moving horizon at `J_h/rh0` and a zero reference line. The normalization uses the configured pulse amplitude `ε`.
2. **Horizon mass**: accumulated `μ_h/μ0` against `t/rh0`, updated as the movie progresses.
3. **Observer signals**: accumulated `φ(T,J_obs)` for all configured observers, so the waveform grows as the pulse propagates.

Arguments (these are **function keyword arguments**, not command-line flags):

| Keyword | Meaning | Example |
|---|---|---|
| `outdir` | Where output files are saved | `outdir=run_dir` |
| `fps` | MP4 playback rate; not integration time step | `fps=30` |
| `every` | Select every N-th **saved** snapshot; final frame always included | `every=3` |
| `radial_max` | Requested maximum physical areal radius `J` in run length units | `radial_max=120.0` |
| `show_window` | Display the final-frame GLMakie figure (requires GUI) | `show_window=true` |

Examples, all on the existing dataset:

```julia
# Higher playback frame rate and reduced frame count:
make_pulse_movie(run; outdir=run_dir, fps=30, every=2, radial_max=80.0)

# Keep every stored frame, and show final frame in a window:
make_pulse_movie(run; outdir=run_dir, fps=18, every=1,
                 radial_max=140.0, show_window=true)
```

`fps` affects how long the MP4 lasts; `every` affects the number of saved frames rendered. Changing them **does not make the physical simulation evolve faster or slower**. If you want more *physical* frames, reduce the simulation setting `snapshot_dt` **and rerun**. If you want more radial detail per frame, reduce `profile_stride` **and rerun**. The MP4 requires at least two saved snapshots.

**Scientific caveat:** the movie is useful for qualitative inspection of infall, reflection, propagation, and changing horizon geometry. Its playback is not a quantitative method to measure QNMs or tails.

## 5. Task C — full observer waveform history

```julia
make_observer_history(run; outdir=run_dir, show_window=true)
```

Saves `observer_history_5plus1.png`, containing:

- Left: signed field `φ(T,J_obs)` versus *linear* `T/rh0` (oscillations and sign changes visible).
- Right: `log10|φ(T,J_obs)|` versus *linear* `T/rh0` (semilog attenuation history).

This is distinct from panel 4 of the dashboard, which uses **log–log** axes and is appropriate for exploratory power-law identification. The history plot reads `run.signals` and does not fit an exponent.

The default observer radii are fixed **areal radii** `J_obs = [10,30,50]`. Because the numerical radius is `r` and the MMIEF areal radius is `J=r+f(t)`, an observer stays at `r_obs(t)=J_obs-f(t)`. The plotted exported time variable is `T=t-f(t)`, not the observer's proper time.

### 5.1 Inspect a single observer interactively, without rerunning

Once `run` has been loaded, you can construct custom GLMakie figures in the REPL:

```julia
using GLMakie
T = run.times .- run.shifts
rh0 = run.settings.rh0
k = 2                                   # second configured observer
J = run.settings.observers[k]
phi = run.signals[k]

fig = Figure(size=(1150,700))
ax1 = Axis(fig[1,1], title="Waveform at J=$J", xlabel="T/rh0", ylabel="φ")
ax2 = Axis(fig[2,1], title="Amplitude at J=$J", xlabel="T/rh0", ylabel="log₁₀|φ|")
lines!(ax1, T ./ rh0, phi)
lines!(ax2, T ./ rh0, log10.(max.(abs.(phi), 1e-300)))
save(joinpath(run_dir,"single_observer_J$(Int(round(J))).png"), fig)
display(fig)
```

This code can select only the observers **already saved**. To add an observer at an entirely new areal radius, specify a new `observers` array in the simulation configuration and rerun: the output does not contain full time-resolved radial fields at arbitrary temporal cadence.

### 5.2 Custom log–log tail view with theoretical slope

```julia
using GLMakie
T = run.times .- run.shifts
rh0 = run.settings.rh0
phi = run.signals[2]
valid = findall(i -> T[i] > 0 && abs(phi[i]) > 1e-290, eachindex(T))
fig = Figure(size=(1000,600))
ax = Axis(fig[1,1], xlabel="log₁₀(T/rh0)", ylabel="log₁₀|φ|",
          title="Exploratory tail amplitude")
lines!(ax, log10.(T[valid] ./ rh0), log10.(abs.(phi[valid])))
save(joinpath(run_dir,"custom_tail_view.png"), fig)
display(fig)
```

This is an **exploratory view only**: a finite window may contain ringdown, zeros, boundary reflections, or noise. Do not infer an asymptotic exponent solely from the plot.

## 6. Run all figures from one persistent REPL

### 6.1 Direct functions (no new helper)

Once the visualization module and dataset have been loaded:

```julia
make_dashboard(run; outdir=run_dir, show_window=false)
make_observer_history(run; outdir=run_dir, show_window=false)
make_pulse_movie(run; outdir=run_dir, fps=24, every=1,
                 radial_max=100.0, show_window=false)
```

The outputs are saved to the existing dataset folder. Repeating calls overwrites the files of the same names. You can switch datasets without switching environments:

```julia
run_dir = joinpath(root,"results","standard_fixed")
run = MMIEFD6.load_run(run_dir)
make_dashboard(run; outdir=run_dir)
```

### 6.2 Dedicated one-command REPL helper

Alternatively, from a **fresh visualization REPL session**, after completing §2.3 through `using MMIEFD6` but **before** manually including the visualization module:

```julia
include(joinpath(root,"scripts","repl_visualizations.jl"))
visualize_run("results/quick_fixed")                  # all outputs
visualize_run("results/quick_fixed"; movie=false)      # PNGs only
visualize_run("results/quick_fixed"; fps=30, every=2)  # configurable MP4
visualize_run(joinpath(root,"results","standard_fixed"); radial_max=120)
```

The helper loads the requested saved run and generates files. It accepts `dashboard`, `history`, `movie`, `fps`, `every`, `radial_max`, and `show_window` keyword arguments. It returns a named tuple of output paths. It does not rerun simulations and can be called repeatedly without re-including the source module.

## 7. Numerical tail fitting from the REPL

The original `scripts/analyze_tails.jl` is a terminal entry point, but the underlying Julia function can be called directly:

```julia
run_dir = joinpath(root,"results","standard_fixed")
run = MMIEFD6.load_run(run_dir)
T = run.times .- run.shifts
for (J, phi) in zip(run.settings.observers, run.signals)
    result = MMIEFD6.fit_tail(T, phi; tlo=180.0, thi=300.0)
    println("J=$J | ", result)
end
```

For a sufficiently sampled window, `result` includes `power`, `log_rms`, `n`, `tlo`, and `thi`. If there are fewer than 10 usable samples, it returns `nothing`. The time window uses **raw simulation units**, not `T/rh0`. If you want to specify window endpoints in units of `rh0`, explicitly multiply each endpoint by `run.settings.rh0`.

The expected fixed Schwarzschild background reference for the six-dimensional spherical test field is `p=10`. A fit producing `p≈10` without resolution/outer-boundary/window convergence is **not** a verification of Price's law. The `p=6` nonlinear small-data result concerns a different dispersive physical setting.

To compute a local power index without manually differentiating noisy logs:

```julia
phi = run.signals[1]
rh0 = run.settings.rh0
p = MMIEFD6.local_power_index(T, phi;
    peak_halfwidth=12rh0, fit_halfwidth=24rh0, tstart=50rh0)
finite = findall(isfinite, p)
using GLMakie
fig = Figure(size=(950,500))
ax = Axis(fig[1,1], xlabel="T/rh0", ylabel="local power index")
lines!(ax, T[finite] ./ rh0, p[finite])
hlines!(ax, [10.0], linestyle=:dash, color=:black)
save(joinpath(run_dir,"custom_LPI.png"), fig)
display(fig)
```

The code above computes only from the already recorded observer data; an empty `finite` set for a short run is normal.

## 8. Exploratory QNM fitting from the REPL

```julia
run_dir = joinpath(root,"results","quick_fixed")
run = MMIEFD6.load_run(run_dir)
T = run.times .- run.shifts
for (J,phi) in zip(run.settings.observers,run.signals)
    result = MMIEFD6.fit_qnm(T, phi;
        tlo=20.0, thi=60.0, length_scale=run.settings.rh0)
    println("J=$J | ", result)
end
```

The result provides an exploratory damped-oscillation frequency, decay rate and error estimate, or `nothing` when no stable single-mode fit is found. It is **not** compared automatically with a six-dimensional theoretical QNM frequency. Check for sufficiently many oscillations, isolation of the ringdown phase, stability across radii and windows, and convergence. There is no standalone plotted-fit script included; the dashboard and single-observer REPL plot supply the waveform context.

## 9. Recovering from common Julia REPL problems

| Symptom | Cause / action in the REPL |
|---|---|
| `ERROR: UndefVarError: root not defined` | At project root run `root = pwd()`; check `isfile(joinpath(root,"Project.toml"))`. |
| `Package MMIEFD6 not found` | From the project root: `Pkg.activate(joinpath(root,"viz")); Pkg.develop(path=root); Pkg.instantiate()`; then `using MMIEFD6`. |
| `Package GLMakie not found` | You are probably still in the core environment. Run `Pkg.activate(joinpath(root,"viz")); Pkg.instantiate()` before `using GLMakie` or including `viz/Visualization.jl`. |
| `UndefVarError: make_dashboard not defined` | Include `viz/Visualization.jl` and run `using .MMIEFD6Visualization`, or use `MMIEFD6Visualization.make_dashboard(...)`. |
| `WARNING: replacing module MMIEFD6Visualization` | Avoid repeatedly including the same module. The functions are already loaded; for a clean source reload, restart the REPL. |
| `metadata.toml not found` | Your `run_dir` has no saved simulation; generate data using `main(["--quick","--fixed"])` in the **core** environment first. |
| `Unsupported ... schema version` | You loaded 3+1 or older-format data; use a run created by the 5+1 solver. |
| `Need at least two saved snapshots` | The chosen run lacks movie frames; use a longer run or smaller `snapshot_dt` in a fresh simulation. |
| GLMakie/OpenGL/GLFW errors | Use local VS Code desktop session with functioning graphics drivers. For remote/headless rendering, this GLMakie renderer is not guaranteed to work. |
| MP4 encoder error | First test the dashboard PNG to isolate OpenGL setup from video encoding; inspect the video encoder error. |
| LPI empty | `tmax` is too short or the signal is under the numerical threshold; the quick run usually ends before the default LPI start. |
| `MethodError` after changing code | Check types and active environment; if you edited a loaded module, restart the REPL for a clean reload. |

### 9.1 Switching back from visualizations to evolution

From your persistent REPL:

```julia
Pkg.activate(root)       # back to core environment
# Simulation functions remain importable through MMIEFD6
p = MMIEFD6.Config(tmax=2.0, rmax=70.0, dr=0.5,
                   output_dt=0.25, snapshot_dt=0.25,
                   fixed_background=true, pulse_amplitude=0.0)
check = MMIEFD6.run_evolution(p)
MMIEFD6.save_run(check, joinpath(root,"results","repl_check"))

Pkg.activate(joinpath(root,"viz")) # back to renderer; no REPL restart needed
run = MMIEFD6.load_run(joinpath(root,"results","repl_check"))
make_dashboard(run; outdir=joinpath(root,"results","repl_check"))
```

*If that smoke dataset has fewer than two frames, skip the movie; dashboard/history still work. The REPL keeps loaded Julia modules in memory while environments are switched; switching environments affects future package resolution, not the already-loaded module definitions.*

## 10. Data file meaning, plotting cadence, and limitations

| Saved file | Main contents | Used by |
|---|---|---|
| `config.toml` | `rh0`, `G6`, simulation settings, observer radii, strides | Every renderer |
| `metadata.toml` | `schema_version=2`, coordinate conventions, data schema | `load_run` |
| `diagnostics.csv` | horizon/outer masses, constraint norms, time and shift | Dashboard and movie |
| `observers.csv` | `(t,T_asymptotic,J_obs,phi)` fixed-radius time series | Dashboard/history/movie/fits |
| `snapshots.csv` | `(t,f,r,J,J2phi)` for saved frames | Dashboard last profile and MP4 |
| `final_state.csv` | Final `(r,a,B,Phi,Pi,phi)` grid state | RunData reconstruction |

The six-dimensional energy and mass scales differ from four dimensions. The geometric mass aspect is

\[
\mu(t,J)=\frac{J^3}{2}\left(1-a^{-2}+J^2B^2\right),
\qquad M=\frac{4\pi}{3G_6}\mu,
\]

where `rh0` is a length scale (the initial horizon radius). Dashboard and movie axes are typically normalized by `rh0`, **not by a four-dimensional Schwarzschild mass**. `J_obs` is areal radius, `r` is grid radius, and `T=t-f(t)` is the exported observation time, not observer proper time.

**Final scientific warning:** these plots are diagnostics from a research prototype. Numerical stability, constraint control, QNM identification, and late-time power-law tails require systematic convergence studies, boundary-reflection checks, and credible extraction windows. For equations and conventions see `docs/EQUATIONS.md`; for remaining limitations see `VALIDATION.md`.


---

## 11. New: dedicated late-time tail extraction and convergence figures

A separate numerical/visual module is provided for fitting and validating
tails at fixed areal radii. In the Julia REPL, after activating the `viz`
environment and including `scripts/repl_tail_analysis.jl`:

```julia
tail_diagnostics_run("results/standard_fixed";
    tlo=180.0,thi=300.0,window_scan=true,show_window=true)
```

This makes per-observer fitted power-law overlays, local power indices,
residuals, a moving-window scan, and reproducible CSV outputs. Once you have
three independently generated resolution runs:

```julia
tail_convergence_runs(["results/h050","results/h025","results/h0125"];
    tlo=180.0,thi=300.0,show_window=true)
```

For the complete REPL tutorial, source-code APIs, creating the three runs,
interpreting the results, and limitations, read
[`LATE_TIME_TAIL_EXTRACTION.md`](LATE_TIME_TAIL_EXTRACTION.md).
