# Late-Time Tail Extraction in 5+1 Dimensions — Julia REPL Guide

This guide applies to the `MMIEFD6` repository (5 spatial + 1 time dimensions).
It uses **the Julia REPL inside VS Code**, not PowerShell. It adds three
independent tasks beyond the original `make_dashboard` plot:

1. A per-observer **log–log tail-amplitude plot with an actual power-law fit**, a
   local power index (LPI) and a logarithmic fit residual.
2. A **fit-window stability scan**, displaying fitted powers and residuals for
   many explicitly chosen windows.
3. A **three-resolution self-convergence comparison**, including waveform
   differences and tail exponents at all fixed-areal-radius observers.

The corresponding code is in `src/TailExtraction.jl` (no Makie dependency) and
`viz/TailVisualization.jl` (GLMakie plotting), and the easy-to-use REPL functions
are in `scripts/repl_visualizations.jl`.

> **Scientific status:** These are numerical diagnostics, not a proof of Price-law
> convergence. An apparent `p≈10` does not establish the asymptotic decay law.
> A `p≈6` comparison is *not* automatically justified in the black-hole-plus-
> scalar coupled setup; it arose from small-data dispersal in a different
> geometry. Use the convergence and outer-boundary checks described below.

## 1. Set up the Julia REPL (once per session)

Open the repository folder in VS Code and select **Julia: Start REPL**.
At the `julia>` prompt:

```julia
using Pkg
cd(raw"C:\Git\tails_five_plus_one")  # change to YOUR repository path
root = pwd()
@assert isfile(joinpath(root, "Project.toml"))

# Core package and solver; run Pkg.instantiate() on first setup.
Pkg.activate(root)
Pkg.instantiate()
using MMIEFD6
```

A new clone contains no simulation results. To generate saved data, use the
existing REPL simulation runner:

```julia
include(joinpath(root, "scripts", "run_simulation.jl"))
main(["--standard", "--fixed"])  # creates results/standard_fixed
# main(["--standard"])           # coupled 5+1 EKG research prototype
```

A *quick* run ends around `T≈65`. It is useful for smoke tests but cannot yield
an exponent in `[180,300]`; such fits require a sufficiently long run. The
standard run ends near `T≈350` in the fixed-background gauge (`f=0`).

To render, activate the **separate visualization environment in the same REPL**:

```julia
Pkg.activate(joinpath(root, "viz"))
Pkg.develop(path=root)  # first setup only
Pkg.instantiate()      # first setup only; installs GLMakie
using MMIEFD6

include(joinpath(root, "scripts", "repl_tail_analysis.jl"))
```

Including `repl_tail_analysis.jl` defines `tail_diagnostics_run` and
`tail_convergence_runs`, and imports the GLMakie modules. You may instead
include `scripts/repl_visualizations.jl` (which also exposes
`visualize_run`), **but do not include both scripts in the same session**.

## 2. Power-law fit overlays + local power index for all observers

```julia
report = tail_diagnostics_run(
    joinpath(root, "results", "standard_fixed");
    tlo=180.0,
    thi=300.0,
    method=:absolute,
    window_scan=true,
    show_window=true,
)
```

This creates the following in `results/standard_fixed/`:

| Output | Contents |
|---|---|
| `tail_extraction_5plus1.png` | One row per fixed-J observer; log–log `|φ|` with a fitted curve, LPI and log-residual plot |
| `tail_fit_summary.csv` | `J,tlo,thi,power,amplitude,log_rms,r2,n,method,time_scale` |
| `tail_local_power_index.csv` | `J,T,T_over_rh0,local_power_index` for saved finite index samples |
| `tail_window_stability_5plus1.png` | Fit exponent and log-RMS as the start of the fit window varies |
| `tail_window_scan.csv` | Complete scan data; useful for reproducibility and replots |

The fitted model is

\[
 |\varphi(T,J_{\rm obs})|\approx A(J_{\rm obs})
 \left(\frac{T}{r_{h0}}\right)^{-p(J_{\rm obs})},\qquad
 T=t-f(t).
\]

The returned Julia value contains the paths and a `fits` table:

```julia
report.diagnostics.figure
report.diagnostics.summary
for row in report.diagnostics.fits
    println("J=", row.J, "  result=", row.fit)
end
```

The time-window endpoints are **raw `T` units**, not normalized `T/rh0`.
The plot axes *are* normalized using `rh0` for readability.

For the **fixed-background** model, the reference line `p=10` is shown. For
coupled runs no theoretical reference line is selected automatically. You may
supply `reference_power=6.0` manually for an exploratory comparison, but its
physical relevance must be separately justified.

### Direct fit without loading GLMakie

In the core environment:

```julia
using MMIEFD6
run = load_run(joinpath(root, "results", "standard_fixed"))
T = run.times .- run.shifts

for (J, ϕ) in zip(run.settings.observers, run.signals)
    fit = fit_power_law(T, ϕ;
                        tlo=180.0, thi=300.0,
                        time_scale=run.settings.rh0)
    println("J=", J, "  fit=", fit)
end
```

The fit reports `power` (p), `amplitude` (A), `log_rms` (natural-log
residual RMS), `r2`, and the number of valid samples. Those are descriptive
fit statistics, **not** confidence intervals for the continuum exponent.

### Alternative for oscillatory fields with zero crossings

If the field crosses zero repeatedly, fitting `log(abs(ϕ))` can lead to dips.
You can choose a local running-maximum envelope as the fit target:

```julia
report = tail_diagnostics_run(
    "results/standard_fixed";
    tlo=180.0, thi=300.0,
    method=:envelope,
    envelope_halfwidth=4.0,
    peak_halfwidth=4.0,
    lpi_halfwidth=20.0,
)
```

A running maximum biases a clean monotonic power law when its width is
nonzero. Check the result against `method=:absolute` and vary its width.
When using `method=:envelope`, the fit residual panel still compares the raw
amplitude to the fitted smooth model, so zero crossings appear as negative
residual excursions.

Use `amplitude_floor=1e-20` (in your actual scalar-field amplitude units) if
needed to exclude known noise-floor data; choose it from numerical convergence
measurements, not to artificially produce the expected exponent.

## 3. Window stability: is the inferred exponent stable?

The `window_scan=true` option above runs a default scan with the same width
`thi-tlo`, moving the start across all times that leave a full window. For a
more deliberate comparison across **multiple widths**, use the plotting
function directly:

```julia
run = MMIEFD6.load_run("results/standard_fixed")

scan = MMIEFD6TailVisualization.make_tail_window_stability(
    run;
    outdir="results/standard_fixed",
    tlo=150.0,
    thi=260.0,
    starts=collect(150.0:10.0:230.0),
    widths=[80.0,110.0],
    show_window=true,
)
```

Each window is `[start,start+width]`. Missing or under-sampled windows are
omitted from the scan. No window is automatically chosen to match a predicted
exponent. Compare the fitted `p` against `T0/rh0` and the fit residual: a
stable local slope together with an improving residual is encouraging, but
can still reflect a preasymptotic regime.

## 4. Three-resolution self-convergence (must run the solver three times)

You **cannot** test convergence by plotting one run at three sampling rates.
Generate genuinely independent runs with the same physical data, boundary
radius, and measurement radii, and mesh spacings in the ratio 4:2:1.

Here is an example using a **fixed Schwarzschild–Tangherlini background**.
It may be computationally expensive at the finest resolution, so run the
`standard_fixed` simulation first to confirm your machine can handle it.

Switch back to the core Julia environment (same REPL is fine):

```julia
Pkg.activate(root)
using MMIEFD6

function save_resolution_run(h)
    p = Config(
        rh0=2.0, G6=1.0,
        rmax=700.0, dr=h, cfl=0.10,
        tmax=350.0, output_dt=0.5,
        snapshot_dt=10.0, profile_stride=8,
        pulse_r0=20.0, pulse_sigma=2.0,
        pulse_amplitude=1e-5,
        cutoff_start=6.0, cutoff_end=9.0,
        ko_strength=0.03,
        observers=[10.0,30.0,50.0],
        fixed_background=true,
    )
    run = MMIEFD6.run_evolution(p)
    name = h == 0.5 ? "h050" : h == 0.25 ? "h025" : "h0125"
    return MMIEFD6.save_run(run,joinpath(root,"results",name))
end

coarse_dir = save_resolution_run(0.50)
medium_dir = save_resolution_run(0.25)
fine_dir   = save_resolution_run(0.125)
```

The fine run is expected to be costly and all runs need enough samples in the
chosen fit interval. Start with cheaper shorter tests for software debugging;
they will not by themselves validate late-time tails.

After these three runs finish, switch back to GLMakie:

```julia
Pkg.activate(joinpath(root,"viz"))

convergence = tail_convergence_runs(
    [coarse_dir, medium_dir, fine_dir];
    tlo=180.0,
    thi=300.0,
    outdir=joinpath(root,"results","tail_convergence"),
    show_window=true,
)
```

This saves:

| Output | Contents |
|---|---|
| `tail_convergence_5plus1.png` | Per-observer waveform differences and fitted power against `dr/rh0` |
| `tail_convergence_summary.csv` | Resolution spacings, shared window, mismatch norms, apparent order, powers |
| `tail_convergence_curves.csv` | Resampled observer time, three fields, and their pairwise differences |

The comparison aligns the three waveforms at the **same fixed areal radius**
`J_obs`, using observer time `T=t-f(t)`. The medium/fine data are interpolated
onto the coarse-run sample times, and there is **no extrapolation**.

For matching refinement factor `r`, the reported apparent order is

\[
 q_{\rm app}=\frac{\log\left(\Vert\varphi_h-\varphi_{h/r}\Vert_2/
                   \Vert\varphi_{h/r}-\varphi_{h/r^2}\Vert_2\right)}{\log r}.
\]

The norm is a discrete RMS over the common coarse observer times. A reliable
order estimate also needs time-step refinement, verified outer-boundary
independence and a regime where the leading truncation error dominates. The
reported order may be undefined if either difference norm vanishes.

The code **rejects** comparisons with different `rh0`, `G6`, `rmax`, initial
pulse parameters, fixed/coupled mode, damping strength or observer radii.
It also requires approximately geometric mesh spacing ratios.

## 5. Interpreting tail curves: analysis checklist

1. Look at signed `ϕ(T,J_obs)` first; identify QNM ringing, later-time decay,
   sign changes and a possible floating-point/noise floor.
2. In `tail_extraction_5plus1.png`, use the log–log fit as a *candidate* slope
   and inspect the fit residual for systematic curvature.
3. Check the LPI. Its estimated `p(T)` should plausibly approach a plateau;
   a noisy spike or a single crossing of `p=10` is not enough.
4. Scan several fit starts and widths; inspect `tail_window_scan.csv`.
5. Compare at least three grid resolutions and change the outer boundary
   while keeping the observers and physics fixed.
6. Only then report a converged, physically interpreted exponent.

In the *fixed* six-dimensional Schwarzschild spherical test-field model the
reference is `p=10` (Cardoso et al., 2003). The nonlinear dispersive `p=6`
result of Bizoń et al. (2009) does **not** transfer automatically to an
initially existing black-hole spacetime.

## 6. Troubleshooting

| Symptom | Explanation / action |
|---|---|
| `UndefVarError: tail_diagnostics_run` | Include `scripts/repl_tail_analysis.jl` in the current Julia session |
| `Package GLMakie not found` | Activate `viz` and run `Pkg.develop(path=root); Pkg.instantiate()` |
| `No saved run` | Generate results using `main(["--standard","--fixed"])` or `save_run` |
| `fit=nothing` or no fitted curve | Window outside run length, too few valid samples, or all amplitudes below cutoff |
| Empty LPI panel | Run too short for `tlo` / LPI half-width; check `maximum(T)` |
| No scan points | Start plus width is beyond available `T` or there are too few samples |
| Resolution comparison rejects configuration | Read `config.toml` of each run; match all physics/boundary parameters |
| Undefined / wildly varying `q_app` | Differences near zero, pre-asymptotic regime, interpolation contamination, or boundary effects |
| Plot appears but no MP4 | Tail-analysis functions export PNG+CSV only; use `visualize_run(...; movie=true)` for the scalar pulse movie |
| Blank GLMakie window | Set `show_window=false` to save PNG without opening a desktop window; verify graphics driver separately |

## 7. References

- V. Cardoso, S. Yoshida, Ó. J. C. Dias, J. P. S. Lemos,
  *Late-Time Tails of Wave Propagation in Higher Dimensional Spacetimes*,
  [arXiv:hep-th/0307122](https://arxiv.org/abs/hep-th/0307122).
- P. Bizoń, T. Chmaj, A. Rostworowski,
  *Late-Time Tails of a Self-Gravitating Massless Scalar Field, Revisited*,
  [arXiv:0812.4333](https://arxiv.org/abs/0812.4333).
- R. L. Marsa and M. W. Choptuik,
  *Black Hole–Scalar Field Interactions in Spherical Symmetry*,
  [arXiv:gr-qc/9607034](https://arxiv.org/abs/gr-qc/9607034).

**Implementation status:** Julia and GLMakie were unavailable in the artifact
creation environment. The additional unit tests contain exact synthetic
power-law and second-order self-convergence examples, but must still be run
using `Pkg.test()` locally. All produced physical tail exponents require
independent convergence and extraction-window validation.
