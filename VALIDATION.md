# Validation status and outstanding scientific checks

**What has been checked in the creation environment**

1. Retrieved the user's previously saved MIT-licensed **source ZIP** corresponding to the linked GitHub repository and ported the modular numerical architecture. Network access in the build container was disabled, so it was not freshly cloned; see `ARCHITECTURE.md` for the modifications.
2. Independent SymPy substitution confirms the 5+1 Tangherlini metric satisfies the **continuum** Hamiltonian and momentum constraints exactly, and \(\mu=\mathrm{const}=r_h^3/2\). The horizon-locking coefficient was also checked by symbolic cancellation at \(aJB=1\).
3. All source paths, package UUID dependencies, test entrypoints, and output schema were statically inspected. **Julia was not installed**, so neither `Pkg.test()` nor the numerical solver nor GLMakie was run. Therefore *no scientific tail or QNM result is reported*.

**Required local checks, in order**

1. Install Julia 1.12.1, run `julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'`. Fix any errors found in your environment before experiments.
2. Run `--smoke --vacuum` and confirm `maxabs(RHS)` close to roundoff, horizon error and mass-aspect drift negligible. Inspect constraint norms, and repeat with `dr=0.5`, `0.25`, `0.125` on the same domain.
3. Test `--quick --fixed` with weak scalar data. Confirm `a`, `B`, and `f=0` stay constant and scalar propagation agrees with an independently derived test-field equation on a Tangherlini background. Repeat at multiple resolutions.
4. Test nonlinear initial data at several pulse amplitudes; monitor `hamiltonian_L2`, `scalar_definition_L2`, `horizon_preprojection_error` and the initial-to-final *mass aspect* budget.
5. Move the outer boundary (`rmax`) substantially farther away while keeping physical observation radius fixed; compare observer signals before the first reflected contamination can return. The outer Sommerfeld condition is approximate in six dimensions.
6. Compare the waveforms at fixed **areal radius**, `J=10,30,50` in the chosen length units. For late-time fits compare `T=t-f`, windows and double/quad precision. A 6D `T^-10` field may disappear below floating-point numerical error.
7. Compare fixed-background monopole tail against Cardoso et al. (hep-th/0307122). Do **not** interpret the `T^-6` nonlinear dispersive Minkowski result of Bizoń et al. (0812.4333) as a guaranteed exponent for black-hole accretion.
8. Fit QNMs using an independent 6D scalar QNM theory/ODE calculation; **do not use 4D frequencies 0.110455/0.104896** or a single short-window fit as evidence of agreement.
9. Assess the effect of the numerical well-balancing subtraction and of Hamiltonian projection on the convergence order, especially near excision and outer boundary.

**Boundary/algorithm limitations**: no SBP/SAT boundary energy estimate, no formal hyperbolicity proof for partially constrained discretization, no adaptive step or constraint-preserving boundary implementation, no exact radiative boundary on a finite domain, no generalized harmonic/BSSN hyperbolicity. The `--standard` and `--long` modes are starting configurations only. The result is a *research implementation proposal*, not a verified Price-law result.


## Tail extraction and convergence module

`src/TailExtraction.jl` and `viz/TailVisualization.jl` add windowed power-law
fits, LPI plots, fit-window scans and three-run waveform comparison. The Julia
test suite includes deterministic synthetic `T^-10` and known second-order
waveform-difference checks. These tests have **not** been executed in the
artifact creation environment (Julia and GLMakie unavailable). They do not
validate the PDE solver, the asymptotic regime, or a measured physical tail.

Before reporting any exponent: verify constraints and mass conservation,
show convergence of the waveform and exponent with at least three spatial/time
resolutions, repeat with a more distant outer boundary, confirm a stable
late-time extraction window, and distinguish absolute-value fits from an
envelope fit near waveform zero crossings.
