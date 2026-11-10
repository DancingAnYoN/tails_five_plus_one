# Port map from 3+1 to 5+1

The original modular architecture was kept but its physical coefficients were ported explicitly, not by global search-and-replace.

| Module | Important 5+1 change |
| --- | --- |
| `Types/Grid` | horizon radius `rh0` replaces 3+1 mass length `M0`; explicit `G6` |
| `Geometry` | Tangherlini `q=(rh0/J)^3`, `mu=J^3(...) /2`, `fdot=2pi G6 J^2(Phi+Pi)^2/(3a^2)` |
| `ScalarField` | `J^4` divergence and `-4 fdot Pi/J` |
| `Constraints` | 5-spatial-dimensional ADM coefficients: `(3/2)(a^3-a)/J`, `2pi G6 Phi Pi/a`, `mu'~2pi G6 J^4` |
| `EinsteinEquations` | angular curvature `3(a-a^-1)(1-beta)/J^2` and extrinsic trace `A+4B` |
| `InitialData` | vacuum A/B and derivative in Tangherlini, asymptotic ingoing `Pi=Phi+2(1+aJB)phi/J` |
| `BoundaryConditions` | asymptotic outgoing `J^2 phi` (only approximate at finite radius) |
| `Evolution` | fixed-background mode plus coupled horizon-locked mode, dimension-specific mass budget |
| `Observers/Diagnostics` | fixed areal J; geometric mu and physical 6D mass separated |
| `DataOutput` | schema 2, `J2phi` snapshots, `rh0`, explicit `G6` |
| `Analysis/Viz` | six-dimensional labels; no imported four-dimensional QNM frequencies |

The finite-difference, RK4 and generic Method of Lines infrastructure are the same as in the original source and are tested independently here.

The user project contained historical 3+1 `legacy/` and QNM comparison scripts; they have **not** been brought into the new package as they encode four-dimensional physics. The GitHub repository has not been created online: this ZIP is a standalone, GitHub-ready source tree.
