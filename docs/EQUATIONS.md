# 5+1 MMIEF Einstein--Klein--Gordon: continuum derivation

**Conventions:** \(d+1=5+1\), spatial \(d=5\), sphere \(S^4\) (\(k=d-1=4\)); metric signature \((-+++++)\), MTW curvature sign \([\nabla_\mu,\nabla_\nu]V^\rho=R^\rho{}_{\sigma\mu\nu}V^\sigma\), \(R_{\mu\nu}=R^\rho{}_{\mu\rho\nu}\), and extrinsic curvature \(K_{ij}=-\frac12\mathcal L_n\gamma_{ij}\). We use \(c=1\), retain \(G_6\), and set \(\Lambda=0\).

## Spherical ADM reduction

\[
 ds^2=-N^2dt^2+a^2(dr+\beta dt)^2+J^2d\Omega_4^2,
 \qquad n^\mu=N^{-1}(1,-\beta,0,\ldots),
\]

with \(A=K^r{}_r\), \(B=K^\theta{}_\theta\). Hence

\[
 K=A+4B,\qquad K_{ij}K^{ij}=A^2+4B^2.
\]

For the slice metric \(a^2dr^2+J^2d\Omega_4^2\), impose \(J=r+f(t)\), so \(J'=1\), \(J''=0\). Its independent Ricci components are

\[
 {}^{(5)}R^r{}_r=\frac{4a'}{a^3J},\qquad
 {}^{(5)}R^\theta{}_\theta=\frac{a'}{a^3J}+\frac{3}{J^2}(1-a^{-2}),
\]

and therefore

\[
 {}^{(5)}R=\frac{8a'}{a^3J}+\frac{12}{J^2}(1-a^{-2}).
\]

These follow from the warped-product formulas with a spherical fiber of dimension \(4\) and are independently checked against vacuum Tangherlini.

## Gauge and angular metric evolution

MMIEF imposes \(J=r+f(t)\) and \(g(\partial_t-\partial_r,\partial_t-\partial_r)=0\), giving

\[
 N=a(1-\beta).
\]

The angular metric equation \(\partial_t(J^2)=-2N J^2B+2\beta JJ'\) gives

\[
 \dot f=\beta-NJB.
\]

Write \(Q=aJB\); solving yields

\[
 \boxed{\beta=\frac{\dot f+Q}{1+Q},\qquad N=\frac{a(1-\dot f)}{1+Q}.}
\]

The radial coordinate light speeds are \(c_-=-1\) and \(c_+=1-2\beta\). Gauge regularity requires \(1+Q>0\) and \(\dot f<1\).

## Scalar field

Define

\[
 \Phi=\partial_r\varphi,\qquad
 \Pi=\frac aN(\partial_t\varphi-\beta\Phi).
\]

Then \(n^\mu\partial_\mu\varphi=\Pi/a\),

\[
 \rho=\frac{\Phi^2+\Pi^2}{2a^2},\qquad j_r=-\frac{\Phi\Pi}{a}.
\]

Using \(\sqrt{-g}=NaJ^4\sqrt{\det\Omega_4}\), the divergence form of \(\Box_g\varphi=0\) is

\[
 \partial_t(J^4\Pi)=\partial_r\left[J^4\bigl(\beta\Pi+(1-\beta)\Phi\bigr)\right].
\]

Therefore the first-order time equations are

\[
 \boxed{\begin{aligned}
 \dot\varphi&=\beta\Phi+(1-\beta)\Pi,\\
 \dot\Phi&=\partial_r[\beta\Phi+(1-\beta)\Pi],\\
 \dot\Pi&=J^{-4}\partial_r\{J^4[\beta\Pi+(1-\beta)\Phi]\}-\frac{4\dot f}{J}\Pi.
 \end{aligned}}
\]

## ADM constraints and reconstruction

The Hamiltonian constraint \({}^{(5)}R+K^2-K_{ij}K^{ij}=16\pi G_6\rho\) becomes

\[
 \frac{8a'}{a^3J}+\frac{12}{J^2}(1-a^{-2})+8AB+12B^2
 =\frac{8\pi G_6}{a^2}(\Phi^2+\Pi^2).
\]

Multiplying by \(a^3J/8\),

\[
 \boxed{a'+\frac{3}{2J}(a^3-a)+a^3J(AB+\tfrac32B^2)
 -\pi G_6 aJ(\Phi^2+\Pi^2)=0.}
\]

The radial momentum constraint gives

\[
 -4B'+\frac{4(A-B)}{J}=-\frac{8\pi G_6}{a}\Phi\Pi,
\]

hence

\[
 \boxed{A=B+JB'-\frac{2\pi G_6J}{a}\Phi\Pi.}
\]

Here `A` is reconstructed algebraically from numerically differentiated `B` at each RK4 substage. The momentum residual is therefore *enforced*, not independently monitored.

## Evolution of metric and angular extrinsic curvature

The radial metric equation follows from \(\partial_t\gamma_{rr}=-2NK_{rr}+(\mathcal L_\beta\gamma)_{rr}\):

\[
 \boxed{\dot a=-a^2(1-\beta)A+(a\beta)'.}
\]

For the mixed-index extrinsic-curvature ADM equation, the massless-scalar trace-reversed Einstein equations imply

\[
 (\partial_t-\mathcal L_\beta)K^i{}_j
 =-D^iD_jN+N[{}^{(5)}R^i{}_j+KK^i{}_j-8\pi G_6 D^i\varphi D_j\varphi].
\]

The angular component has no direct spatial scalar gradient. Since \(D^\theta D_\theta N=N'/(a^2J)\),

\[
 \boxed{\dot B=\beta B'+a(1-\beta)B(A+4B)
 +\frac{3(1-\beta)}{J^2}(a-a^{-1})+\frac{\beta'}{aJ}.}
\]

The code evolves \(a,B,\varphi,\Phi,\Pi\) and the global coordinate shift \(f\), with `a` Hamiltonian-projected at intermediate stages.

## Horizon expansion and locking

An outgoing null normal to a spherical orbit is \(\ell^\mu=n^\mu+a^{-1}\delta^\mu_r\). Its expansion is

\[
 \theta_{(\ell)}=\frac4J\ell(J)=4\left(\frac1{aJ}-B\right).
\]

Thus the horizon condition is \(Q=aJB=1\). Fixing its computational coordinate \(r_h\) requires \(\partial_t Q(t,r_h)=0\). Substituting the evolution equations and both constraints at \(Q=1\) yields

\[
 \left.\partial_tQ\right|_h=\pi G_6 J_h(\Phi_h+\Pi_h)^2
 -\frac{3a_h^2}{2J_h}\dot f,
\]

so

\[
 \boxed{\dot f=\left.\frac{2\pi G_6}{3}\frac{J^2}{a^2}(\Phi+\Pi)^2\right|_h.}
\]

This sign is nonnegative for the usual scalar stress-energy tensor. Stationary vacuum gives \(\dot f=0\). The factor \(2\pi G_6/3\) replaces the original \(4\pi G_4\) coefficient in 3+1.

## Mass aspect and Hamiltonian projection

The generalized six-dimensional Misner--Sharp mass aspect is defined by

\[
 \boxed{\mu=\frac{J^3}{2}(1-\nabla_\mu J\nabla^\mu J)
 =\frac{J^3}{2}(1-a^{-2}+J^2B^2).}
\]

The physical quasi-local mass and horizon mass are

\[
 M_{\rm MS}=\frac{4\pi}{3G_6}\mu,\qquad
 \mu_h=\frac{J_h^3}{2},\qquad
 M_h=\frac{2\pi}{3G_6}J_h^3.
\]

Using the constraints to eliminate \(a'\) and \(B'\), obtain the radial mass relation

\[
 \boxed{\mu'=2\pi G_6J^4\left[\frac{\Phi^2+\Pi^2}{2a^2}
 +\frac{JB}{a}\Phi\Pi\right].}
\]

The constrained RK4 code integrates this relation from \(\mu(r_h)=J_h^3/2\) and sets

\[
 a^{-2}=1-2\mu/J^3+(JB)^2.
\]

This maintains the Hamiltonian constraint up to numerical quadrature/differentiation errors, but it does not establish constraint preservation of the full initial-boundary-value system.

## Exact Schwarzschild--Tangherlini benchmark

For \(q=(r_{h0}/J)^3\), \(f=0\), the ingoing EF-type slicing of the Tangherlini vacuum has

\[
 \boxed{a=\sqrt{1+q},\quad B=\frac{q}{J\sqrt{1+q}},\quad
 A=-\frac{3q(1+q/2)}{J(1+q)^{3/2}},\quad
 \beta=\frac{q}{1+q}.}
\]

The lapse is \(N=(1+q)^{-1/2}\), and \(\mu=r_{h0}^3/2\) is constant. Analytic substitution gives zero Hamiltonian and momentum residual. The vacuum balance routine subtracts the **discrete** vacuum RHS to maintain stationarity to floating-point precision, a separate numerical device whose truncation effects must be convergence-tested.

## Boundary and data limitations

The ingoing initial data approximates \(\partial_t(J^2\varphi)=\partial_r(J^2\varphi)\), appropriate to leading far-field asymptotics in five spatial dimensions, and uses a smooth cutoff to avoid scalar support at the horizon. At the outer boundary the code similarly imposes an **approximate** outgoing condition \((\partial_t+c_+\partial_r)(J^2\varphi)\simeq0\). Neither is an exact local spherical-wave boundary condition at finite radius in six dimensions. For strong matter, gauge collapse \(\dot f\geq1\), excision breakdown and outer reflections are possible. See `VALIDATION.md`.
