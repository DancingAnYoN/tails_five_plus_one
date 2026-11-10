#!/usr/bin/env python3
"""Independent optional SymPy checks of six-dimensional MMIEF identities.
Run: python scripts/check_continuum.py
Requires sympy (not required by Julia solver or tests).
"""
import sympy as s
J, rh, a, B, A, G, Phi, Pi, fdot, ap, bp = s.symbols(
    'J rh a B A G Phi Pi fdot ap bp', nonzero=True, real=True)
q=(rh/J)**3
av=s.sqrt(1+q)
Bv=q/(J*av)
Av=-3*q*(1+q/2)/(J*av**3)
H= s.diff(av,J)+s.Rational(3,2)/J*(av**3-av) + av**3*J*(Av*Bv+s.Rational(3,2)*Bv**2)
M=s.diff(Bv,J)+(Bv-Av)/J
mu=J**3/s.Integer(2)*(1-1/av**2+(J*Bv)**2)
assert s.simplify(H)==0
assert s.simplify(M)==0
assert s.simplify(mu-rh**3/2)==0
# Symbolically differentiate Q=a J B, then impose ADM equations/constraints at Q=1.
Q=a*J*B
beta=(fdot+Q)/(1+Q)
bprime=s.diff(beta,J)+s.diff(beta,a)*ap+s.diff(beta,B)*bp
adot=-a*a*(1-beta)*A+ap*beta+a*bprime
Bdot=beta*bp+a*(1-beta)*B*(A+4*B)+3*(1-beta)*(a-1/a)/J**2+bprime/(a*J)
Qdot=s.expand(adot*J*B + a*fdot*B + a*J*Bdot)
bpsub=(A-B)/J + 2*s.pi*G*Phi*Pi/a
apsub=-3*(a**3-a)/(2*J)-a**3*J*(A*B+s.Rational(3,2)*B**2)+s.pi*G*a*J*(Phi**2+Pi**2)
Qdot=s.factor(Qdot.subs({ap:apsub,bp:bpsub}).subs(B,1/(a*J)))
target=s.pi*G*J*(Phi+Pi)**2 - 3*a*a*fdot/(2*J)
assert s.simplify(Qdot-target)==0,(Qdot,target)
# Derivative of mass aspect from the constraints.
mu=J**3/s.Integer(2)*(1-a**-2+(J*B)**2)
mup=s.diff(mu,J)+s.diff(mu,a)*apsub+s.diff(mu,B)*bpsub
expected=2*s.pi*G*J**4*((Phi**2+Pi**2)/(2*a**2)+J*B*Phi*Pi/a)
assert s.simplify(mup-expected)==0
print('PASS: Tangherlini H/M, constant mass, MMIEF horizon-locking and mu prime.')
