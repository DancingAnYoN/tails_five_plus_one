module EinsteinEquations
using ..Types: IA,IB
export gravity_rhs!
"""5+1 mixed angular ADM evolution and radial metric evolution (J'=1)."""
function gravity_rhs!(du,u,J,A,beta,dB,dBeta,dABeta)
    @inbounds for i in eachindex(J)
        a=u[i,IA]; B=u[i,IB]; Ai=A[i]; b=beta[i]; Ji=J[i]
        du[i,IA]=-a^2*(1-b)*Ai+dABeta[i]
        du[i,IB]=b*dB[i]+a*(1-b)*B*(Ai+4B) +
                 3*(1-b)*(a-1/a)/Ji^2+dBeta[i]/(a*Ji)
    end
    return nothing
end
end
