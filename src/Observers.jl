module Observers
using Printf
using ..Types: Config, IFIELD
using ..Geometry: fdot_horizon
export sample_fixed_J
"""Cubic interpolation at a fixed physical areal radius, J_obs=r+f.
Never silently clamp a fixed-J observer that has left the grid.
"""
function sample_fixed_J(u::Matrix{Float64}, r::Vector{Float64}, f::Float64,
                        Jobs::Float64, p::Config; time::Float64=NaN)
    x=Jobs-f
    if !(isfinite(x) && r[1]<=x<=r[end])
        error(@sprintf("Fixed-J observer left grid: t=%.6g, J_obs=%.6g, f=%.6g, r_obs=%.6g, domain=[%.6g, %.6g], fdot=%.6g",
                       time,Jobs,f,x,r[1],r[end],fdot_horizon(u,f,p)))
    end
    n=length(r)
    i0=clamp(floor(Int,(x-r[1])/p.dr),1,n-3)
    ans=0.0
    @inbounds for j in i0:i0+3
        w=1.0
        for k in i0:i0+3
            if k!=j
                w*=(x-r[k])/(r[j]-r[k])
            end
        end
        ans+=w*u[j,IFIELD]
    end
    return ans
end
end
