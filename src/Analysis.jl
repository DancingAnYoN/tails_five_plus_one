module Analysis
export fit_qnm, local_power_index, fit_tail

"""Exploratory single damped-sinusoid fit. Frequencies scale with the
user-supplied *length* scale, NOT physical ADM mass in 5+1.
Fails if window lacks enough oscillations; not a QNM certification.
"""
function fit_qnm(t::Vector{Float64},y::Vector{Float64};
                 tlo::Float64,thi::Float64,length_scale::Float64=1.0,
                 omega_range::Tuple{Float64,Float64}=(0.04,2.0),
                 gamma_range::Tuple{Float64,Float64}=(0.01,1.0),
                 nsearch::Int=50)
    length(t)==length(y) || throw(ArgumentError("Inconsistent time-series lengths"))
    length_scale>0 || throw(ArgumentError("length_scale must be positive"))
    idx=findall(i->tlo<=t[i]<=thi,eachindex(t))
    length(idx)>=25 || return nothing
    tt=t[idx] .- tlo; yy=y[idx]
    eyy=sum(abs2,yy)
    eyy>eps(Float64) || return nothing
    best=(Inf,0.0,0.0,0.0,0.0)
    for gammaL in range(gamma_range[1],gamma_range[2],length=nsearch),
        omegaL in range(omega_range[1],omega_range[2],length=nsearch)
        g=gammaL/length_scale;w=omegaL/length_scale
        scc=0.0;sss=0.0;scs=0.0;syc=0.0;sys=0.0
        for k in eachindex(tt)
            e=exp(-g*tt[k]);C=e*cos(w*tt[k]);S=e*sin(w*tt[k])
            scc+=C*C;sss+=S*S;scs+=C*S;syc+=C*yy[k];sys+=S*yy[k]
        end
        den=scc*sss-scs*scs
        den<=1e-25 && continue
        C=(syc*sss-sys*scs)/den;S=(sys*scc-syc*scs)/den
        err=sum((yy[k]-exp(-g*tt[k])*(C*cos(w*tt[k])+S*sin(w*tt[k])))^2
                for k in eachindex(tt))
        if err<best[1]
            best=(err,gammaL,omegaL,C,S)
        end
    end
    isfinite(best[1]) || return nothing
    err,gl,wl,C,S=best
    return (gammaL=gl,omegaL=wl,C=C,S=S,tlo=tlo,thi=thi,
            length_scale=length_scale,relative_fit_error=sqrt(err/eyy))
end

"""Running maximum of |phi| within +/- peak_halfwidth, followed by local
least-squares slope in log(time), returned as p=-d log|phi|/d log(time).
"""
function local_power_index(t::Vector{Float64},y::Vector{Float64};
                           peak_halfwidth::Float64=12.0,
                           fit_halfwidth::Float64=24.0,tstart::Float64=100.0)
    n=length(t); n==length(y) || error("Incompatible t and y")
    issorted(t) || error("Time series must be sorted")
    n==0 && return Float64[]
    env=zeros(n)
    for i in 1:n
        lo=searchsortedfirst(t,t[i]-peak_halfwidth)
        hi=searchsortedlast(t,t[i]+peak_halfwidth)
        env[i]=maximum(abs,view(y,lo:hi))
    end
    out=fill(NaN,n)
    flooramp=max(maximum(env)*1e-22,1e-29)
    for i in eachindex(t)
        (t[i]>=tstart && env[i]>=flooramp) || continue
        lo=searchsortedfirst(t,max(tstart,t[i]-fit_halfwidth))
        hi=searchsortedlast(t,t[i]+fit_halfwidth)
        hi-lo>=12 || continue
        sx=0.0;sy=0.0;sxx=0.0;sxy=0.0;count=0
        for j in lo:hi
            if t[j]>0 && env[j]>=flooramp
                x=log(t[j]);z=log(env[j])
                sx+=x;sy+=z;sxx+=x*x;sxy+=x*z;count+=1
            end
        end
        den=count*sxx-sx*sx
        if count>12 && den>1e-12
            out[i]=-(count*sxy-sx*sy)/den
        end
    end
    return out
end

"""OLS decay index p in |phi|~T^(-p), absolute-value fit.
Use only in a late-time window with no unresolved oscillations/zero crossings.
"""
function fit_tail(t::Vector{Float64},y::Vector{Float64};
                  tlo::Float64,thi::Float64)
    length(t)==length(y) || error("Time/amplitude length mismatch")
    ids=[i for i in eachindex(t) if tlo<=t[i]<=thi && t[i]>0 &&
         isfinite(y[i]) && abs(y[i])>1e-300]
    length(ids)>=10 || return nothing
    xx=log.(t[ids]);zz=log.(abs.(y[ids]))
    xm=sum(xx)/length(xx);zm=sum(zz)/length(zz)
    den=sum((x-xm)^2 for x in xx)
    den>0 || return nothing
    slope=sum((x-xm)*(z-zm) for (x,z) in zip(xx,zz))/den
    resid=sqrt(sum((z-zm-slope*(x-xm))^2 for (x,z) in zip(xx,zz))/length(xx))
    return (power=-slope,log_rms=resid,n=length(ids),tlo=tlo,thi=thi)
end
end
