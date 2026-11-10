module TailExtraction

using Printf
using ..DataOutput: RunData
export fit_power_law, predict_power_law, amplitude_envelope,
       scan_tail_windows, observer_tail_fits, compare_resolutions,
       save_tail_fits, save_window_scan, save_convergence

"""Running maximum of |phi| over a symmetric *coordinate-time* interval.
Pure power laws are slightly biased by nonzero window width; use the `:absolute`
method when the signal has no late-time sign changes.
"""
function amplitude_envelope(t::AbstractVector{<:Real}, y::AbstractVector{<:Real};
                            halfwidth::Real=0.0)
    length(t)==length(y) || throw(ArgumentError("Time and field lengths differ"))
    halfwidth>=0 || throw(ArgumentError("halfwidth must be nonnegative"))
    issorted(t) || throw(ArgumentError("Time samples must be sorted"))
    n=length(t)
    a=Vector{Float64}(undef,n)
    for i in eachindex(t)
        lo=searchsortedfirst(t,t[i]-halfwidth)
        hi=searchsortedlast(t,t[i]+halfwidth)
        # Ignore non-finite samples instead of letting one NaN contaminate a window.
        good=(abs(Float64(y[k])) for k in lo:hi if isfinite(y[k]))
        value=maximum(good;init=-Inf)
        a[i]=value==-Inf ? NaN : value
    end
    return a
end

"""Fit |phi| = A*(T/time_scale)^(-p) on [tlo,thi] using log-space OLS.

`method=:absolute` fits nonzero |phi| samples. `method=:envelope` fits a
running-maximum envelope, useful near oscillatory zero crossings but potentially
biased for short windows. `amplitude_floor` is an *absolute* data-unit cutoff.
Returns `nothing` on inadequate data, otherwise a named tuple including the
amplitude, exponent, residual RMS in natural-log units, and R². This is a
point estimate, not a physical uncertainty interval.
"""
function fit_power_law(t::AbstractVector{<:Real}, y::AbstractVector{<:Real};
                       tlo::Real, thi::Real, time_scale::Real=1.0,
                       method::Symbol=:absolute, envelope_halfwidth::Real=0.0,
                       amplitude_floor::Real=0.0, min_samples::Int=12)
    length(t)==length(y) || throw(ArgumentError("Time and field lengths differ"))
    tlo>0 && thi>tlo || throw(ArgumentError("Require 0 < tlo < thi"))
    time_scale>0 || throw(ArgumentError("time_scale must be positive"))
    amplitude_floor>=0 || throw(ArgumentError("amplitude_floor must be nonnegative"))
    min_samples>=3 || throw(ArgumentError("min_samples must be >= 3"))
    method in (:absolute,:envelope) || throw(ArgumentError("method must be :absolute or :envelope"))
    issorted(t) || throw(ArgumentError("Time samples must be sorted"))
    a=method==:absolute ? abs.(Float64.(y)) :
      amplitude_envelope(t,y;halfwidth=envelope_halfwidth)
    ids=[i for i in eachindex(t) if tlo<=t[i]<=thi && isfinite(t[i]) &&
         isfinite(a[i]) && a[i]>max(amplitude_floor,0.0)]
    length(ids)>=min_samples || return nothing
    x=log.(Float64.(t[ids])./time_scale)
    z=log.(a[ids])
    xm=sum(x)/length(x); zm=sum(z)/length(z)
    sxx=sum((u-xm)^2 for u in x)
    sxx>1e-16 || return nothing
    slope=sum((x[k]-xm)*(z[k]-zm) for k in eachindex(x))/sxx
    intercept=zm-slope*xm
    ssr=sum((z[k]-intercept-slope*x[k])^2 for k in eachindex(x))
    sst=sum((zk-zm)^2 for zk in z)
    return (power=-slope, log_amplitude=intercept, amplitude=exp(intercept),
            log_rms=sqrt(ssr/length(x)), r2=sst>0 ? 1-ssr/sst : NaN,
            n=length(x), tlo=Float64(tlo), thi=Float64(thi),
            time_scale=Float64(time_scale), method=method)
end

"""Power-law model evaluated at positive physical observer time T."""
function predict_power_law(fit, T::Real)
    T>0 || return NaN
    return exp(fit.log_amplitude-fit.power*log(T/fit.time_scale))
end

"""Compute manual, rather than automatically selected, overlapping fit windows.
Each output row specifies the observer, [tlo,thi], width, fitted exponent,
log-RMS, R², and sample count. Windows with insufficient samples are omitted.
"""
function scan_tail_windows(t::AbstractVector{<:Real}, y::AbstractVector{<:Real};
                           starts::AbstractVector{<:Real}, widths::AbstractVector{<:Real},
                           time_scale::Real=1.0, method::Symbol=:absolute,
                           envelope_halfwidth::Real=0.0,
                           amplitude_floor::Real=0.0, min_samples::Int=12)
    out=NamedTuple[]
    for width in widths
        width>0 || throw(ArgumentError("Scan widths must be positive"))
        for tlo in starts
            tlo>0 || continue
            thi=tlo+width
            fit=fit_power_law(t,y;tlo=tlo,thi=thi,time_scale=time_scale,
                              method=method,envelope_halfwidth=envelope_halfwidth,
                              amplitude_floor=amplitude_floor,min_samples=min_samples)
            fit===nothing && continue
            push!(out,(tlo=Float64(tlo),thi=Float64(thi),width=Float64(width),
                       power=fit.power,log_rms=fit.log_rms,r2=fit.r2,n=fit.n))
        end
    end
    return out
end

function observer_tail_fits(data::RunData; tlo::Real, thi::Real,
                            method::Symbol=:absolute, envelope_halfwidth::Real=0.0,
                            amplitude_floor::Real=0.0, min_samples::Int=12)
    T=data.times .- data.shifts  # fixed-areal-radius observer time
    return [(J=J,fit=fit_power_law(T,y;tlo=tlo,thi=thi,
                    time_scale=data.settings.rh0,method=method,
                    envelope_halfwidth=envelope_halfwidth,
                    amplitude_floor=amplitude_floor,min_samples=min_samples))
            for (J,y) in zip(data.settings.observers,data.signals)]
end

function _interp_linear(x::Vector{Float64}, y::Vector{Float64}, q::Float64)
    # All times should already be within [first(x),last(x)]. No extrapolation.
    first(x)<=q<=last(x) || throw(ArgumentError("Interpolation outside saved time range"))
    k=searchsortedlast(x,q)
    k==length(x) && return y[end]
    k>=1 || throw(ArgumentError("Interpolation before first time sample"))
    h=x[k+1]-x[k]
    h>0 || throw(ArgumentError("Repeated/non-increasing observer time"))
    z=(q-x[k])/h
    return (1-z)*y[k]+z*y[k+1]
end

function _validate_resolution_runs(a::RunData,b::RunData,c::RunData)
    p=a.settings
    for other in (b,c)
        q=other.settings
        # Only the grid/time resolution and output cadence may differ.
        for field in (:rh0,:G6,:rmax,:pulse_r0,:pulse_sigma,
                      :pulse_amplitude,:cutoff_start,:cutoff_end,:fixed_background,
                      :ko_strength)
            isapprox(getfield(p,field),getfield(q,field);rtol=1e-12,atol=1e-14) ||
                throw(ArgumentError("Runs have different $(field); cannot claim a resolution comparison"))
        end
        length(p.observers)==length(q.observers) &&
        all(isapprox.(p.observers,q.observers;rtol=1e-12,atol=1e-12)) ||
            throw(ArgumentError("Runs have different fixed-J observers"))
    end
    h=(p.dr,b.settings.dr,c.settings.dr)
    h[1]>h[2]>h[3]>0 || throw(ArgumentError("Supply runs in coarse, medium, fine order"))
    r1=h[1]/h[2];r2=h[2]/h[3]
    isapprox(r1,r2;rtol=0.02) ||
        throw(ArgumentError("Grid spacing ratios must agree within 2% (got $r1, $r2)"))
    return (h=h,ratio=(r1+r2)/2)
end

"""Three-resolution observed waveform self-convergence at the *same fixed J*.

All three runs must use matched physics, observers and approximately geometric
mesh refinement. Samples at coarse-run observer times are interpolated from
medium/fine runs, within a shared [tlo,thi] only. Reports L2 waveform
*difference* norms and apparent order q = ln(||c-m||/||m-f||)/ln(r).
Do not interpret q alone as evidence of physical tail or continuum convergence;
check output cadence, outer boundary, and fit windows independently.
"""
function compare_resolutions(coarse::RunData, middle::RunData, fine::RunData;
                             tlo::Real, thi::Real, observer_index::Int=1,
                             method::Symbol=:absolute,
                             envelope_halfwidth::Real=0.0,
                             amplitude_floor::Real=0.0)
    x=_validate_resolution_runs(coarse,middle,fine)
    tlo>0 && thi>tlo || throw(ArgumentError("Require 0 < tlo < thi"))
    1<=observer_index<=length(coarse.signals) ||
        throw(ArgumentError("Invalid observer_index"))
    runs=(coarse,middle,fine)
    T=map(run->run.times .- run.shifts,runs)
    for tt in T
        length(tt)>=3 && all(isfinite,tt) && all(diff(tt).>0) ||
            throw(ArgumentError("Observer time T must be finite and strictly increasing"))
    end
    lo=max(Float64(tlo),maximum(first(tt) for tt in T))
    hi=min(Float64(thi),minimum(last(tt) for tt in T))
    hi>lo || throw(ArgumentError("No common observer-time overlap in requested window"))
    times=[Float64(t) for t in T[1] if lo<=t<=hi]
    length(times)>=12 || throw(ArgumentError("Need at least 12 coarse-run samples in common window"))
    waves=map(eachindex(runs)) do k
        [_interp_linear(T[k],runs[k].signals[observer_index],t) for t in times]
    end
    u,v,w=waves
    cm=u.-v; mf=v.-w
    ecm=sqrt(sum(abs2,cm)/length(cm))
    emf=sqrt(sum(abs2,mf)/length(mf))
    order=(ecm>0 && emf>0) ? log(ecm/emf)/log(x.ratio) : NaN
    fits=map(1:3) do k
        fit_power_law(T[k],runs[k].signals[observer_index];
                      tlo=lo,thi=hi,time_scale=coarse.settings.rh0,
                      method=method,envelope_halfwidth=envelope_halfwidth,
                      amplitude_floor=amplitude_floor)
    end
    return (J=coarse.settings.observers[observer_index], tlo=lo, thi=hi,
            dr=x.h, ratio=x.ratio, times=times, wave_coarse=u,
            wave_middle=v,wave_fine=w,difference_cm=cm,difference_mf=mf,
            norm_cm=ecm,norm_mf=emf,observed_order=order,
            tail_fits=fits)
end

function save_tail_fits(path::AbstractString, rows)
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,"J,tlo,thi,power,amplitude,log_rms,r2,n,method,time_scale")
        for row in rows
            f=row.fit
            f===nothing && continue
            @printf(io,"%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%s,%.17g\n",
                    row.J,f.tlo,f.thi,f.power,f.amplitude,f.log_rms,f.r2,
                    f.n,string(f.method),f.time_scale)
        end
    end
    return path
end

function save_window_scan(path::AbstractString, rows)
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,"J,tlo,thi,width,power,log_rms,r2,n")
        for row in rows
            @printf(io,"%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d\n",
                    row.J,row.tlo,row.thi,row.width,row.power,row.log_rms,row.r2,row.n)
        end
    end
    return path
end

function save_convergence(path::AbstractString, results)
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,"J,tlo,thi,dr_coarse,dr_medium,dr_fine,refinement_ratio,L2_coarse_medium,L2_medium_fine,apparent_order,p_coarse,p_medium,p_fine")
        for result in results
            pow=[f===nothing ? NaN : f.power for f in result.tail_fits]
            @printf(io,"%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    result.J,result.tlo,result.thi,result.dr[1],result.dr[2],result.dr[3],
                    result.ratio,result.norm_cm,result.norm_mf,result.observed_order,
                    pow[1],pow[2],pow[3])
        end
    end
    return path
end

end # module TailExtraction
