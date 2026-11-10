module MMIEFD6TailVisualization

using MMIEFD6
using GLMakie
using Printf
export make_tail_diagnostics, make_tail_window_stability, make_tail_convergence

const PALETTE=[:royalblue,:darkorange,:seagreen,:purple,:firebrick]

function _valid_time(data)
    T=data.times .- data.shifts
    length(T)>1 && all(isfinite,T) && all(diff(T).>0) ||
        throw(ArgumentError("Observer times T=t-f must be finite and increasing"))
    return T
end

"""Plot |phi| against T on log-log axes with actual fit overlays, the local
power index, and log-amplitude residuals for every fixed-J observer.
The returned paths are PNG and CSV. Requires stored waveform samples.
"""
function make_tail_diagnostics(data::MMIEFD6.RunData;
        outdir::AbstractString="results", tlo::Real, thi::Real,
        method::Symbol=:absolute, envelope_halfwidth::Real=0.0,
        amplitude_floor::Real=0.0, peak_halfwidth::Real=0.0,
        lpi_halfwidth::Real=12.0, reference_power=nothing,
        show_window::Bool=false)
    GLMakie.activate!()
    T=_valid_time(data)
    rh=data.settings.rh0
    ref=reference_power===nothing && data.settings.fixed_background ? 10.0 : reference_power
    fits=MMIEFD6.observer_tail_fits(data;tlo=tlo,thi=thi,method=method,
         envelope_halfwidth=envelope_halfwidth,amplitude_floor=amplitude_floor)
    nobs=length(fits)
    nobs>0 || throw(ArgumentError("No fixed-J observers in this run"))
    fig=Figure(size=(1710,max(600,nobs*345+140)),fontsize=15)
    Label(fig[0,1:3],"5+1 MMIEF — fixed-J late-time tail extraction",fontsize=23,font=:bold)
    lpi_rows=Tuple{Float64,Float64,Float64}[]
    for (k,row) in enumerate(fits)
        J=row.J; fit=row.fit; y=data.signals[k];col=PALETTE[mod1(k,length(PALETTE))]
        status=fit===nothing ? "insufficient samples" :
          @sprintf("p=%.5f · RMS(log)=%.3g · N=%d",fit.power,fit.log_rms,fit.n)
        a1=Axis(fig[k,1],title=@sprintf("J=%.2f | %s",J,status),
                xlabel="log₁₀(T/rₕ₀)",ylabel="log₁₀ |φ|")
        a2=Axis(fig[k,2],title="Local power index",xlabel="T/rₕ₀",ylabel="p(T)")
        a3=Axis(fig[k,3],title="Log residual in fit interval",xlabel="T/rₕ₀",ylabel="ln |φ| − ln |fit|")
        valid=findall(i->isfinite(T[i]) && T[i]>0 && isfinite(y[i]) &&
                       abs(y[i])>max(amplitude_floor,0.0),eachindex(T))
        if !isempty(valid)
            lines!(a1,log10.(T[valid]./rh),log10.(abs.(y[valid]));
                   color=col,label="|φ| (stored)")
        end
        if method==:envelope
            env=MMIEFD6.amplitude_envelope(T,y;halfwidth=envelope_halfwidth)
            idx=findall(i->T[i]>0 && isfinite(env[i]) && env[i]>0,eachindex(T))
            !isempty(idx) && lines!(a1,log10.(T[idx]./rh),log10.(env[idx]);
                                     color=:gray,linestyle=:dot,label="max envelope")
        end
        vlines!(a1,log10.([tlo,thi]./rh);color=:gray,linestyle=:dash)
        lpi=MMIEFD6.local_power_index(T,y;peak_halfwidth=Float64(peak_halfwidth),
                   fit_halfwidth=Float64(lpi_halfwidth),tstart=Float64(tlo))
        for i in eachindex(lpi)
            isfinite(lpi[i]) && push!(lpi_rows,(J,Float64(T[i]),Float64(lpi[i])))
        end
        ii=findall(isfinite,lpi)
        if !isempty(ii)
            lines!(a2,T[ii]./rh,lpi[ii];color=col)
        end
        hlines!(a3,[0.0];color=:gray,linestyle=:dash)
        if ref!==nothing
            hlines!(a2,[Float64(ref)];color=:black,linestyle=:dash,label=@sprintf("reference p=%.2g",ref))
            axislegend(a2,position=:rt)
        end
        if fit!==nothing
            selected=findall(i->tlo<=T[i]<=thi && T[i]>0 && isfinite(y[i]) &&
                                abs(y[i])>max(amplitude_floor,0.0),eachindex(T))
            if !isempty(selected)
                lines!(a1,log10.(T[selected]./rh),
                    [ (fit.log_amplitude-fit.power*log(T[i]/rh))/log(10) for i in selected ];
                    color=:firebrick,linewidth=3,label="selected fit")
                lines!(a3,T[selected]./rh,
                    [log(abs(y[i]))-(fit.log_amplitude-fit.power*log(T[i]/rh)) for i in selected];
                    color=col)
            end
        end
        !isempty(valid) && axislegend(a1,position=:lb)
    end
    Label(fig[nobs+1,1:3],
          "T=t−f; fits use selected windows, not automatic asymptotic detection. The LPI is an envelope estimate and can be biased near zeros/noise. No convergence claim is implied.",
          fontsize=12,tellwidth=false)
    mkpath(outdir)
    figpath=joinpath(outdir,"tail_extraction_5plus1.png")
    csvpath=joinpath(outdir,"tail_fit_summary.csv")
    lpipath=joinpath(outdir,"tail_local_power_index.csv")
    save(figpath,fig)
    MMIEFD6.save_tail_fits(csvpath,fits)
    open(lpipath,"w") do io
        println(io,"J,T,T_over_rh0,local_power_index")
        for (J,time,power) in lpi_rows
            @printf(io,"%.17g,%.17g,%.17g,%.17g\n",J,time,time/rh,power)
        end
    end
    show_window && display(fig)
    println("Saved tail fit dashboard: ",figpath,"\nSaved tail fits and LPI: ",csvpath,", ",lpipath)
    return (figure=figpath,summary=csvpath,lpi=lpipath,fits=fits)
end

"""Scan manually specified late-time windows (without choosing an 'optimal'
window). Rows are observers, columns show fitted p and log-RMS against the
window start, and line colors indicate window width.
"""
function make_tail_window_stability(data::MMIEFD6.RunData;
        outdir::AbstractString="results", tlo::Real, thi::Real,
        starts=nothing, widths=nothing, nstarts::Int=16,
        method::Symbol=:absolute, envelope_halfwidth::Real=0.0,
        amplitude_floor::Real=0.0, reference_power=nothing,
        show_window::Bool=false)
    GLMakie.activate!()
    T=_valid_time(data);rh=data.settings.rh0
    tlo>0 && thi>tlo || throw(ArgumentError("Require 0 < tlo < thi"))
    nstarts>=2 || throw(ArgumentError("nstarts must be at least 2"))
    widths===nothing && (widths=[Float64(thi-tlo)])
    widths=Float64.(collect(widths))
    isempty(widths) && throw(ArgumentError("At least one scan width required"))
    minw=minimum(widths)
    starts===nothing && (starts=collect(range(Float64(tlo),
                                  stop=max(Float64(tlo),T[end]-minw),length=nstarts)))
    starts=Float64.(collect(starts))
    ref=reference_power===nothing && data.settings.fixed_background ? 10.0 : reference_power
    nobs=length(data.settings.observers)
    fig=Figure(size=(1330,max(590,nobs*310+140)),fontsize=15)
    Label(fig[0,1:2],"5+1 MMIEF — sensitivity to fit window",fontsize=22,font=:bold)
    rows=NamedTuple[]
    for k in 1:nobs
        J=data.settings.observers[k];y=data.signals[k]
        axp=Axis(fig[k,1],title=@sprintf("J=%.2f | fitted p",J),
                 xlabel="Window start T₀ / rₕ₀",ylabel="Power-law exponent p")
        axr=Axis(fig[k,2],title="Logarithmic fit residual",
                 xlabel="Window start T₀ / rₕ₀",ylabel="RMS ln amplitude")
        scan=MMIEFD6.scan_tail_windows(T,y;starts=starts,widths=widths,
                   time_scale=rh,method=method,
                   envelope_halfwidth=envelope_halfwidth,
                   amplitude_floor=amplitude_floor)
        for row in scan
            push!(rows,merge((J=J,),row))
        end
        for (widx,w) in enumerate(widths)
            vals=[row for row in scan if row.width==w]
            isempty(vals) && continue
            xs=[row.tlo/rh for row in vals]
            col=PALETTE[mod1(widx,length(PALETTE))]
            label=@sprintf("width=%.2f rₕ₀",w/rh)
            scatterlines!(axp,xs,[row.power for row in vals];color=col,marker=:circle,label=label)
            scatterlines!(axr,xs,[row.log_rms for row in vals];color=col,marker=:circle,label=label)
        end
        ref!==nothing && hlines!(axp,[Float64(ref)];color=:black,linestyle=:dash)
        !isempty(scan) && axislegend(axp,position=:rt)
    end
    Label(fig[nobs+1,1:2],
          "Overlapping windows are strongly correlated. A flat p(T₀) does not prove asymptotic convergence; vary grid spacing and boundary radius.",
          fontsize=12,tellwidth=false)
    mkpath(outdir)
    figpath=joinpath(outdir,"tail_window_stability_5plus1.png")
    csvpath=joinpath(outdir,"tail_window_scan.csv")
    save(figpath,fig)
    MMIEFD6.save_window_scan(csvpath,rows)
    show_window && display(fig)
    println("Saved window scan: ",figpath,"\nSaved scan table: ",csvpath)
    return (figure=figpath,scan=csvpath,rows=rows)
end

"""Compare three physics-matched resolution runs at fixed areal radii.
A log-amplitude difference plot and tail-power-vs-dr plot are saved alongside
one summary CSV and one interpolated waveform-difference CSV.
"""
function make_tail_convergence(coarse::MMIEFD6.RunData,
        middle::MMIEFD6.RunData,fine::MMIEFD6.RunData;
        outdir::AbstractString="results/tail_convergence",tlo::Real,thi::Real,
        method::Symbol=:absolute, envelope_halfwidth::Real=0.0,
        amplitude_floor::Real=0.0, reference_power=nothing,
        show_window::Bool=false)
    GLMakie.activate!()
    rh=coarse.settings.rh0
    ref=reference_power===nothing && coarse.settings.fixed_background ? 10.0 : reference_power
    nobs=length(coarse.settings.observers)
    results=[MMIEFD6.compare_resolutions(coarse,middle,fine;
                tlo=tlo,thi=thi,observer_index=k,method=method,
                envelope_halfwidth=envelope_halfwidth,
                amplitude_floor=amplitude_floor)
             for k in 1:nobs]
    fig=Figure(size=(1390,max(580,300*nobs+140)),fontsize=15)
    Label(fig[0,1:2],"5+1 MMIEF — three-resolution waveform comparison",fontsize=22,font=:bold)
    for (k,res) in enumerate(results)
        status=isfinite(res.observed_order) ? @sprintf("apparent order q=%.3f",res.observed_order) : "order undefined"
        a=Axis(fig[k,1],xlabel="T/rₕ₀",ylabel="log₁₀ |difference|",
               title=@sprintf("J=%.2f | %s",res.J,status))
        b=Axis(fig[k,2],xlabel="Δr/rₕ₀ (coarse → fine)",ylabel="Fitted p",
               title="Exponent across resolutions")
        for (y,label,col) in ((res.difference_cm,"coarse − medium",:royalblue),
                              (res.difference_mf,"medium − fine",:darkorange))
            vals=log10.(max.(abs.(y),1e-300))
            lines!(a,res.times./rh,vals;color=col,label=label)
        end
        axislegend(a,position=:rt)
        for (idx,fit) in enumerate(res.tail_fits)
            fit===nothing && continue
            scatter!(b,[res.dr[idx]/rh],[fit.power];markersize=14,
                     color=PALETTE[idx],label=string("Δr=",@sprintf("%.4g",res.dr[idx])))
        end
        ref!==nothing && hlines!(b,[Float64(ref)];color=:black,linestyle=:dash)
        any(fit!==nothing for fit in res.tail_fits) && axislegend(b,position=:rt)
    end
    Label(fig[nobs+1,1:2],
          "Apparent order from L2 waveform differences on shared coarse T samples; time interpolation and boundary effects are NOT removed.",
          fontsize=12,tellwidth=false)
    mkpath(outdir)
    figpath=joinpath(outdir,"tail_convergence_5plus1.png")
    summary=joinpath(outdir,"tail_convergence_summary.csv")
    curves=joinpath(outdir,"tail_convergence_curves.csv")
    save(figpath,fig)
    MMIEFD6.save_convergence(summary,results)
    open(curves,"w") do io
        println(io,"J,T,phi_coarse,phi_medium,phi_fine,diff_coarse_medium,diff_medium_fine")
        for res in results
            for k in eachindex(res.times)
                @printf(io,"%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    res.J,res.times[k],res.wave_coarse[k],res.wave_middle[k],
                    res.wave_fine[k],res.difference_cm[k],res.difference_mf[k])
            end
        end
    end
    show_window && display(fig)
    println("Saved convergence figure: ",figpath,"\nSaved CSV reports: ",summary,", ",curves)
    return (figure=figpath,summary=summary,curves=curves,results=results)
end

end # module
