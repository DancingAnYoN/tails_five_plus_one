module MMIEFD6Visualization
using MMIEFD6
using GLMakie
using Printf
export make_dashboard,make_pulse_movie,make_observer_history

@inline point(x,y)=Point2f(Float32(x),Float32(y))

"""Render a 5+1 dashboard from stored data only (does not rerun the solver)."""
function make_dashboard(data::MMIEFD6.RunData;
                        outdir::AbstractString="results",show_window::Bool=false)
    GLMakie.activate!()
    mkpath(outdir)
    p=data.settings;rh=p.rh0
    T=data.times .- data.shifts
    fig=Figure(size=(1510,1040),fontsize=17)
    Label(fig[0,1:2],"5+1 Einstein--Klein--Gordon | MMIEF",fontsize=25,font=:bold)
    ax1=Axis(fig[1,1],title="Last scalar profile",xlabel="J / rₕ₀",ylabel="J² φ")
    ax2=Axis(fig[1,2],title="Mass aspect",xlabel="t / rₕ₀",ylabel="μ / μ₀")
    ax3=Axis(fig[2,1],title="Fixed-J waveform history",xlabel="T / rₕ₀",ylabel="φ")
    ax4=Axis(fig[2,2],title="Tail amplitude (not a verified power law)",xlabel="log₁₀(T/rₕ₀)",ylabel="log₁₀ |φ|")
    ax5=Axis(fig[3,1],title="Local power index",xlabel="T/rₕ₀",ylabel="-d ln(env)/d ln T")
    ax6=Axis(fig[3,2],title="Constraint diagnostics",xlabel="t/rₕ₀",ylabel="log₁₀ residual")
    if !isempty(data.frame_profiles)
        f=data.frame_shift[end]
        rgrid=data.r[1:p.profile_stride:end]
        scatterX=(rgrid .+ f)./rh
        lines!(ax1,scatterX,Float64.(data.frame_profiles[end]);color=:royalblue)
        vlines!(ax1,[(p.rh0+f)/rh],color=:firebrick,linestyle=:dash,label="Apparent horizon")
    end
    mu0=rh^3/2
    muH=data.masses_h .* (3*p.G6/(4π))
    muO=data.masses_out .* (3*p.G6/(4π))
    lines!(ax2,data.times./rh,muH./mu0,label="Horizon",color=:firebrick)
    lines!(ax2,data.times./rh,muO./mu0,label="Outer",color=:slateblue,linestyle=:dash)
    axislegend(ax2,position=:rb)
    colors=[:dodgerblue,:darkorange,:seagreen,:purple,:firebrick]
    for (k,J) in enumerate(p.observers)
        phi=data.signals[k]
        label=@sprintf("J = %.1f",J)
        col=colors[mod1(k,length(colors))]
        lines!(ax3,T./rh,phi,label=label,color=col)
        keep=findall(i->T[i]>0 && abs(phi[i])>1e-299,eachindex(T))
        lines!(ax4,log10.(T[keep]./rh),log10.(abs.(phi[keep])),color=col,label=label)
        lpi=MMIEFD6.local_power_index(T,phi;
                     peak_halfwidth=12rh,fit_halfwidth=24rh,tstart=50rh)
        finite=findall(isfinite,lpi)
        lines!(ax5,T[finite]./rh,lpi[finite],color=col,label=label)
    end
    hlines!(ax5,[10.0],color=:black,linestyle=:dash,label="linear reference p=10")
    axislegend(ax3,position=:rt)
    axislegend(ax4,position=:lb)
    lines!(ax6,data.times./rh,log10.(max.(data.ham_l2,1e-30)),label="Hamiltonian",color=:royalblue)
    lines!(ax6,data.times./rh,log10.(max.(data.scalar_constraint_l2,1e-30)),
           label="Scalar definition",color=:seagreen)
    axislegend(ax6,position=:rt)
    Label(fig[4,1:2],"Fixed: linear test field (H residual includes scalar source). Coupled: constrained research prototype. No tail exponent is claimed without convergence tests.",
          fontsize=13,tellwidth=false)
    out=joinpath(outdir,"dashboard.png")
    save(out,fig)
    println("Saved 5+1 dashboard: ",out)
    show_window && display(fig)
    return out
end

"""Snapshot MP4: J² phi, areal radius J, and moving apparent horizon."""
function make_pulse_movie(data::MMIEFD6.RunData;
                          outdir::AbstractString="results",fps::Int=24,
                          every::Int=1,radial_max::Real=100.0,
                          show_window::Bool=false)
    GLMakie.activate!()
    mkpath(outdir)
    nt=length(data.frame_times)
    nt>=2 || error("Need at least two saved snapshots")
    every>=1 && fps>=1 || error("Invalid animation parameters")
    p=data.settings;rh=p.rh0
    frames=collect(1:every:nt)
    frames[end]!=nt && push!(frames,nt)
    rgrid=data.r[1:p.profile_stride:end]
    scale=max(p.pulse_amplitude*rh^2,eps())
    maxamp=max(maximum(maximum(abs,prof) for prof in data.frame_profiles)/scale,1e-6)
    fig=Figure(size=(1360,800),fontsize=19)
    ttl=Observable("5+1 MMIEF scalar pulse")
    Label(fig[0,1:2],ttl,fontsize=24,font=:bold)
    ax1=Axis(fig[1:2,1],xlabel="J / rₕ₀",ylabel="J²φ / (rₕ₀² ε)",title="Scalar profile")
    ax2=Axis(fig[1,2],xlabel="t / rₕ₀",ylabel="μₕ / μ₀",title="Apparent horizon mass aspect")
    ax3=Axis(fig[2,2],xlabel="T / rₕ₀",ylabel="φ(T,Jₒᵦₛ)",title="Stored observer signals")
    xs=Observable(Float32[]);ys=Observable(Float32[]);hs=Observable(Float32[p.rh0/rh])
    lines!(ax1,xs,ys,color=:royalblue,linewidth=3)
    vlines!(ax1,hs,color=:firebrick,linestyle=:dash)
    hlines!(ax1,[0.0],color=:gray)
    xlims!(ax1,0,min(radial_max/rh,(data.r[end]+minimum(data.frame_shift))/rh))
    ylims!(ax1,-1.15maxamp,1.15maxamp)
    mu0=rh^3/2
    muH=data.masses_h .* (3*p.G6/(4π))
    hts=Observable(Float32[]);hvals=Observable(Float32[])
    lines!(ax2,hts,hvals,color=:firebrick,linewidth=3)
    xlims!(ax2,0,max(data.times[end]/rh,1e-6))
    ylims!(ax2,max(0.0,minimum(muH)/mu0-0.02),maximum(muH)/mu0+max(0.01,0.1*(maximum(muH)-minimum(muH))/mu0))
    T=data.times .- data.shifts
    obx=[Observable(Float32[]) for _ in p.observers]
    oby=[Observable(Float32[]) for _ in p.observers]
    cols=[:royalblue,:darkorange,:seagreen,:purple]
    for j in eachindex(p.observers)
        lines!(ax3,obx[j],oby[j],label=@sprintf("J=%.0f",p.observers[j]),color=cols[mod1(j,length(cols))])
    end
    axislegend(ax3,position=:rt)
    function update!(k::Int)
        t=data.frame_times[k]; f=data.frame_shift[k]
        xs[]=Float32.((rgrid .+ f)./rh)
        ys[]=Float32.(data.frame_profiles[k]./scale)
        hs[]=Float32[(p.rh0+f)/rh]
        idx=clamp(searchsortedlast(data.times,t+1e-10),1,length(data.times))
        hts[]=Float32.(data.times[1:idx]./rh)
        hvals[]=Float32.(muH[1:idx]./mu0)
        for j in eachindex(obx)
            obx[j][]=Float32.(T[1:idx]./rh)
            oby[j][]=Float32.(data.signals[j][1:idx])
        end
        ttl[]=@sprintf("5+1 MMIEF | t/rh0=%.2f | Jh/rh0=%.6f",t/rh,(p.rh0+f)/rh)
    end
    update!(1)
    out=joinpath(outdir,"scalar_pulse_5plus1.mp4")
    record(fig,out,frames;framerate=fps) do k
        update!(k)
    end
    update!(nt)
    save(joinpath(outdir,"scalar_pulse_final.png"),fig)
    show_window && display(fig)
    println("Saved 5+1 pulse movie: ",out)
    return out
end

function make_observer_history(data::MMIEFD6.RunData;
                               outdir::AbstractString="results",show_window::Bool=false)
    GLMakie.activate!();mkpath(outdir)
    T=data.times .- data.shifts;rh=data.settings.rh0
    fig=Figure(size=(1340,840),fontsize=17)
    a1=Axis(fig[1,1],xlabel="T/rh0",ylabel="φ",title="Full fixed-J history")
    a2=Axis(fig[1,2],xlabel="T/rh0",ylabel="log10 |φ|",title="Log amplitude history")
    for (k,J) in enumerate(data.settings.observers)
        phi=data.signals[k]
        label=@sprintf("J=%.0f",J)
        lines!(a1,T./rh,phi,label=label)
        lines!(a2,T./rh,log10.(max.(abs.(phi),1e-300)),label=label)
    end
    axislegend(a1,position=:rt)
    axislegend(a2,position=:rt)
    out=joinpath(outdir,"observer_history_5plus1.png")
    save(out,fig)
    show_window && display(fig)
    println("Saved fixed-J history: ",out)
    return out
end
end
