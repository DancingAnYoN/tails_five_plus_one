using Test
using MMIEFD6
using MMIEFD6.Types: IA,IB,IPHI,IPI,IFIELD
using MMIEFD6.Geometry: schwarzschild_Bprime,initial_mu
using MMIEFD6.Constraints: mass_rhs,project_hamiltonian!,hamiltonian_residual
using MMIEFD6.Evolution: EvolutionWorkspace,evolve_one_step!,rhs!,Workspace,setup_vacuum_balance!
using MMIEFD6.Observers: sample_fixed_J
using MMIEFD6.ScalarField: gaussian

@testset "5+1 MMIEF physical identities" begin
    rh=2.0; G6=1.0
    for J in (2.0,3.0,6.0,50.0)
        s=schwarzschild_fields(J,rh)
        q=(rh/J)^3
        @test isapprox(s.a,sqrt(1+q);rtol=1e-14)
        @test isapprox(s.beta,q/(1+q);rtol=1e-14)
        @test isapprox(s.a*J*s.B,q;rtol=1e-14)
        @test isapprox(misner_sharp_mass(J,s.a,s.B),rh^3/2;rtol=1e-12)
        @test isapprox(s.a*(1-s.beta),inv(s.a);rtol=1e-12)
        ap=-3*q/(2J*s.a)
        bp=schwarzschild_Bprime(J,rh)
        @test abs(hamiltonian_residual(J,s.a,ap,s.B,s.A,0.0,0.0,G6))<1e-10
        @test abs(bp+(s.B-s.A)/J)<1e-10
    end
    s=schwarzschild_fields(rh,rh)
    @test isapprox(s.a*rh*s.B,1.0;atol=1e-12)
    @test isapprox(physical_mass(rh^3/2,G6),2π*rh^3/(3G6);rtol=1e-12)
end

@testset "Initial values, constrained vacuum and horizon locking" begin
    p=Config(rmax=70.0,dr=0.5,tmax=0.2,pulse_amplitude=0.0)
    r,u=initial_data(p)
    for i in eachindex(r)
        s=schwarzschild_fields(r[i],p.rh0)
        @test isapprox(u[i,IA],s.a;atol=1e-11)
        @test isapprox(u[i,IB],s.B;atol=1e-11)
        @test u[i,IFIELD]==0.0
    end
    @test fdot_horizon(u,0.0,p)==0.0
    @test maximum(abs,vacuum_test(dr=0.5,rmax=70.0))<1e-11
    # Vacuum Hamiltonian projection leaves exact geometry unmodified.
    project_hamiltonian!(u,0.0,p,r)
    @test abs(u[1,IA]*p.rh0*u[1,IB]-1)<1e-10
    @test isapprox(misner_sharp_mass(r[end],u[end,IA],u[end,IB]),initial_mu(p);atol=1e-9)
    @test isapprox(sample_fixed_J(u,r,0.0,10.0,p),0.0;atol=1e-12)
    @test_throws ErrorException sample_fixed_J(u,r,100.0,10.0,p)
    u[1,IPHI]=2e-4;u[1,IPI]=1e-4
    exact=(2π*p.G6/3)*(p.rh0/u[1,IA])^2*(3e-4)^2
    @test isapprox(fdot_horizon(u,0.0,p),exact;rtol=1e-12)
end

@testset "Scalar and finite-radius output" begin
    p=Config(rmax=70.0,dr=0.5,cfl=0.1,tmax=0.2,
             output_dt=0.1,snapshot_dt=0.1,pulse_amplitude=1e-6,
             fixed_background=true)
    r,u=initial_data(p)
    @test maximum(abs,u[:,IFIELD])>0
    @test all(isfinite,u)
    @test maximum(abs,[misner_sharp_mass(r[i],u[i,IA],u[i,IB])-initial_mu(p)
                       for i in eachindex(r)])<1e-10
    sample_phi=sample_fixed_J(u,r,0.0,30.0,p)
    @test isfinite(sample_phi)
    dat=run_evolution(p;verbose=false)
    @test dat.settings.fixed_background
    @test dat.shifts[end]==0.0
    @test all(isfinite,dat.final_state)
    @test isapprox(dat.final_state[20,IA],u[20,IA];atol=1e-12)
    @test isapprox(dat.final_state[20,IB],u[20,IB];atol=1e-12)
    @test length(dat.times)>=2
    mktempdir() do path
        save_run(dat,path)
        loaded=load_run(path)
        @test length(loaded.times)==length(dat.times)
        @test loaded.settings.rh0==p.rh0
        @test loaded.settings.fixed_background
        @test isapprox(loaded.final_state[20,IA],dat.final_state[20,IA];rtol=1e-12)
        @test length(loaded.frame_profiles)==length(dat.frame_profiles)
    end
end

@testset "Diagnostic tail estimator" begin
    t=collect(10.0:1.0:400.0)
    y=1.2e-2 .* t.^(-10.0)
    slope=fit_tail(t,y;tlo=120.0,thi=350.0)
    @test slope!==nothing
    @test abs(slope.power-10.0)<1e-10
    @test slope.log_rms<1e-10
    lpi=local_power_index(t,y;peak_halfwidth=0.0,fit_halfwidth=15.0,tstart=120.0)
    @test abs(lpi[250]-10.0)<1e-9
end

@testset "Constrained six-dimensional smoke vacuum" begin
    p=Config(rmax=70.0,dr=0.5,cfl=0.1,tmax=0.2,
             output_dt=0.1,snapshot_dt=0.1,pulse_amplitude=0.0)
    run=run_evolution(p;verbose=false)
    @test maximum(abs,run.final_state[:,IPHI])<1e-12
    @test maximum(abs,run.final_state[:,IPI])<1e-12
    @test abs(run.shifts[end])<1e-12
    @test all(isfinite,run.final_state)
    @test isapprox(run.masses_h[end],run.masses_h[1];atol=1e-10)
end

@testset "Dedicated late-time tail extraction (synthetic analytic checks)" begin
    t=collect(10.0:1.0:360.0)
    scale=2.0
    y=0.125 .* (t ./ scale).^(-10)
    f=fit_power_law(t,y;tlo=140.0,thi=300.0,time_scale=scale)
    @test f !== nothing
    @test isapprox(f.power,10.0;atol=1e-10)
    @test isapprox(f.amplitude,0.125;atol=1e-10)
    @test f.log_rms < 1e-10
    @test isapprox(predict_power_law(f,220.0),0.125*(220/scale)^(-10);rtol=1e-10)
    e=amplitude_envelope(t,y;halfwidth=0.0)
    @test e==abs.(y)
    fe=fit_power_law(t,y;tlo=140.0,thi=300.0,time_scale=scale,
                     method=:envelope,envelope_halfwidth=0.0)
    @test isapprox(fe.power,10.0;atol=1e-10)
    @test fit_power_law(t,y;tlo=400.0,thi=500.0) === nothing
    @test_throws ArgumentError fit_power_law(t,y;tlo=0.0,thi=20.0)
    @test_throws ArgumentError fit_power_law(t,y;tlo=10.0,thi=20.0,method=:unknown)
    scans=scan_tail_windows(t,y;starts=[140.0,150.0],widths=[60.0,90.0],time_scale=scale)
    @test length(scans)==4
    @test all(abs(row.power-10.0)<1e-10 for row in scans)
    mktempdir() do dir
        summary=joinpath(dir,"tail.csv")
        scan=joinpath(dir,"scan.csv")
        save_tail_fits(summary,[(J=10.0,fit=f),(J=30.0,fit=nothing)])
        save_window_scan(scan,[merge((J=10.0,),row) for row in scans])
        @test length(readlines(summary))==2
        @test length(readlines(scan))==5
    end
end

@testset "Three-level observed waveform convergence (synthetic order two)" begin
    t=collect(20.0:2.0:180.0)
    function synthetic_run(h)
        cfg=Config(dr=h,fixed_background=true,
                   observers=[10.0,30.0,50.0])
        run=RunData(cfg,[2.0,3.0],zeros(2,5))
        append!(run.times,t)
        append!(run.shifts,zeros(length(t)))
        for (i,J) in enumerate(cfg.observers)
            # Shared continuum waveform plus a smooth leading O(dr^2) defect.
            cont=0.001 .* (t ./ 2).^(-4) .* (1 + J/100)
            defect=1e-7 .* cos.(t ./ 17) .* (t ./ 100).^(-2)
            append!(run.signals[i],cont .+ h^2 .* defect)
        end
        return run
    end
    a,b,c=(synthetic_run(h) for h in (0.5,0.25,0.125))
    report=compare_resolutions(a,b,c;tlo=50.0,thi=160.0,observer_index=2)
    @test isapprox(report.ratio,2.0;atol=1e-13)
    @test isapprox(report.observed_order,2.0;atol=1e-8)
    @test report.J==30.0
    @test all(isfinite,report.difference_cm)
    @test length(report.times)>12
    mktempdir() do dir
        save_convergence(joinpath(dir,"convergence.csv"),[report])
        @test length(readlines(joinpath(dir,"convergence.csv")))==2
    end
    @test_throws ArgumentError compare_resolutions(c,b,a;tlo=50.0,thi=160.0)
    different=synthetic_run(0.25)
    different.times .= t .+ 1000.0
    @test_throws ArgumentError compare_resolutions(a,different,c;tlo=50.0,thi=160.0)
end
