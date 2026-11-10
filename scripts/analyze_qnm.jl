# Exploratory damped sinusoid, no built-in 4D frequency benchmark.
# julia --project=. scripts/analyze_qnm.jl results/quick 20 60
using MMIEFD6
using Printf
function main(args::Vector{String}=ARGS)
    length(args)==3 || error("Usage: analyze_qnm.jl RUN_DIR TLO THI")
    run=load_run(args[1]);tlo=parse(Float64,args[2]);thi=parse(Float64,args[3])
    T=run.times .- run.shifts
    for (J,phi) in zip(run.settings.observers,run.signals)
        f=fit_qnm(T,phi;tlo=tlo,thi=thi,length_scale=run.settings.rh0)
        if f===nothing
            @printf("J=%.2f: no stable single-mode fit available\n",J)
        else
            @printf("J=%.2f | rh0 omega=%.5g | rh0 Gamma=%.5g | rel RMS=%.4g\n",
                J,f.omegaL,f.gammaL,f.relative_fit_error)
        end
    end
    println("Exploratory only: window, radius, resolution and independent theory checks required.")
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main(ARGS)
end
