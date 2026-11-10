# No plotting dependency. Fits slopes only within user-selected windows.
# Example: julia --project=. scripts/analyze_tails.jl results/standard_fixed 180 300
using MMIEFD6
using Printf
function main(args::Vector{String}=ARGS)
    length(args)==3 || error("Usage: analyze_tails.jl RUN_DIR TLO THI (coordinate time units)")
    run=load_run(args[1]);tlo=parse(Float64,args[2]);thi=parse(Float64,args[3])
    tlo<thi || error("TLO must be less than THI")
    T=run.times .- run.shifts
    @printf("5+1 exploratory tail slopes; T in length units; interval [%.2f,%.2f]\n",tlo,thi)
    @printf("Reference fixed-background l=0: p=10; dispersive small-data nonlinear: p=6 (NOT automatically valid with a black hole).\n")
    for (J,phi) in zip(run.settings.observers,run.signals)
        fit=fit_tail(T,phi;tlo=tlo,thi=thi)
        if fit===nothing
            @printf("J=%.2f | insufficient samples\n",J)
        else
            @printf("J=%.2f | p=%.5g | log RMS=%.4g | N=%d\n",J,fit.power,fit.log_rms,fit.n)
        end
    end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main(ARGS)
end
