<#
.SYNOPSIS
    Render all currently supported 5+1 MMIEF figures from an existing run.
.DESCRIPTION
    Calls the existing Julia dashboard, observer-history, and movie scripts.
    Does not rerun the PDE or change any solver parameters.
.EXAMPLE
    .\scripts\render_all.ps1 -RunDir results/quick_fixed
.EXAMPLE
    .\scripts\render_all.ps1 -RunDir results/standard_fixed -SkipMovie
#>
param(
    [string]$RunDir = "results/quick_fixed",
    [switch]$SkipMovie
)

$ErrorActionPreference = "Stop"
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
Push-Location $projectRoot
try {
    if (-not (Get-Command julia -ErrorAction SilentlyContinue)) {
        throw "Julia not found on PATH. Start Julia from VS Code or add Julia to PATH."
    }
    if (-not (Test-Path -Path $RunDir -PathType Container)) {
        throw "Run directory '$RunDir' not found. Generate data with scripts/run_simulation.jl first."
    }
    foreach ($name in @("config.toml", "metadata.toml", "diagnostics.csv", "observers.csv", "snapshots.csv", "final_state.csv")) {
        if (-not (Test-Path -Path (Join-Path $RunDir $name) -PathType Leaf)) {
            throw "Required data file missing: $(Join-Path $RunDir $name)"
        }
    }

    Write-Host "Rendering dashboard for $RunDir"
    & julia --project=viz scripts/render_dashboard.jl $RunDir
    if ($LASTEXITCODE -ne 0) { throw "Dashboard renderer exited with code $LASTEXITCODE" }

    Write-Host "Rendering fixed-J observer histories for $RunDir"
    & julia --project=viz scripts/render_observer_history.jl $RunDir
    if ($LASTEXITCODE -ne 0) { throw "Observer-history renderer exited with code $LASTEXITCODE" }

    if (-not $SkipMovie) {
        Write-Host "Rendering scalar pulse MP4 for $RunDir"
        & julia --project=viz scripts/render_movie.jl $RunDir
        if ($LASTEXITCODE -ne 0) { throw "Movie renderer exited with code $LASTEXITCODE" }
    }

    Write-Host "Output files saved under: $RunDir"
} finally {
    Pop-Location
}
