module FiniteDifferences
export deriv4!, ko6!
"""Fourth-order accurate 5-point derivative, including 4th-order one-sided endpoints."""
function deriv4!(d::AbstractVector{Float64}, y::AbstractVector{Float64}, h::Float64)
    n = length(y)
    n >= 7 || error("At least 7 grid points are required")
    z = 1.0 / (12h)
    @inbounds begin
        d[1] = z*(-25y[1] + 48y[2] - 36y[3] + 16y[4] - 3y[5])
        d[2] = z*(-3y[1] - 10y[2] + 18y[3] - 6y[4] + y[5])
        for i in 3:n-2
            d[i] = z*(y[i-2] - 8y[i-1] + 8y[i+1] - y[i+2])
        end
        d[n-1] = z*(3y[n] + 10y[n-1] - 18y[n-2] + 6y[n-3] - y[n-4])
        d[n] = z*(25y[n] - 48y[n-1] + 36y[n-2] - 16y[n-3] + 3y[n-4])
    end
    return d
end

"""Sixth-derivative KO filter: negative semidefinite for Fourier modes."""
function ko6!(du::Matrix{Float64}, u::Matrix{Float64}, h::Float64, strength::Float64)
    strength == 0 && return
    n = size(u, 1)
    weight = strength / (64h)
    @inbounds for col in 1:size(u, 2), i in 4:n-3
        du[i,col] += weight*(u[i-3,col] - 6u[i-2,col] + 15u[i-1,col]
                            - 20u[i,col] + 15u[i+1,col] - 6u[i+2,col] + u[i+3,col])
    end
end

end # module FiniteDifferences
