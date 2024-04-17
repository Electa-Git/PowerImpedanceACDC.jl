export Controller, PI_control, VSE, VI

abstract type Controller end


@with_kw mutable struct PI_control <: Controller
    bandwidth :: Union{Int, Float64}  = 0
    ζ :: Union{Int, Float64}  = 0
    Kₚ :: Union{Int, Float64}  = 0                     # proportional gain
    Kᵢ :: Union{Int, Float64}  = 0                     # integral gain
    ref :: Array{Union{Int, Float64}}  = [0]           # reference value
    n_f :: Int = 0                                     # Butterworth Filter order
    ω_f  :: Union{Int, Float64}  = 0                   # Filter cutoff frequency
end

@with_kw mutable struct VSE <: Controller
    H :: Union{Int, Float64}  = 0                       # Virtual Inertia 
    K_d :: Union{Int, Float64}  = 0                     # Damping coefficient
    K_ω    :: Union{Int, Float64}  = 0                  # Droop coefficient
    ref :: Array{Union{Int, Float64}}  = [0]            # reference value
    n_f :: Int = 0                                      # Butterworth Filter order
    ω_f  :: Union{Int, Float64}  = 0                    # Butterworth Filter cutoff frequency
end

@with_kw mutable struct VI <: Controller                 
    Rᵥ :: Union{Int, Float64}  = 0                       # Virtual Resistance
    Lᵥ :: Union{Int, Float64}  = 0                       # Virtual inductance
    n_f :: Int  = 0                                      # Butterworth Filter order
    ω_f  :: Union{Int, Float64}  = 0                     # Butterworth Filter cutoff frequency
end