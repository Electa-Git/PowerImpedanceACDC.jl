export Controller, PI_control, VSE, CCQSEM

abstract type Controller end


@with_kw mutable struct PI_control <: Controller
    bandwidth :: Union{Int, Float64}  = 0
    ζ :: Union{Int, Float64}  = 0
    Kₚ :: Union{Int, Float64}  = 0                     # proportional gain [pu]
    Kᵢ :: Union{Int, Float64}  = 0                     # integral gain [pu/s]
    ref :: Array{Union{Int, Float64}}  = [0]           # reference value
    n_f :: Int = 0                                     # Butterworth Filter order [-]
    ω_f  :: Union{Int, Float64}  = 0                   # Filter cutoff frequency [rad/s]
end

@with_kw mutable struct VSE <: Controller
    H :: Union{Int, Float64}  = 0                       # Virtual Inertia [s]
    K_d :: Union{Int, Float64}  = 0                     # Damping coefficient [-]
    K_ω    :: Union{Int, Float64}  = 0                  # Droop coefficient [-]
    ref :: Array{Union{Int, Float64}}  = [0]            # Active power reference [MW]
    ref_ω :: Union{Int, Float64}  = 0                   # Angular frequency reference [pu]
    n_f :: Int = 0                                      # Butterworth Filter order [-]
    ω_f  :: Union{Int, Float64}  = 0                    # Butterworth Filter cutoff frequency [rad/s]
end

@with_kw mutable struct CCQSEM <: Controller                 
    Rᵥ :: Union{Int, Float64}  = 0                       # Virtual Resistance [pu]
    Lᵥ :: Union{Int, Float64}  = 0                       # Virtual inductance [pu]
    ref_vd :: Union{Int, Float64}  = 0                   # Vd voltage reference [pu]
    ref_vq :: Union{Int, Float64}  = 0                   # Vq voltage reference [pu]
    n_f :: Int  = 0                                      # Butterworth Filter order [-]
    ω_f  :: Union{Int, Float64}  = 0                     # Butterworth Filter cutoff frequency [rad/s]


end