export Controller, PI_control

abstract type Controller end

@with_kw mutable struct PI_control <: Controller
    bandwidth :: Union{Int, Float64}  = 0
    ζ :: Union{Int, Float64}  = 0
    Kₚ :: Union{Int, Float64}  = 0                     # proportional gain
    Kᵢ :: Union{Int, Float64}  = 0                     # integral gain
    ref :: Array{Union{Int, Float64}}  = [0]           # reference value
end
