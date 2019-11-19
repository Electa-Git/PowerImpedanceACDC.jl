isdefined(Base, :__precompile__) && __precompile__()

module HVDCstability

    # default struct values
    using Parameters, DataStructures

    # files manipulation
    using FileIO, DelimitedFiles

    # plotting
    using Plots, LaTeXStrings
    using Compat, Compat.Dates  # compatibility across Julia

    # symbolic and numerical calculations
    using SymEngine
    using LinearAlgebra
    using NLsolve, ForwardDiff  # solve diffs and nonlinear equations

    # Power flow
    using PowerModels, PowerModelsACDC
    using Ipopt
    using JuMP

    # file includes
    include("Network/globals.jl")

    include("Network/compat.jl")
    include("Network/Components/AbstractElement.jl")
    #include("GUI/Interactions.jl")
    include("Network/Impedance_estimation/determine_impedance.jl")
    include("Network/Impedance_estimation/stability.jl")
end
