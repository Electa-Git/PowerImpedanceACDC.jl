isdefined(Base, :__precompile__) && __precompile__()

module HVDCstability

    # default struct values
    using ForwardDiff: isodd
    using Parameters, DataStructures

    # files manipulation
    using FileIO, DelimitedFiles

    # plotting
    using PGFPlotsX
    using Compat, Dates # compatibility across Julia

    # symbolic and numerical calculations
    using SymEngine
    using LinearAlgebra
    using NLsolve, ForwardDiff  # solve diffs and nonlinear equations

    # Power flow
    using PowerModels, PowerModelsACDC
    using Ipopt
    using JuMP

    # New Nyquist plots from Thomas
    using Plots
    using Combinatorics

    # file includes
    include("Network/globals.jl")

    include("Network/compat.jl")
    include("Network/Components/AbstractElement.jl")
    #include("GUI/Interactions.jl")
    include("Network/Solvers/solvers.jl")

end
