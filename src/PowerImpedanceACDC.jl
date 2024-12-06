isdefined(Base, :__precompile__) && __precompile__()

module PowerImpedanceACDC

    # default struct values
    using Parameters
    using DataStructures

    # files manipulation
    using FileIO, DelimitedFiles

    # time delays
    using ControlSystemsBase, RobustAndOptimalControl

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
    using Munkres

    # file includes
    include("Network/globals.jl")
    include("Network/compat.jl")
    include("Network/Components/AbstractElement.jl")
    include("Network/Solvers/solvers.jl")

end
