module PowerImpedanceACDC
"""
Copyright (C) 2024  Jan Kircheis  
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 
"""
    # Default struct values
    using Parameters
    using DataStructures

    # Time delays
    using ControlSystemsBase, RobustAndOptimalControl

    # Plotting
    using Plots

    # Symbolic and numerical calculations
    using SymEngine
    using LinearAlgebra
    using NLsolve, ForwardDiff,NonlinearSolve, SteadyStateDiffEq  # solve diffs and nonlinear equations
    
    # Power flow
    using PowerModels, PowerModelsACDC
    using Ipopt
    using JuMP

    # For tools
    using Peaks

    # Miscellaneous
    using Munkres
    using Markdown

    # Including all components    
    include("Network/Components/AbstractElement.jl")

    # Impedance
    include("Network/Components/impedance/impedance.jl")

    # Transfromer
    include("Network/Components/transformer/transformer.jl")

    # Shunt reactor
    include("Network/Components/shunt_reactor/shunt_reactor.jl")

    # Cables and transmission lines
    include("Network/Components/transmission_line/transmission_line.jl")
    include("Network/Components/transmission_line/cable.jl")
    include("Network/Components/transmission_line/overhead_line.jl")
    include("Network/Components/transmission_line/mixed_OHL_cable.jl")

    # Grid or source
    include("Network/Components/source/source.jl")
    include("Network/Components/source/dc_source.jl")
    include("Network/Components/source/ac_source.jl")

    # Converter
    include("Network/Components/converter/converter.jl")
    include("Network/Components/converter/controller.jl")
    include("Network/Components/converter/MMC.jl")
    include("Network/Components/converter/TLC.jl")

    # Machines
    include("Network/Components/machine/machine.jl")
    include("Network/Components/machine/SynchronousMachine.jl")


    # Including network
    include("Network/Network.jl")

    # New power flow
    include("Network/power_flow.jl")
    
    # Including network solvers
    include("Network/Solvers/make_abcd.jl")
    include("Network/Solvers/make_y.jl")
    include("Network/Solvers/make_z.jl")
    include("Network/Solvers/determine_impedance.jl")
    include("Network/Solvers/make_y_matrix.jl")
    include("Network/Solvers/stability.jl")

    # Including tools
    include("Tools/abcd_parameters.jl")
    include("Tools/kron.jl")
    include("Tools/nyquistplot.jl")
    include("Tools/small_gain.jl")
    include("Tools/stabilitymargin.jl")
    include("Tools/EVD.jl")
    include("Tools/bodeplot.jl")
    include("Tools/passivity.jl")

    
    
end
