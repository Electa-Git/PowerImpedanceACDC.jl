# For debugging 
# include("../../src/HVDCstability.jl")
# using .HVDCstability

# For better match with sweep: In MMC.jl, VGd = Vm and VGq = 0

using DelimitedFiles
using SymEngine
using Plots
using LinearAlgebra

s = symbols("s")

#Reactive power compensation


net = @network begin

    Zg1 = impedance(z = 1 , pins = 3, transformation = true)
    Zg2 = impedance(z = 2 , pins = 3, transformation = true)
    # Zg3 = impedance(z = 4 , pins = 3, transformation = true)
    # Z_shunt = impedance(z = 5 , pins = 3, transformation = true)
    # Z_shunt2 = impedance(z = 10 , pins = 3, transformation = true)

    # Since the impedances are all passive, the voltage magnitude does not matter.
    g1 = ac_source(V = 5, pins = 3, transformation = true)
    g2 = ac_source(V = 5, pins = 3, transformation = true)

    # ABC
    # g1[1.1] ⟷ Zg1[1.1] == Zg3[1.1] == Z_shunt[1.1] == Node1A
    # g1[1.2] ⟷ Zg1[1.2] == Zg3[1.2] == Z_shunt[1.2] == Node1B
    # g1[1.3] ⟷ Zg1[1.3] == Zg3[1.3] == Z_shunt[1.3] == Node1C

    # Zg1[2.1] ⟷ Zg2[1.1] == Z_shunt2[1.1] == Node2A
    # Zg1[2.2] ⟷ Zg2[1.2] == Z_shunt2[1.2] == Node2B
    # Zg1[2.3] ⟷ Zg2[1.3] == Z_shunt2[1.3] == Node2C

    # Zg2[2.1] ⟷ g2[1.1] == Zg3[2.1] == Node3A
    # Zg2[2.2] ⟷ g2[1.2] == Zg3[2.2] == Node3B
    # Zg2[2.3] ⟷ g2[1.3] == Zg3[2.3] == Node3C

    # Z_shunt[2.1] == gnda
    # Z_shunt[2.2] == gndb
    # Z_shunt[2.3] == gndc

    # Z_shunt2[2.1] == gnda
    # Z_shunt2[2.2] == gndb
    # Z_shunt2[2.3] == gndc

    # g1[2.1] ⟷ gnda
    # g1[2.2] ⟷ gndb
    # g1[2.3] ⟷ gndc

    # g2[2.1] ⟷ gnda
    # g2[2.2] ⟷ gndb
    # g2[2.3] ⟷ gndc
    # # Single phase
    # g1[1.1] ⟷ Zg1[1.1] == Zg3[1.1] == Z_shunt[1.1] == Node1A
    # Zg1[2.1] ⟷ Zg2[1.1] == Z_shunt2[1.1] == Node2A
    # Zg2[2.1] ⟷ g2[1.1] == Zg3[2.1] == Node3A
    # Z_shunt[2.1] == gnda
    # Z_shunt2[2.1] == gnda
    # g1[2.1] ⟷ gnda
    # g2[2.1] ⟷ gnda

    # dq
    # g1[1.1] ⟷ Zg1[1.1] == Zg3[1.1] == Z_shunt[1.1] == Node1D
    # g1[1.2] ⟷ Zg1[1.2] == Zg3[1.2] == Z_shunt[1.2] == Node1Q

    # Zg1[2.1] ⟷ Zg2[1.1] == Z_shunt2[1.1] == Node2D
    # Zg1[2.2] ⟷ Zg2[1.2] == Z_shunt2[1.2] == Node2Q

    # Zg2[2.1] ⟷ g2[1.1] == Zg3[2.1] == Node3D
    # Zg2[2.2] ⟷ g2[1.2] == Zg3[2.2] == Node3Q

    # Z_shunt[2.1] == gndd
    # Z_shunt[2.2] == gndq

    # Z_shunt2[2.1] == gndd
    # Z_shunt2[2.2] == gndq

    # g1[2.1] ⟷ gndd
    # g1[2.2] ⟷ gndq

    # g2[2.1] ⟷ gndd
    # g2[2.2] ⟷ gndq

    # debugging
    g1[1.1] ⟷ Zg1[1.1] == Node1D
    g1[1.2] ⟷ Zg1[1.2] == Node1Q

    Zg1[2.1] ⟷ Zg2[1.1] == Node2D
    Zg1[2.2] ⟷ Zg2[1.2] == Node2Q

    Zg2[2.1] ⟷ g2[1.1] == Node3D
    Zg2[2.2] ⟷ g2[1.2] == Node3Q

    g1[2.1] ⟷ gndd
    g1[2.2] ⟷ gndq

    g2[2.1] ⟷ gndd
    g2[2.2] ⟷ gndq

end

# Y_matrix_AC = HVDCstability.make_y_matrix(net, elim_elements = [:Z_shunt,:Z_shunt2], input_pins = Any[:Node1A,:Node1B,:Node1C,:Node2A,:Node2B,:Node2C,:Node3A,:Node3B,:Node3C], omega_range = (-3, 4, 1000));
# Y_matrix_AC = HVDCstability.make_y_matrix(net,input_pins = Any[:Node1A,:Node2A,:Node3A], omega_range = (-3, 4, 1000));
# Y_matrix_AC = HVDCstability.make_y_matrix(net, elim_elements = [:Z_shunt,:Z_shunt2], input_pins = Any[:Node1D,:Node1Q,:Node2D,:Node2Q,:Node3D,:Node3Q], omega_range = (0, 1, 10));
Y_matrix_AC = HVDCstability.make_y_matrix(net, elim_elements = [:Z_shunt,:Z_shunt2,:Zg2,:Zg3], input_pins = Any[:Node1D,:Node1Q,:Node2D,:Node2Q], omega_range = (0, 1, 10));

Y_BUS_ABCD = []
omega= 2*pi* 10 .^range(0, 1, length= 10)
ω₀ = 100*pi 
    for i in eachindex(omega)

        ABCD1 = HVDCstability.transformation_dq(eval_abcd(net.elements[:Zg1].element_value, 1im*omega[i] + 1im*ω₀),
        eval_abcd(net.elements[:Zg1].element_value, 1im*omega[i] - 1im*ω₀))
        ABCD2 = HVDCstability.transformation_dq(eval_abcd(net.elements[:Zg2].element_value, 1im*omega[i] + 1im*ω₀),
        eval_abcd(net.elements[:Zg2].element_value, 1im*omega[i] - 1im*ω₀))

        ABCDeq = Matrix{Complex}(ABCD1 * ABCD2)
        # ABCDeq = ABCD1
        push!(Y_BUS_ABCD, HVDCstability.abcd_to_y(ABCDeq))
    end
