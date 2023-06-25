using SymEngine,HVDCstability
s = symbols("s")
net = @network begin
    # vs = dc_source(V = 5, pins = 3, transformation = true)
    # z1 = impedance(z = s+2, pins = 3, transformation = true)
    # z2 = impedance(z = s, pins = 3, transformation = true)
    # z3 = impedance(z = s, pins = 3, transformation = true)

    # vs[1.1] ⟷ z1[1.1] ⟷ Node1d
    # vs[1.2] ⟷ z1[1.2] ⟷ Node1q
    # z1[2.1] ⟷ z2[1.1] ⟷ z3[1.1] ⟷ Node2d
    # z1[2.2] ⟷ z2[1.2] ⟷ z3[1.2] ⟷ Node2q
    # vs[2.1] ⟷ z2[2.1] ⟷ z3[2.1] ⟷ gndd
    # vs[2.2] ⟷ z2[2.2] ⟷ z3[2.2] ⟷ gndq
    vs = dc_source(V = 5, pins = 3, transformation = true)
    z1 = impedance(z = s+2, pins = 3, transformation = true)

    vs[1.1] ⟷ z1[1.1] ⟷ Node1d
    vs[1.2] ⟷ z1[1.2] ⟷ Node1q
    vs[2.1] ⟷ gndd
    vs[2.2] ⟷ gndq
    z1[2.1] ⟷ gndd
    z1[2.2] ⟷ gndq
end
imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1d,:Node1q],
                output_pins= Any[:gndd,:gndq], omega_range = (-2, 3, 1000),
                parameters_type = :ABCD)
bode(imp, omega = omega)
