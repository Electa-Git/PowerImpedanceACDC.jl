
s = symbols("s")
net = @network begin
    vs = dc_source(voltage = 5)
    z1 = impedance(z = s+2, pins = 1)
    z2 = impedance(z = s, pins = 1)
    z3 = impedance(z = s, pins = 1)

    vs[1.1] ⟷ z1[1.1] ⟷ Node1
    z1[2.1] ⟷ z2[1.1] ⟷ z3[1.1]
    vs[2.1] ⟷ z2[2.1] ⟷ z3[2.1] ⟷ gnd
end
imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1],
                                        output_pins= Any[:gnd])
# bode(imp, omega = omega)
