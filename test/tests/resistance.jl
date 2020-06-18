#using SymEngine

# test grid for determine_impedance fnc
net = @network begin
    vs  = ac_source(V= 100,pins=3)
    r1 = impedance(z = 10, pins = 1)
    r4 = impedance(z = 10, pins = 1)
    r5 = impedance(z = 10, pins = 1)
    r2 = impedance(z = 20, pins = 1)
    r3 = impedance(z = 30, pins = 1)

    vs[1.1] ⟷ r1[1.1] ⟷ Node1
    vs[1.2] ⟷ r2[1.1] ⟷ Node2
    vs[1.3] ⟷ r3[1.1] ⟷ Node3
    vs[2.1] ⟷ vs[2.2] ⟷ vs[2.3] ⟷ gnd 
    r1[2.1] ⟷ r4[1.1] ⟷ r5[1.1] ⟷ Node4
    r4[2.1] ⟷ r5[2.1] ⟷ gnd1
    r2[2.1] ⟷ gnd1
    r3[2.1] ⟷ gnd1
end

imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1,:Node2,:Node3],
                            output_pins = Any[:gnd1,:gnd1,:gnd1], omega_range = (0.8981798683,5.798179868,500))
bode(imp, omega = omega, axis_type = :loglin)
