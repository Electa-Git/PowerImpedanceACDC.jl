s = symbols("s")
net = @network begin
    vs = ac_source(V = 380, pins = 3)
    at = autotransformer(S=600, Xₕₗ = 0.15, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)
    #z = impedance(pins = 3, z = 25e-3*s)

    vs[1.1] ⟷ at[1.1] ⟷ Node1
    vs[1.2] ⟷ at[1.2] ⟷ Node2
    vs[1.3] ⟷ at[1.3] ⟷ Node3
    #t[2.1] ⟷ z[1.1]
    #t[2.2] ⟷ z[1.2]
    #t[2.3] ⟷ z[1.3]
    vs[2.1] ⟷ vs[2.2] ⟷ vs[2.3] ⟷ at[2.1] ⟷ at[2.2] ⟷ at[2.3] ⟷ gnd
end

imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1, :Node2, :Node3],
                                        output_pins= Any[:gnd, :gnd, :gnd], omega_range = (-2, 6, 1000))
 bode(imp, omega = omega)
