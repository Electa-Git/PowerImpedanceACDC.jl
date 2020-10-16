
s = symbols("s")
net = @network begin
    vs = ac_source(V = 10, pins = 3)
    t = transformer(V₁ᵒ = 2.4e3, V₁ˢ = 51.87, V₂ᵒ = 240, P₁ᵒ = 171.1, P₁ˢ = 642.1,
                    I₁ᵒ = 0.48, I₁ˢ = 20.83, Cₛ = 12e-6, Cₜ = 7e-6, pins = 3)
    #z = impedance(pins = 3, z = 25e-3*s)

    vs[1.1] ⟷ t[1.1] ⟷ Node1
    vs[1.2] ⟷ t[1.2] ⟷ Node2
    vs[1.3] ⟷ t[1.3] ⟷ Node3
    #t[2.1] ⟷ z[1.1]
    #t[2.2] ⟷ z[1.2]
    #t[2.3] ⟷ z[1.3]
    vs[2.1] ⟷ vs[2.2] ⟷ vs[2.3] ⟷gnd
    t[2.1] ⟷ gnd1
    t[2.2] ⟷ gnd2
    t[2.3] ⟷ gnd3
end

imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1],
                                        output_pins= Any[:gnd1], omega_range = (-2, 6, 1000))
 bode(imp, omega = omega)
