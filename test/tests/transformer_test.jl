#circuit containing a single phase transformer and an inductive load
s = symbols("s")
net = @network begin
    vs = ac_source(V = 50, pins = 1) #single phase ac source
    t = transformer(V₁ᵒ = 2.4e3, V₁ˢ = 51.87, V₂ᵒ = 240, P₁ᵒ = 171.1, P₁ˢ = 642.1,
                    I₁ᵒ = 0.48, I₁ˢ = 20.83, Cₛ = 12e-6, Cₜ = 7e-6) #single phase transformer
    z = impedance(pins = 1, z = 25e-3*s) #inductive load

    vs[1.1] ⟷ t[1.1] ⟷ Node1
    t[2.1] ⟷ z[1.1]
    vs[2.1] ⟷ z[2.1] ⟷ gnd
end
imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1],
                                        output_pins= Any[:gnd], omega_range = (-2, 6, 1000))
 bode(imp, omega = omega)
