using SymEngine
s = symbols("s")
net = @network begin
    vs = ac_source(V = 400, pins = 3)
    a = fautotransformer( S=600, Zₕₗ= 0.15im, Zₕₜ= 0.02im, Zₜₗ=0.02im, V₁ᵒ=400, V₂ᵒ=225)
    z = impedance(z = 25e-3*s, pins=3)

    vs[1.1] ⟷ a[1.1] ⟷ Node1
    vs[1.2] ⟷ a[1.2] ⟷ Node2
    vs[1.3] ⟷ a[1.3] ⟷ Node3
    a[2.1] ⟷ z[1.1]
    a[2.2] ⟷ z[1.2]
    a[2.3] ⟷ z[1.3]
    vs[2.1] ⟷ vs[2.2] ⟷ vs[2.3] ⟷ gnd
    z[2.1] ⟷ z[2.2] ⟷ z[2.3] ⟷ gnd

end
#imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1, :Node2, :Node3],
#                                        output_pins= Any[:gnd], omega_range = (-2, 6, 1000))
# bode(imp, omega = omega)
