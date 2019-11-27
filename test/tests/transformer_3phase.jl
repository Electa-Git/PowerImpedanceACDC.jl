
s = symbols("s")
net = @network begin
    vs = ac_source(V = 10, pins = 3)
    t = transformer(V₁ᵒ = 2.4e3, V₁ˢ = 51.87, V₂ᵒ = 240, P₁ᵒ = 171.1, P₁ˢ = 642.1,
                    I₁ᵒ = 0.48, I₁ˢ = 20.83, Cₛ = 12e-6, Cₜ = 7e-6, pins = 3)
    z = impedance(pins = 3, z = 25e-3*s)

    vs[1.1] ⟷ t[1.1]
    vs[1.2] ⟷ t[1.2]
    vs[1.3] ⟷ t[1.3]
    t[2.1] ⟷ z[1.1]
    t[2.2] ⟷ z[1.2]
    t[2.3] ⟷ z[1.3]
    vs[2.1] ⟷ vs[2.2] ⟷ vs[2.3] ⟷ z[2.1] ⟷ z[2.2] ⟷ z[2.3] ⟷ gnd
end
