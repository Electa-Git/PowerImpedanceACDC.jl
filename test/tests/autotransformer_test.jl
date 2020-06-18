using SymEngine
s = symbols("s")
net = @network begin
    vs = ac_source(V = 10, pins = 3)
    at = autotransformer(Xₕₗ = 0.15, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)
    z = impedance(pins = 3, z = 25e-3*s)

    vs[1.1] ⟷ at[1.1]
    vs[1.2] ⟷ at[1.2]
    vs[1.3] ⟷ at[1.3]
    at[2.1] ⟷ z[1.1]
    at[2.2] ⟷ z[1.2]
    at[2.3] ⟷ z[1.3]
    vs[2.1] ⟷ vs[2.2] ⟷ vs[2.3] ⟷ z[2.1] ⟷ z[2.2] ⟷ z[2.3] ⟷ gnd
end
