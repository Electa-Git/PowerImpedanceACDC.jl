# include("../../src/HVDCstability.jl")
# using .HVDCstability
using SymEngine

# semiconducting layers test
net = @network begin
    ac = ac_source(V = 380, pins = 3)
    c = cable(length = 170e3, positions = [(-0.5,1.9), (0,1.9), (0.5,1.9)],
              C1 = Conductor(rₒ = 31.75e-3, ρ = 2.47e-8),
              C2 = Conductor(rᵢ = 57.55e-3, rₒ = 60.65e-3, ρ = 2.83e-8),
              I1 = Insulator(rᵢ = 31.75e-3, a = 33.75e-3, b = 56.25e-3, rₒ = 57.55e-3, ϵᵣ = 2.3),
              I2 = Insulator(rᵢ = 60.65e-3, rₒ = 65.95e-3, ϵᵣ = 2.3))
    ac[1.1] ⟷ c[1.1] ⟷ Node1
    ac[1.2] ⟷ c[1.2] ⟷ Node2
    ac[1.3] ⟷ c[1.3] ⟷ Node3
    ac[2.1] ⟷ ac[2.2] ⟷ ac[2.3] ⟷ gnd
    c[2.1] ⟷ gnd1
    c[2.2] ⟷ gnd2
    c[2.3] ⟷ gnd3
end
imp, omega = determine_impedance(net, elim_elements = [:ac], input_pins = Any[:Node1, :Node2, :Node3],
                            output_pins = Any[:gnd1, :gnd2, :gnd3], omega_range = (-3,5,1000))
bode(imp, omega = omega, axis_type = :loglog)

# cross-bonding cable test
s = symbols(:s)
net = @network begin
    ac = ac_source(V = 380, pins = 3)
    c = crossbonded_cable(C1 = cable(length = 834, positions = [(-0.5,1.9), (0,1.9), (0.5,1.9)],
              C1 = Conductor(rₒ = 31.75e-3, ρ = 2.47e-8),
              C2 = Conductor(rᵢ = 57.55e-3, rₒ = 60.65e-3, ρ = 2.83e-8),
              I1 = Insulator(rᵢ = 31.75e-3, a = 33.75e-3, b = 56.25e-3, rₒ = 57.55e-3, ϵᵣ = 2.3),
              I2 = Insulator(rᵢ = 60.65e-3, rₒ = 65.95e-3, ϵᵣ = 2.3)), nₛ = 3, mₛ = 4,
              Zᶜᵇ = 1.1)
    ac[1.1] ⟷ c[1.1] ⟷ Node1
    ac[1.2] ⟷ c[1.2] ⟷ Node2
    ac[1.3] ⟷ c[1.3] ⟷ Node3
    ac[2.1] ⟷ ac[2.2] ⟷ ac[2.3] ⟷ gnd
    c[2.1] ⟷ gnd1
    c[2.2] ⟷ gnd2
    c[2.3] ⟷ gnd3
end
imp, omega = determine_impedance(net, elim_elements = [:ac], input_pins = Any[:Node1, :Node2, :Node3],
                            output_pins = Any[:gnd1, :gnd2, :gnd3], omega_range = (-3,5,1000))
bode(imp, omega = omega, axis_type = :loglog)
