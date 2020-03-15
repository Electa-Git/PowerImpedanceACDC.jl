# include("../../src/HVDCstability.jl")
# using .HVDCstability

# 3 phase OHL case
net = @network begin
    vs = dc_source(V = 5)
    l = overhead_line(length = 200e3,
                    conductors = Conductors(organization = :flat,
                        nᵇ = 3, Rᵈᶜ = 0.063, rᶜ = 0.015, yᵇᶜ = 30,
                        Δxᵇᶜ = 10,  dˢᵃᵍ = 10),
                    groundwires = Groundwires(nᵍ    = 2, Rᵍᵈᶜ    = 0.92,
                        rᵍ      = 0.0062, Δxᵍ     = 6.5, Δyᵍ     = 7.5, dᵍˢᵃᵍ   = 10),
                    earth_parameters = (1,1,100),
                    transformation = true)
    vs[1.1] ⟷ l[1.1] ⟷ Node1
    vs[2.1] ⟷ l[1.2] ⟷ Node2
    l[2.1] ⟷ gnd
    l[2.2] ⟷ gnd1
end
imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1, :Node2],
                output_pins = Any[:gnd, :gnd1], omega_range = (-1, 6, 500))
bode(imp, omega = omega, titles = ["Z_{dd}" "Z_{dq}"; "Z_{qd}" "Z_{qq}"])

# 3 phase cable case
net = @network begin
    ac = ac_source(V = 380, pins = 2)
    c = cable(length = 170e3, positions = [(-0.5,1.9), (0,1.9), (0.5,1.9)],
              C1 = Conductor(rₒ = 31.75e-3, ρ = 2.47e-8),
              C2 = Conductor(rᵢ = 57.55e-3, rₒ = 60.65e-3, ρ = 2.83e-8),
              I1 = Insulator(rᵢ = 31.75e-3, a = 33.75e-3, b = 56.25e-3, rₒ = 57.55e-3, ϵᵣ = 2.3),
              I2 = Insulator(rᵢ = 60.65e-3, rₒ = 65.95e-3, ϵᵣ = 2.3),
              transformation = true)
    ac[1.1] ⟷ c[1.1] ⟷ Node1
    ac[1.2] ⟷ c[1.2] ⟷ Node2
    ac[2.1] ⟷ ac[2.2] ⟷ gnd
    c[2.1] ⟷ gnd1
    c[2.2] ⟷ gnd2
end
imp, omega = determine_impedance(net, elim_elements = [:ac], input_pins = Any[:Node1, :Node2],
                            output_pins = Any[:gnd1, :gnd2], omega_range = (-3,5,1000))
bode(imp, omega = omega, axis_type = :loglog)

# impedance case
using SymEngine
s = symbols(:s)
net = @network begin
    vs = dc_source(V = 5)
    z = impedance(z = 1/1e-3/s, pins = 3, transformation = true)
    vs[1.1] ⟷ z[1.1] ⟷ Node1
    vs[2.1] ⟷ z[1.2] ⟷ Node2
    z[2.1] ⟷ gnd
    z[2.2] ⟷ gnd1
end
imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1, :Node2],
                output_pins = Any[:gnd, :gnd1], omega_range = (-1, 6, 500))
bode(imp, omega = omega, titles = ["Z_{dd}" "Z_{dq}"; "Z_{qd}" "Z_{qq}"])
