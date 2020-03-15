# include("../../src/HVDCstability.jl")
# using .HVDCstability

using SymEngine

s = symbols(:s)
ω = 2π*50

net = @network begin
    ac = ac_source(V = 230)
    x_sys = impedance(pins = 1, z = 0.04 + s*0.3/ω)
    x_l1 = impedance(pins = 1, z = 0.835 + 4*s/ω)
    x_l2 = impedance(pins = 1, z = 0.835 + 4*s/ω)
    c_1 = impedance(pins = 1, z = 1/(0.0013 * s / ω))
    c_2 = impedance(pins = 1, z = 1/(0.0013 * s / ω))

    ac[1.1] ⟷ x_sys[1.1]
    x_sys[2.1] ⟷ c_1[1.1] ⟷ x_l1[1.1] ⟷ bus1
    x_l1[2.1] ⟷ x_l2[1.1] ⟷ bus2
    x_l2[2.1] ⟷ c_2[1.1] ⟷ bus3
    ac[2.1] ⟷ c_1[2.1] ⟷ c_2[2.1] ⟷ gnd
end

imp, omega = determine_impedance(net, input_pins = Any[:bus1],
    output_pins = Any[:gnd], elim_elements = Symbol[], parameters_type = :Y,
    omega_range = (0, 4.5, 5000))

bode(imp, omega = omega / ω, axis_type = :linlin)
