# using HVDCstability

net = @network begin
        c = mmc(energy = PI_control(Kₚ = 120, Kᵢ = 400),
                occ = PI_control(ζ = 0.7, bandwidth = 1000),
                ccc = PI_control(ζ = 0.7, bandwidth = 300),
                zcc = PI_control(ζ = 0.7, bandwidth = 300),
                power = PI_control(Kₚ = 2.0020e-07, Kᵢ = 1.0010e-04),
                pll = PI_control(Kₚ = 2e-3, Kᵢ = 2)
                )
        dc = dc_source(voltage = 640e3)

        dc[1.1] ⟷ c[1.1] ⟷ Node1
        dc[2.1] ⟷ c[1.2] ⟷ Node2
        c[2.1] ⟷ gnd
        c[2.2] ⟷ gnd1
end

imp_dc, omega_dc = determine_impedance(net, elim_elements = [:dc], input_pins = Any[:Node1,:Node2],
        output_pins = Any[:gnd,:gnd1], omega_range = (0, 4, 1000))
bode(imp_dc, omega = omega_dc)

imp_ac, omega_ac = determine_impedance(net, elim_elements = [:dc,:ac], input_pins = Any[:gnd, :gnd1],
        output_pins = Any[:Node1,:Node2], omega_range = (0, 4, 1000))
bode(imp_ac, omega = omega_ac)
