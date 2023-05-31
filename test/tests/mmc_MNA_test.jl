using HVDCstability
net = @network begin
        c = mmc(energy = PI_control(Kₚ = 120, Kᵢ = 400),
                occ = PI_control(ζ = 0.7, bandwidth = 1000),
                ccc = PI_control(ζ = 0.7, bandwidth = 300),
                zcc = PI_control(ζ = 0.7, bandwidth = 300),
                power = PI_control(Kₚ = 2.0020e-07, Kᵢ = 1.0010e-04)#,
                )
        dc = dc_source(V = 640)
        ac = ac_source(V = 320, pins = 3, transformation = true)

        dc[1.1] ⟷ c[1.1] ⟷ Node1
        dc[2.1] ⟷ gnd
        c[2.1] ⟷ ac[1.1] ⟷ Noded
        c[2.2] ⟷ ac[1.2] ⟷ Nodeq
        ac[2.1] ⟷ ac[2.2] ⟷ gnd
        # ac[2.1] ⟷ gndd
        # ac[2.2] ⟷ gndq
end

# imp_dc, omega_dc = determine_impedance(net, elim_elements = [:ac], input_pins = Any[:Noded,:Nodeq],
#         output_pins = Any[:gnd], omega_range = (0, 4, 1000))
imp_dc, omega_dc = determine_impedance(net, elim_elements = [:dc], input_pins = Any[:Node1],
        output_pins = Any[:Noded], omega_range = (0, 4, 1000))
# bode(imp_dc, omega = omega_dc)
