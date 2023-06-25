using HVDCstability
net = @network begin
        c1 = mmc(energy = PI_control(Kₚ = 120, Kᵢ = 400),
                occ = PI_control(ζ = 0.7, bandwidth = 1000),
                ccc = PI_control(ζ = 0.7, bandwidth = 300),
                zcc = PI_control(ζ = 0.7, bandwidth = 300),
                dc = PI_control(Kₚ = 0.1, Kᵢ = 2),
                q = PI_control(Kₚ = 2.0020e-07, Kᵢ = 1.0010e-04)
                )
        c2 = mmc(energy = PI_control(Kₚ = 120, Kᵢ = 400),
                occ = PI_control(ζ = 0.7, bandwidth = 1000),
                ccc = PI_control(ζ = 0.7, bandwidth = 300),
                zcc = PI_control(ζ = 0.7, bandwidth = 300),
                p = PI_control(Kₚ = 2.0020e-07, Kᵢ = 1.0010e-04),
                q = PI_control(Kₚ = 2.0020e-07, Kᵢ = 1.0010e-04)
                )
        g1 = ac_source(V = 220, pins = 3, transformation = true)
        g4 = ac_source(V = 220, pins = 3, transformation = true)

        dc1 = impedance(z = 1, pins = 1)

        g1[2.1] == gndd
        g1[2.2] == gndq

        g1[1.1] == c1[2.1] == Node1d
        g1[1.2] == c1[2.2] == Node1q

        c1[1.1] == dc1[1.1] == Node2
        c2[1.1] == dc1[2.1] == Node3

        c2[2.1] == g4[1.1] == Node4d
        c2[2.2] == g4[1.2] == Node4q

        g4[2.1] == gndd
        g4[2.2] == gndq


end

imp_dc, omega_dc = determine_impedance(net, elim_elements = [:g4], input_pins = Any[:Node4d,:Node4q],
        output_pins = Any[:gndd,:gndq], omega_range = (0, 4, 1000))

bode(imp_dc, omega = omega_dc)
