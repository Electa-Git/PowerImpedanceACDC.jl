include("../../src/HVDCstability.jl")
using .HVDCstability

@time net = @network begin
    gen1 = ac_source(pins = 3, P_min = 50, P = 1200, P_max = 1500, Q = 0, Q_max = 500, Q_min = -500,
                    V = 380, transformation = true)
    gen2 = ac_source(pins = 3, P_min = -1500, P = -600, P_max = -50, Q = 0, Q_max = 500, Q_min = -500,
                    V = 380, transformation = true)
    gen3 = ac_source(pins = 3, P_min = -1500, P = -200, P_max = -50, Q = 0, Q_max = 500, Q_min = -500,
                    V = 380, transformation = true)
    gen4 = ac_source(pins = 3, P_min = -1500, P = -400, P_max = -50, Q = 0, Q_max = 500, Q_min = -500,
                    V = 380, transformation = true)

    tl1 = overhead_line(length = 20e3,
        conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
        Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        earth_parameters = (1,1,100), transformation = true)
    tl2 = overhead_line(length = 200e3,
        conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
        Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        earth_parameters = (1,1,100), transformation = true)
    tl3 = overhead_line(length = 200e3,
        conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
        Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        earth_parameters = (1,1,100), transformation = true)
    tl4 = overhead_line(length = 200e3,
        conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
        Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        earth_parameters = (1,1,100), transformation = true)


    dc_line1 = cable(length = 100e3, positions = [(-0.5,1), (0.5,1)],
          C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
          C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
          C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
          I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
          I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
          I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)
    dc_line2 = cable(length = 150e3, positions = [(-0.5,1), (0.5,1)],
            C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
            C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
            C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
            I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
            I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
            I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)
    dc_line3 = cable(length = 200e3, positions = [(-0.5,1), (0.5,1)],
            C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
              C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
              C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
              I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
              I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
              I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)
    dc_line4 = cable(length = 200e3, positions = [(-0.5,1), (0.5,1)],
            C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
            C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
            C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
            I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
            I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
            I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)
    dc_line5 = cable(length = 100e3, positions = [(-0.5,1), (0.5,1)],
          C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
          C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
          C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
          I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
          I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
          I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)

    c1 = mmc(Vᵈᶜ = 640, Vₘ = 320,
            P_max = -500, P_min = -1500, P = -1200, Q = 0, Q_max = 500, Q_min = -500, P_dc = 1000,
            Rᵣ = 0.3429, Lᵣ = 62.9e-3, Rₐᵣₘ = 0.6017, Lₐᵣₘ = 30.6e-3, N = 200, Cₐᵣₘ = 4.23e-3,
            occ = PI_control(ζ = 0.7, bandwidth = 1000),
            ccc = PI_control(ζ = 0.7, bandwidth = 300),
            dc = PI_control(Kₚ = 0.01, Kᵢ = 2),
            energy = PI_control(Kₚ = 120, Kᵢ = 400),
            zcc = PI_control(ζ = 0.7, bandwidth = 300)
            )
    c2 = mmc(Vᵈᶜ = 640, Vₘ = 320,
            P_max = 1500, P_min = 300, P = 600, Q = 0, Q_max = 300, Q_min = -300, P_dc = -1000,
            Rᵣ = 0.3429, Lᵣ = 62.9e-3, Rₐᵣₘ = 0.6017, Lₐᵣₘ = 30.6e-3, N = 200, Cₐᵣₘ = 4.23e-3,
            occ = PI_control(ζ = 0.7, bandwidth = 1000),
            ccc = PI_control(ζ = 0.7, bandwidth = 300),
            power = PI_control(Kₚ = 2.0020e-07, Kᵢ = 1.0010e-04),
            energy = PI_control(Kₚ = 120, Kᵢ = 400),
            zcc = PI_control(ζ = 0.7, bandwidth = 300)
            )
    c3 = mmc(Vᵈᶜ = 640, Vₘ = 320,
            P_max = 1500, P_min = 300, P = 200, Q = 0, Q_max = 300, Q_min = -300, P_dc = -1000,
            Rᵣ = 0.3429, Lᵣ = 62.9e-3, Rₐᵣₘ = 0.6017, Lₐᵣₘ = 30.6e-3, N = 200, Cₐᵣₘ = 4.23e-3,
            occ = PI_control(ζ = 0.7, bandwidth = 1000),
            ccc = PI_control(ζ = 0.7, bandwidth = 300),
            power = PI_control(Kₚ = 2.0020e-07, Kᵢ = 1.0010e-04),
            energy = PI_control(Kₚ = 120, Kᵢ = 400),
            zcc = PI_control(ζ = 0.7, bandwidth = 300)
            )
    c4 = mmc(Vᵈᶜ = 640, Vₘ = 320,
            P_max = 1500, P_min = 300, P = 400, Q = 0, Q_max = 300, Q_min = -300, P_dc = -1000,
            Rᵣ = 0.3429, Lᵣ = 62.9e-3, Rₐᵣₘ = 0.6017, Lₐᵣₘ = 30.6e-3, N = 200, Cₐᵣₘ = 4.23e-3,
            occ = PI_control(ζ = 0.7, bandwidth = 1000),
            ccc = PI_control(ζ = 0.7, bandwidth = 300),
            power = PI_control(Kₚ = 2.0020e-07, Kᵢ = 1.0010e-04),
            energy = PI_control(Kₚ = 120, Kᵢ = 400),
            zcc = PI_control(ζ = 0.7, bandwidth = 300)
            )

    # connections
    gen1[2.1] ⟷ gen1[2.2] ⟷ gnd1
    gen1[1.1] ⟷ tl1[1.1]
    gen1[1.2] ⟷ tl1[1.2]
    tl1[2.1] ⟷ c1[2.1]
    tl1[2.2] ⟷ c1[2.2]

    gen2[1.1] ⟷ tl2[1.1]
    gen2[1.2] ⟷ tl2[1.2]
    tl2[2.1] ⟷ c2[2.1]
    tl2[2.2] ⟷ c2[2.2]
    gen2[2.1] ⟷ gen2[2.2] ⟷ gnd2

    gen3[1.1] ⟷ tl3[1.1]
    gen3[1.2] ⟷ tl3[1.2]
    tl3[2.1] ⟷ c3[2.1]
    tl3[2.2] ⟷ c3[2.2]
    gen3[2.1] ⟷ gen3[2.2] ⟷ gnd3

    gen4[1.1] ⟷ tl4[1.1]
    gen4[1.2] ⟷ tl4[1.2]
    tl4[2.1] ⟷ c4[2.1]
    tl4[2.2] ⟷ c4[2.2]
    gen4[2.1] ⟷ gen4[2.2] ⟷ gnd4

    c1[1.1] ⟷ dc_line1[1.1]
    c2[1.1] ⟷ dc_line1[2.1]

    c1[1.1] ⟷ dc_line2[1.1]
    c3[1.1] ⟷ dc_line2[2.1]

    c1[1.1] ⟷ dc_line3[1.1]
    c4[1.1] ⟷ dc_line3[2.1]

    c2[1.1] ⟷ dc_line4[1.1]
    c3[1.1] ⟷ dc_line4[2.1]

    c3[1.1] ⟷ dc_line5[1.1]
    c4[1.1] ⟷ dc_line5[2.1]
end

@time imp, omega = check_stability(net, net.elements[:c1], direction = :dc)
bode(imp, omega = omega, titles = ["Z_{MMC1}" "Z_{eq}" "Y_{MMC1} Z_{eq}"])
