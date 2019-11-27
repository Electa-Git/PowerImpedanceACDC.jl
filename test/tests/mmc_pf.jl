include("../../src/HVDCstability.jl")
using .HVDCstability

net = @network begin
    gen1 = ac_source(pins = 3, P_max = 250, P_min = 100, P = 100, Q = 0, Q_max = 500, Q_min = -500,
                    V = 320)
    gen2 = ac_source(pins = 3, P_max = 150, P_min = -100, P = -100, Q = 0, Q_max = 300, Q_min = -300,
                    V = 320)

    tl1 = overhead_line(length = 200e3,
        conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                        Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        earth_parameters = (1,1,100))
    tl2 = overhead_line(length = 200e3,
        conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                        Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        earth_parameters = (1,1,100))

    dc_line = overhead_line(length = 227e3, conductors = Conductors(nᵇ = 2, nˢᵇ = 2, organization = :flat,
    Rᵈᶜ = 0.06266, rᶜ = 0.01436, yᵇᶜ = 27.5, Δxᵇᶜ = 11.8, dˢᵇ = 0.4572, dˢᵃᵍ = 10),
    earth_parameters = (1,1,100),
    groundwires = Groundwires(nᵍ = 2, Δxᵍ = 6.5, Δyᵍ = 7.5, Rᵍᵈᶜ = 0.9196, rᵍ = 0.0062, dᵍˢᵃᵍ = 10))

    c1 = mmc(Vᵈᶜ = 320, Vₘ = 320,
            P_max = 250, P_min = 100, P = -100, Q = 0, Q_max = 500, Q_min = -500, P_dc = 100,
            energy = PI_control(Kₚ = 120, Kᵢ = 400),
            occ = PI_control(ζ = 0.7, bandwidth = 1000),
            ccc = PI_control(ζ = 0.7, bandwidth = 300),
            zcc = PI_control(ζ = 0.7, bandwidth = 300),
            dc = PI_control(Kₚ = 2.0020e-07, Kᵢ = 1.0010e-04)
            )
    c2 = mmc(Vᵈᶜ = 320, Vₘ = 320,
            P_max = 150, P_min = -150, P = 100, Q = 0, Q_max = 300, Q_min = -300, P_dc = -100,
            energy = PI_control(Kₚ = 120, Kᵢ = 400),
            occ = PI_control(ζ = 0.7, bandwidth = 1000),
            ccc = PI_control(ζ = 0.7, bandwidth = 300),
            zcc = PI_control(ζ = 0.7, bandwidth = 300),
            power = PI_control(Kₚ = 2.0020e-07, Kᵢ = 1.0010e-04)
            )

    # connections
    gen1[2.1] ⟷ gen1[2.2] ⟷ gen1[2.3] ⟷ gnd1
    gen1[1.1] ⟷ tl1[1.1]
    gen1[1.2] ⟷ tl1[1.2]
    gen1[1.3] ⟷ tl1[1.3]

    tl1[2.1] ⟷ c1[2.1]
    tl1[2.2] ⟷ c1[2.2]
    tl1[2.3] ⟷ c1[2.3]

    c1[1.1] ⟷ dc_line[1.1]
    c1[1.2] ⟷ dc_line[1.2]

    c2[1.1] ⟷ dc_line[2.1]
    c2[1.2] ⟷ dc_line[2.2]

    c2[2.1] ⟷ tl2[1.1]
    c2[2.2] ⟷ tl2[1.2]
    c2[2.3] ⟷ tl2[1.3]

    gen2[1.1] ⟷ tl2[2.1]
    gen2[1.2] ⟷ tl2[2.2]
    gen2[1.3] ⟷ tl2[2.3]
    gen2[2.1] ⟷ gen2[2.2] ⟷ gen2[2.3] ⟷ gnd2
end
result = power_flow(net)
