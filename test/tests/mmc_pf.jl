#include("../../src/HVDCstability.jl")
#using .HVDCstability

@time net = @network begin
    gen1 = ac_source(pins = 3, P_min = 50, P = 1000, P_max = 1500, Q = 0, Q_max = 500, Q_min = -500,
                    V = 320, transformation = true)
    gen2 = ac_source(pins = 3, P_min = -1500, P = -1000, P_max = -50, Q = 0, Q_max = 500, Q_min = -500,
                    V = 320, transformation = true)

    tl1 = overhead_line(length = 200e3,
        conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                        Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        earth_parameters = (1,1,100), transformation = true)
    tl2 = overhead_line(length = 200e3,
        conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                        Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        earth_parameters = (1,1,100), transformation = true)

    dc_line = cable(length = 100e3, positions = [(-0.5,1), (0.5,1)],
          C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
          C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
          C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
          I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
          I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
          I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)

    c1 = mmc(Vᵈᶜ = 640, Vₘ = 230,
            P_max = -50, P_min = -1500, P = -1000, Q = 0, Q_max = 500, Q_min = -500, P_dc = 100,
            occ = PI_control(ζ = 0.7, bandwidth = 1000),
            ccc = PI_control(ζ = 0.7, bandwidth = 300),
            dc = PI_control(Kₚ = 0.01, Kᵢ = 2),
            energy = PI_control(Kₚ = 120, Kᵢ = 400),
            zcc = PI_control(ζ = 0.7, bandwidth = 300)
            )
    c2 = mmc(Vᵈᶜ = 640, Vₘ = 230,
            P_max = 1500, P_min = 50, P = 1000, Q = 0, Q_max = 300, Q_min = -300, P_dc = -100,
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

    c1[1.1] ⟷ dc_line[1.1]
    c2[1.1] ⟷ dc_line[2.1]

    c2[2.1] ⟷ tl2[1.1]
    c2[2.2] ⟷ tl2[1.2]

    gen2[1.1] ⟷ tl2[2.1]
    gen2[1.2] ⟷ tl2[2.2]
    gen2[2.1] ⟷ gen2[2.2] ⟷ gnd2
end

# imp_ac, omega_ac = determine_impedance(net, elim_elements = [:c1],
#         input_pins = Any[(:c1, Symbol(2.1))],
#         output_pins = Any[(:c1, Symbol(2.2))], omega_range = (0, 4, 1000))
# bode(imp_ac, omega = omega_ac)

@time imp, omega = check_stability(net, net.elements[:c1], direction = :dc, omega_range = (0,4,1000))
p = bode(imp, omega = omega, titles = ["Z_{MMC1}" "Z_{eq}" "Y_{MMC1} Z_{eq}"])

n = nyquist(imp, title = "Y_{MMC1} Z_{eq}")
# imp, omega = check_stability(net, net.elements[:c2], direction = :dc)
# bode(imp, omega = omega, titles = ["Z_{MMC2}" "Z_{eq}" "Y_{MMC2} Z_{eq}"])

@time imp, omega = check_stability(net, net.elements[:c1], direction = :ac, omega_range = (0,4,1000))
p = bode(imp, omega = omega, titles = ["Z_{MMC1, AC}" "Z_{eq, AC}" "Y_{MMC1, AC} Z_{eq, AC}"])

n = nyquist(imp, title = "Y_{MMC1, AC} Z_{eq, AC}")
