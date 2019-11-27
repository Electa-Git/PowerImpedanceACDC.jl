net = @network begin
    vs = ac_source(V = 230, P = 80, pins = 3)
    tl = overhead_line(length = 200e3,
                        conductors = Conductors(organization = :flat,
                            nᵇ = 3, Rᵈᶜ = 0.063,
                            rᶜ = 0.015,  yᵇᶜ = 30,
                            Δxᵇᶜ = 10, dˢᵃᵍ    = 10),
                        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92,
                            rᵍ      = 0.0062, Δxᵍ     = 6.5, Δyᵍ     = 7.5, dᵍˢᵃᵍ   = 10),
                        earth_parameters = (1,1,100))

    vs[2.1] ⟷ vs[2.2] ⟷ vs[2.3] ⟷ gnd
    vs[1.1] ⟷ tl[1.1] ⟷ Node11
    vs[1.2] ⟷ tl[1.2] ⟷ Node12
    vs[1.3] ⟷ tl[1.3] ⟷ Node13
    tl[2.1] ⟷ gnd1
    tl[2.2] ⟷ gnd2
    tl[2.3] ⟷ gnd3
end

imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node11, :Node12, :Node13],
                            output_pins = Any[:gnd1, :gnd2, :gnd3], omega_range = (0, 5, 1000))
bode(imp, omega_range = (0, 5, 1000))
