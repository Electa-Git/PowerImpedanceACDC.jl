net = @network begin
    vs = ac_source(V = 230, P = 80, pins = 3)
    tl = overhead_line(length = 17e3,
                        conductors = Conductors(organization = :vertical,
                            nᵇ = 3, Rᵈᶜ = 0.05, rᶜ = 0.016,  yᵇᶜ = 30, Δxᵇᶜ = 7, Δyᵇᶜ= 9, dˢᵃᵍ = 10),
                        groundwires = Groundwires(nᵍ = 1, Rᵍᵈᶜ = 0.1, rᵍ = 0.011, Δyᵍ = 27, dᵍˢᵃᵍ = 10),
                        earth_parameters = (1,1,25))

    vs[2.1] ⟷ vs[2.2] ⟷ vs[2.3] ⟷ gnd
    vs[1.1] ⟷ tl[1.1] ⟷ Node11
    vs[1.2] ⟷ tl[1.2] ⟷ Node12
    vs[1.3] ⟷ tl[1.3] ⟷ Node13
    tl[2.1] ⟷ gnd1
    tl[2.2] ⟷ gnd2
    tl[2.3] ⟷ gnd3
end

#imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node11, :Node12, :Node13],
                            #output_pins = Any[:gnd1, :gnd2, :gnd3], omega_range = (0, 5, 1000))
#bode(imp, omega_range = (0, 5, 1000))
