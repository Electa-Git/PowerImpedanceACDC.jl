net = @network begin
    vs = dc_source(voltage = 5)
    l = transmission_line(length = 227e3, conductors = Conductors(nᵇ = 2, nˢᵇ = 2, organization = :flat,
    Rᵈᶜ = 0.06266, rᶜ = 0.01436, yᵇᶜ = 27.5, Δxᵇᶜ = 11.8, dˢᵇ = 0.4572, dˢᵃᵍ = 10),
    earth_parameters = (1,1,100),
    groundwires = Groundwires(nᵍ = 2, Δxᵍ = 6.5, Δyᵍ = 7.5, Rᵍᵈᶜ = 0.9196, rᵍ = 0.0062, dᵍˢᵃᵍ = 10))
    vs[1.1] ⟷ l[1.1] ⟷ Node1
    vs[2.1] ⟷ l[1.2] ⟷ Node2
    l[2.1] ⟷ gnd
    l[2.2] ⟷ gnd1
end
imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1, :Node2],
                            output_pins = Any[:gnd, :gnd1], omega_range = (-1, 6, 5000))
bode(imp, omega_range = (-1, 6, 5000))
