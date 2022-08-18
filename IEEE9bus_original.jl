using HVDCstability,SymEngine
s = symbols("s")
@time net = @network begin
    
    g1 = ac_source(pins = 3, transformation = true)
    g2 = ac_source(pins = 3, transformation = true)
    g3 = ac_source(pins = 3, transformation = true)

    tl78 = overhead_line(length = 76e3,
            conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                            Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
            groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
            earth_parameters = (1,1,100), transformation = true)

    tl89 = overhead_line(length = 106e3,
            conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                            Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
            groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
            earth_parameters = (1,1,100), transformation = true)  

    tl75 = overhead_line(length = 170e3,
            conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                            Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
            groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
            earth_parameters = (1,1,100), transformation = true)

    tl96 = overhead_line(length = 179e3,
            conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                            Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
            groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
            earth_parameters = (1,1,100), transformation = true) 

    tl54 = overhead_line(length = 89e3,
            conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                            Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
            groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
            earth_parameters = (1,1,100), transformation = true)

    tl64 = overhead_line(length = 97e3,
            conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                            Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
            groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
            earth_parameters = (1,1,100), transformation = true)
            
            
    t1 = transformer(V₁ᵒ = 2.4e3, V₁ˢ = 51.87, V₂ᵒ = 240, P₁ᵒ = 171.1, P₁ˢ = 642.1,
            I₁ᵒ = 0.48, I₁ˢ = 20.83, Cₛ = 12e-6, Cₜ = 7e-6, pins = 3, transformation = true)

    t2 = transformer(V₁ᵒ = 2.4e3, V₁ˢ = 51.87, V₂ᵒ = 240, P₁ᵒ = 171.1, P₁ˢ = 642.1,
            I₁ᵒ = 0.48, I₁ˢ = 20.83, Cₛ = 12e-6, Cₜ = 7e-6, pins = 3, transformation = true)

    t3 = transformer(V₁ᵒ = 2.4e3, V₁ˢ = 51.87, V₂ᵒ = 240, P₁ᵒ = 171.1, P₁ˢ = 642.1,
            I₁ᵒ = 0.48, I₁ˢ = 20.83, Cₛ = 12e-6, Cₜ = 7e-6, pins = 3, transformation = true)

    l5 = impedance(z = 423.2+1/3.36s/s, pins = 3, transformation = true) # Nominal P 125 MW, nominal Q 50 MVar
    l6 = impedance(z = 587.78+1/5.61s/s, pins = 3, transformation = true) # Nominal P 125 MW, nominal Q 50 MVar
    l8 = impedance(z = 529+1/4.81s/s, pins = 3, transformation = true) # Nominal P 125 MW, nominal Q 50 MVar


    g1[1.1] == gnd1d
    g1[1.2] == gnd1q

    g1[2.1] == t1[1.1] == Bus1d
    g1[2.2] == t1[1.2] == Bus1q

    t1[2.1] == tl64[1.1] == tl54[1.1] == Bus4d
    t1[2.2] == tl64[1.2] == tl54[1.2] == Bus4q

    tl64[2.1] == l6[1.1] == tl96[1.1] == Bus6d
    tl64[2.2] == l6[1.2] == tl96[1.2] == Bus6q

    l6[2.1] == gnd6d
    l6[2.2] == gnd6q

    tl54[2.1] == l5[1.1] == tl75[1.1] == Bus5d
    tl54[2.2] == l5[1.2] == tl75[1.2] == Bus5q

    l5[2.1] == gnd5d
    l5[2.2] == gnd5q

    tl96[2.1] == tl89[1.1] == t3[1.1] == Bus9d
    tl96[2.2] == tl89[1.2] == t3[1.2] == Bus9q

    tl75[2.1] == tl78[1.1] == t2[1.1] == Bus7d
    tl75[2.2] == tl78[1.2] == t2[1.2] == Bus7q

    tl78[2.1] == tl89[2.1] == l8[1.1] == Bus8d
    tl78[2.2] == tl89[2.2] == l8[1.2] == Bus8q

    l8[2.1] == gnd8d
    l8[2.2] == gnd8q

    t2[2.1] == g2[1.1] == Bus2d
    t2[2.2] == g2[1.2] == Bus2q

    g2[2.1] == gnd2d
    g2[2.2] == gnd2q

    t3[2.1] == g3[1.1] == Bus3d
    t3[2.2] == g3[1.2] == Bus3q

    g3[2.1] == gnd3d
    g3[2.2] == gnd3q

end

#@time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:gen1], input_pins=Any[:TestNode1,:TestNode2,:TestNode3], 
#output_pins=Any[:gnd,:gnd,:gnd], omega_range = (0,4,1000))

@time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:g1], input_pins=Any[:Bus1d,:Bus1q], 
output_pins=Any[:gnd3d,:gnd3q], omega_range = (0,4,1000))

# imp_ac_d = []
# imp_ac_q = []

# for n in eachindex(imp_ac)
#     push!(imp_ac_d,imp_ac[n,1][1])
#     push!(imp_ac_q,imp_ac[n,1][2])
# end
p = bode(imp_ac, omega = omega_ac)

save_plot(p)
