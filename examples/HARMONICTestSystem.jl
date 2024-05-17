using DelimitedFiles,SymEngine, HVDCstability, Plots
s = symbols("s")
transmissionVoltage = 380 / sqrt(3)


analysis = false
@time net = @network begin

        voltageBase = transmissionVoltage
    
        # g1 = synchronousmachine(V = 1* transmissionVoltage, Vᵃᶜ_base = 380.0, P_min = -2000, P_max = 2000, Q_max = 2000, Q_min = -2000, P = 900, Q = -100)
        g1 = ac_source(V = transmissionVoltage, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)
        g3 = ac_source(V = transmissionVoltage, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)
        t3 = transformer(n = 380/380 , Lₚ = 0.0269/2, Rₚ = 0, Rₛ = 0, Lₛ = 0.0269/2, Cₜ = 7e-9, Cₛ = 12e-9, pins = 3, transformation = true)

        # Gezelle - Horta OHL
        # tl78 = overhead_line(length = 30e3,
        #         conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
        #                         Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        #         groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        #         earth_parameters = (1,1,100), transformation = true)
        # Gezelle - Horta cable
        tl78 = cable(length = 30e3, positions = [(-0.5,1.9), (0,1.9), (0.5,1.9)],
                C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
                I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
                C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
                I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
                C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
                I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true, earth_parameters = (1,1,100))

        # offshore_cable = cable(length = 30e3, positions = [(-0.5,1.9), (0,1.9), (0.5,1.9)],
        #         C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
        #         I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
        #         C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
        #         I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
        #         C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
        #         I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true, earth_parameters = (1,1,100))

    
        # Horta - Mercator
        tl89 = overhead_line(length = 50e3,
                conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                                Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
                groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
                earth_parameters = (1,1,100), transformation = true)  
    
        # Horta - Avelgem OHL
        tl75 = overhead_line(length = 50e3,
                conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                                Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
                groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
                earth_parameters = (1,1,100), transformation = true)

        # Horta - Avelgem cable
        # tl75 = cable(length = 50e3, positions = [(-0.5,1.9), (0,1.9), (0.5,1.9)],
        #         C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
        #         I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
        #         C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
        #         I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
        #         C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
        #         I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true, earth_parameters = (1,1,100))

        # tl75_p = overhead_line(length = 50e3,
        #         conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
        #                         Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        #         groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        #         earth_parameters = (1,1,100), transformation = true)
    
        # Mercator - Van Eyck
        tl96 = overhead_line(length = 90e3,
                conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                                Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
                groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
                earth_parameters = (1,1,100), transformation = true) 
    
        # Tihange - Avelgem 
        tl54 = overhead_line(length = 120e3,
                conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                                Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
                groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
                earth_parameters = (1,1,100), transformation = true)
    
        # Van Eyck - Tihange
        tl64 = overhead_line(length = 70e3,
                conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                                Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
                groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
                earth_parameters = (1,1,100), transformation = true)
        
        # New line for meshing
        # tl58 = overhead_line(length = 70e3,
        #         conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
        #                         Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        #         groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        #         earth_parameters = (1,1,100), transformation = true)
        # tl58 = cable(length = 70e3, positions = [(-0.5,1.9), (0,1.9), (0.5,1.9)],
        #         C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
        #         I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
        #         C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
        #         I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
        #         C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
        #         I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true, earth_parameters = (1,1,100))
        

        l5 = impedance(z = 1*(192.5 * 4.6s)/(192.5 + 4.6s), pins = 3, transformation = true) 
        l6 = impedance(z = 1*(192.5 * 4.6s)/(192.5 + 4.6s), pins = 3, transformation = true)
        l8 = impedance(z = 1*(192.5 * 4.6s)/(192.5 + 4.6s), pins = 3, transformation = true)

        # l5[1.1] == l6[1.1] == Node1A
        # l5[1.2] == l6[1.2] == Node1B
        # l5[1.3] == l6[1.3] == Node1C

        # l8[1.1] == l6[2.1] == Node2A
        # l8[1.2] == l6[2.2] == Node2B
        # l8[1.3] == l6[2.3] == Node2C

        # l5[2.1] == l8[2.1] == Node3A
        # l5[2.2] == l8[2.2] == Node3B
        # l5[2.3] == l8[2.3] == Node3C

        

                                
        # MMC2 controls P&Q. It is connected to bus 7.
        g2 = ac_source(V = transmissionVoltage, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)
        g4 = ac_source(V = transmissionVoltage, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)
       

        g1[2.1] == gndd
        g1[2.2] == gndq
        g2[2.1] == gndd
        g2[2.2] == gndq
        g4[2.1] == gndd
        g4[2.2] == gndq

        g1[1.1] == tl64[1.1] == tl54[1.1] == "Bus4d"
        g1[1.2] == tl64[1.2] == tl54[1.2] == "Bus4q"

        # Without SG

        # tl64[1.1] == tl54[1.1] == Bus4d
        # tl64[1.2] == tl54[1.2] == Bus4q

        tl64[2.1] == l6[1.1] == tl96[1.1] == Bus6d
        tl64[2.2] == l6[1.2] == tl96[1.2] == Bus6q

        l6[2.1] == gndd
        l6[2.2] == gndq

        l5[2.1] == gndd
        l5[2.2] == gndq

        # tl64[2.1] == tl96[1.1] == Bus6d
        # tl64[2.2] == tl96[1.2] == Bus6q

        # tl54[2.1] == tl75[1.1] == c4[2.1] == Bus5d
        # tl54[2.2] == tl75[1.2] == c4[2.2] == Bus5q

        # # With HVDC #2, original
        tl54[2.1] == l5[1.1] == tl75[1.1] == g4[1.1] == Bus5d
        tl54[2.2] == l5[1.2] == tl75[1.2] == g4[1.2] == Bus5q
        # With HVDC #2 and new line
        # tl54[2.1] == l5[1.1] == tl75[1.1] == c4[2.1] == tl58[1.1] == Bus5d
        # tl54[2.2] == l5[1.2] == tl75[1.2] == c4[2.2] == tl58[1.2] == Bus5q
        # With HVDC #2, new line and parallel line
        # tl54[2.1] == l5[1.1] == tl75[1.1] == tl75_p[1.1] == c4[2.1] == tl58[1.1] == Bus5d
        # tl54[2.2] == l5[1.2] == tl75[1.2] == tl75_p[1.2] == c4[2.2] == tl58[1.2] == Bus5q
        # Without HVDC #2
        # tl54[2.1] == l5[1.1] == tl75[1.1] == Bus5d
        # tl54[2.2] == l5[1.2] == tl75[1.2] == Bus5q

        tl96[2.1] == tl89[1.1] == t3[2.1] == Bus9d
        tl96[2.2] == tl89[1.2] == t3[2.2] == Bus9q

        # Offshore cable between HVDC #1 and Bus7

        # offshore_cable[1.1] == c2[2.1] == BusOffshored
        # offshore_cable[1.2] == c2[2.2] == BusOffshoreq


        # tl75[2.1] == tl78[1.1] == offshore_cable[2.1] == Bus7d
        # tl75[2.2] == tl78[1.2] == offshore_cable[2.2] == Bus7q


        # With HVDC #1
        tl75[2.1] == tl78[1.1] == g2[1.1] == Bus7d
        tl75[2.2] == tl78[1.2] == g2[1.2] == Bus7q
        # tl75[2.1] == tl75_p[2.1] == tl78[1.1] == c2[2.1] == Bus7d
        # tl75[2.2] == tl75_p[2.2] == tl78[1.2] == c2[2.2] == Bus7q
        # Without HVDC #1
        # tl75[2.1] == tl78[1.1] == Bus7d
        # tl75[2.2] == tl78[1.2] == Bus7q

        # Original, without new line
        tl78[2.1] == tl89[2.1] == l8[1.1] == Bus8d
        tl78[2.2] == tl89[2.2] == l8[1.2] == Bus8q
        # With new line
        # tl78[2.1] == tl89[2.1] == l8[1.1] == tl58[2.1] == Bus8d
        # tl78[2.2] == tl89[2.2] == l8[1.2] == tl58[2.2] == Bus8q

        l8[2.1] == gndd
        l8[2.2] == gndq

        # tl78[2.1] == tl89[2.1] == Bus8d
        # tl78[2.2] == tl89[2.2] == Bus8q

        t3[1.1] == g3[1.1] == Bus3d
        t3[1.2] == g3[1.2] == Bus3q

        g3[2.1] == gndd
        g3[2.2] == gndq


end

if analysis


# Impedance calculation at bus 7 and bus 5

@time Z_BUS_AC_7, omega = determine_impedance(net, elim_elements=[:g2], input_pins=Any[:Bus7d,:Bus7q], output_pins=Any[:gndd,:gndq], omega_range = (0,4,1000))
@time Z_BUS_AC_5, omega = determine_impedance(net, elim_elements=[:g4], input_pins=Any[:Bus5d,:Bus5q], output_pins=Any[:gndd,:gndq], omega_range = (0,4,1000))

Z_BUS_dd_7=[]
Z_BUS_dq_7=[]
Z_BUS_qd_7=[]
Z_BUS_qq_7=[]
Z_BUS_dd_5=[]
Z_BUS_dq_5=[]
Z_BUS_qd_5=[]
Z_BUS_qq_5=[]
for x in 1:size(Z_BUS_AC_7,1)
        temp1 = Z_BUS_AC_7[x,1]
        temp2 = Z_BUS_AC_5[x,1]
        push!(Z_BUS_dd_7,temp1[1,1])
        push!(Z_BUS_dq_7,temp1[1,2])
        push!(Z_BUS_qd_7,temp1[2,1])
        push!(Z_BUS_qq_7,temp1[2,2])
        push!(Z_BUS_dd_5,temp2[1,1])
        push!(Z_BUS_dq_5,temp2[1,2])
        push!(Z_BUS_qd_5,temp2[2,1])
        push!(Z_BUS_qq_5,temp2[2,2])
end


p = bodeplot(Z_BUS_dq_7, omega, legend=[""])
# save_plot(p, "files/HARMONIC_networkImpedance")

end