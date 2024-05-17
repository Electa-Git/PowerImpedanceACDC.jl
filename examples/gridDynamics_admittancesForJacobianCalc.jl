using DelimitedFiles,SymEngine, HVDCstability
s = symbols("s")
transmissionVoltage = 380 / sqrt(3)
pHVDC1 = 600
pHVDC2 = 600
qC1 = 0
qC2 = 100
qC3 = 0
qC4 = 100

# The P and Q defined here are what is injected into the network. 
# The setpoint of the reactive power controller is minus the value set here. This is adjusted internally, no action here needed.
@time net = @network begin

        voltageBase = transmissionVoltage
    
        sg1 = synchronousmachine(V = 1* transmissionVoltage, Vᵃᶜ_base = 380.0, P_min = -2000, P_max = 2000, Q_max = 2000, Q_min = -2000, P = 900)
        g3 = ac_source(V = 380/sqrt(3), P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)
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

        l5 = impedance(z = 1.0*(192.5 * 4.6s)/(192.5 + 4.6s), pins = 3, transformation = true) 
        l6 = impedance(z = 1.0*(192.5 * 4.6s)/(192.5 + 4.6s), pins = 3, transformation = true)
        l8 = impedance(z = 1.0*(192.5 * 4.6s)/(192.5 + 4.6s), pins = 3, transformation = true)

        # HVDC link 1
        # PU gains
        # MMC1 controls the DC voltage, and is situated at the remote end.
        c1 = mmc(Vᵈᶜ = 800, vDCbase = 800, vACbase_LL_RMS = 380, vPCC_scaling = 380/380, Vₘ = transmissionVoltage, Lᵣ = 60e-3, Rᵣ = 0.535,
                P_max = 1500, P_min = -1500, P = -pHVDC1, Q = qC1, Q_max = 500, Q_min = -500,
                occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
                ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
                pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
                q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [qC1]),
                dc = PI_control(Kₚ = 5, Kᵢ = 15, ref = [1.0])
                )
        # MMC2 controls P&Q. It is connected to bus 7.
        c2 = mmc(Vᵈᶜ = 800, vDCbase = 800, vACbase_LL_RMS = 380, vPCC_scaling = 380/380, Vₘ = transmissionVoltage, Lᵣ = 60e-3, Rᵣ = 0.535,
                P_max = 1000, P_min = -1000, P = pHVDC1, Q = qC2, Q_max = 1000, Q_min = -1000,
                occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
                ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
                pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
                p = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [pHVDC1]),
                vac_supp = PI_control(Kₚ = 20, ω_f = 100, ref = [transmissionVoltage * sqrt(2)]),
                q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [qC2])
                # vac = PI_control(Kₚ = 0, Kᵢ = 100, ref = [transmissionVoltage * sqrt(2)])
                )

        dc_line = cable(length = 100e3, positions = [(-0.5,1), (0.5,1)],
            C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
            C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
            C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
            I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
            I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
            I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)

        g4 = ac_source(V = transmissionVoltage, P = pHVDC1, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)

        # TL at the remote end
        tl1 = overhead_line(length = 50e3,
                conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                                Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
                groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
                earth_parameters = (1,1,100), transformation = true)

        c1[2.1] ⟷ tl1[2.1]
        c1[2.2] ⟷ tl1[2.2]

        g4[1.1] ⟷ tl1[1.1]
        g4[1.2] ⟷ tl1[1.2]

        g4[2.1] ⟷ gndd
        g4[2.2] ⟷ gndq
        
        c1[1.1] ⟷ dc_line[1.1]
        c2[1.1] ⟷ dc_line[2.1]

        # HVDC link 2
        # MMC3 controls the DC voltage, and is situated at the remote end.
        # # PU gains
        c3 = mmc(Vᵈᶜ = 800, vDCbase = 800, vACbase_LL_RMS = 380, vPCC_scaling = 380/380, Vₘ = transmissionVoltage, Lᵣ = 60e-3, Rᵣ = 0.535,
                P_max = 1500, P_min = -1500, P = -pHVDC2, Q = qC3, Q_max = 500, Q_min = -500,
                occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
                ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
                pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
                q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [qC3]),
                dc = PI_control(Kₚ = 5, Kᵢ = 15, ref = [1.0])
                )
        # MMC4 controls P&Q. It is connected to bus 5.
        c4 = mmc(Vᵈᶜ = 800, vDCbase = 800, vACbase_LL_RMS = 380, vPCC_scaling = 380/380, Vₘ = transmissionVoltage, Lᵣ = 60e-3, Rᵣ = 0.535,
                P_max = 1000, P_min = -1000, P = pHVDC2, Q = qC4, Q_max = 1000, Q_min = -1000,
                occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
                ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
                pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
                p = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [pHVDC2]),
                vac_supp = PI_control(Kₚ = 20, ω_f = 100, ref = [transmissionVoltage * sqrt(2)]),
                q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [qC4])
                # vac = PI_control(Kₚ = 0, Kᵢ = 100, ref = [transmissionVoltage * sqrt(2)])
                )
        
        
        dc_line_2 = cable(length = 50e3, positions = [(-0.5,1), (0.5,1)],
            C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
            C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
            C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
            I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
            I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
            I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)

        g5 = ac_source(V = transmissionVoltage, P = pHVDC2, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)

        # TL at the remote end
        tl2 = overhead_line(length = 50e3,
                conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                                Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
                groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
                earth_parameters = (1,1,100), transformation = true)

        g5[1.1] ⟷ tl2[1.1]
        g5[1.2] ⟷ tl2[1.2]

        c3[2.1] ⟷ tl2[2.1]
        c3[2.2] ⟷ tl2[2.2]

        g5[2.1] ⟷ gndd
        g5[2.2] ⟷ gndq
        
        c3[1.1] ⟷ dc_line_2[1.1]
        c4[1.1] ⟷ dc_line_2[2.1]

        # With SG
        sg1[2.1] == gndd
        sg1[2.2] == gndq

        sg1[1.1] == tl64[1.1] == tl54[1.1] == Bus4d
        sg1[1.2] == tl64[1.2] == tl54[1.2] == Bus4q

        # Without SG

        # tl64[1.1] == tl54[1.1] == Bus4d
        # tl64[1.2] == tl54[1.2] == Bus4q

        tl64[2.1] == l6[1.1] == tl96[1.1] == Bus6d
        tl64[2.2] == l6[1.2] == tl96[1.2] == Bus6q

        l6[2.1] == gndd
        l6[2.2] == gndq

        # With HVDC #2
        tl54[2.1] == l5[1.1] == tl75[1.1] == c4[2.1] == Bus5d
        tl54[2.2] == l5[1.2] == tl75[1.2] == c4[2.2] == Bus5q
        # Without HVDC #2
        # tl54[2.1] == l5[1.1] == tl75[1.1] == Bus5d
        # tl54[2.2] == l5[1.2] == tl75[1.2] == Bus5q

        l5[2.1] == gndd
        l5[2.2] == gndq

        tl96[2.1] == tl89[1.1] == t3[2.1] == Bus9d
        tl96[2.2] == tl89[1.2] == t3[2.2] == Bus9q

        # With HVDC #1
        tl75[2.1] == tl78[1.1] == c2[2.1] == Bus7d
        tl75[2.2] == tl78[1.2] == c2[2.2] == Bus7q
        # Without HVDC #1
        # tl75[2.1] == tl78[1.1] == Bus7d
        # tl75[2.2] == tl78[1.2] == Bus7q

        tl78[2.1] == tl89[2.1] == l8[1.1] == Bus8d
        tl78[2.2] == tl89[2.2] == l8[1.2] == Bus8q

        l8[2.1] == gndd
        l8[2.2] == gndq

        t3[1.1] == g3[1.1] == Bus3d
        t3[1.2] == g3[1.2] == Bus3q

        g3[2.1] == gndd
        g3[2.2] == gndq


end


@time imp_ac1, omega_ac1 = determine_impedance(net, elim_elements=[:tl75,:tl78], input_pins=Any[:Bus7d,:Bus7q], 
output_pins=Any[:gndd,:gndq], omega_range = (-3,3,1000))

@time imp_ac2, omega_ac2 = determine_impedance(net, elim_elements=[:c2], input_pins=Any[:Bus7d,:Bus7q], 
output_pins=Any[:gndd,:gndq], omega_range = (-3,3,1000))

MMC2_Vm = net.elements[:c2].element_value.Vₘ / net.elements[:c2].element_value.vACbase
MMC2_θ = net.elements[:c2].element_value.θ
MMC2_P = -net.elements[:c2].element_value.P / net.elements[:c2].element_value.Sbase
MMC2_Q = net.elements[:c2].element_value.Q / net.elements[:c2].element_value.Sbase

MMC2_vd0 = MMC2_Vm * cos(MMC2_θ) 
MMC2_vq0 = -MMC2_Vm * sin(MMC2_θ)
MMC2_id0 = ((MMC2_vd0 * MMC2_P + MMC2_vq0 * MMC2_Q) / (MMC2_vd0^2 + MMC2_vq0^2)) 
MMC2_iq0 = ((MMC2_vq0 * MMC2_P - MMC2_vd0 * MMC2_Q) / (MMC2_vd0^2 + MMC2_vq0^2))

writedlm("./files/imp_HVDC1.csv",  imp_ac1, ',')
writedlm("./files/w_HVDC1.csv",  omega_ac1, ',')
writedlm("./files/op_HVDC1.csv", [MMC2_vd0 MMC2_vq0 MMC2_id0 MMC2_iq0], ',')
writedlm("./files/op_HVDC1_rest.csv", [MMC2_vd0 MMC2_vq0 -MMC2_id0 -MMC2_iq0], ',')
writedlm("./files/imp_HVDC1_rest.csv",  imp_ac2, ',')
writedlm("./files/w_HVDC1_rest.csv",  omega_ac2, ',')

@time imp_ac3, omega_ac3 = determine_impedance(net, elim_elements=[:tl75,:tl54,:l5], input_pins=Any[:Bus5d,:Bus5q], 
output_pins=Any[:gndd,:gndq], omega_range = (-3,3,1000))

@time imp_ac4, omega_ac4 = determine_impedance(net, elim_elements=[:c4], input_pins=Any[:Bus5d,:Bus5q], 
output_pins=Any[:gndd,:gndq], omega_range = (-3,3,1000))

MMC4_Vm = net.elements[:c4].element_value.Vₘ / net.elements[:c4].element_value.vACbase
MMC4_θ = net.elements[:c4].element_value.θ
MMC4_P = -net.elements[:c4].element_value.P / net.elements[:c4].element_value.Sbase
MMC4_Q = net.elements[:c4].element_value.Q / net.elements[:c4].element_value.Sbase

MMC4_vd0 = MMC4_Vm * cos(MMC4_θ)
MMC4_vq0 = -MMC4_Vm * sin(MMC4_θ)
MMC4_id0 = ((MMC4_vd0 * MMC4_P + MMC4_vq0 * MMC4_Q) / (MMC4_vd0^2 + MMC4_vq0^2)) 
MMC4_iq0 = ((MMC4_vq0 * MMC4_P - MMC4_vd0 * MMC4_Q) / (MMC4_vd0^2 + MMC4_vq0^2)) 



writedlm("./files/imp_HVDC2.csv",  imp_ac3, ',')
writedlm("./files/w_HVDC2.csv",  omega_ac3, ',')
writedlm("./files/op_HVDC2.csv", [MMC4_vd0 MMC4_vq0 MMC4_id0 MMC4_iq0], ',')
writedlm("./files/op_HVDC2_rest.csv", [MMC4_vd0 MMC4_vq0 -MMC4_id0 -MMC4_iq0], ',')
writedlm("./files/imp_HVDC2_rest.csv",  imp_ac4, ',')
writedlm("./files/w_HVDC2_rest.csv",  omega_ac4, ',')


