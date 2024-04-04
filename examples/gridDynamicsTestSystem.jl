using DelimitedFiles,SymEngine, HVDCstability, Plots
s = symbols("s")
transmissionVoltage = 380 / sqrt(3)
pHVDC1 = 600
pHVDC2 = 600
qC1 = 0
qC2 = 100
qC3 = 0
qC4 = 100

analysis = false
# The P and Q defined here are what is injected into the network. 
# The setpoint of the reactive power controller is minus the value set here. This is adjusted internally, no action here needed.
@time net = @network begin
    
        sg1 = synchronousmachine(V = 1* transmissionVoltage, Vᵃᶜ_base = 380.0, P_min = -2000, P_max = 2000, Q_max = 2000, Q_min = -2000, P = 900, Q = -100)
        # New definition with a 1:1 trafo
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
                                
        # HVDC link 1
        # PU gains
        # MMC1 controls the DC voltage, and is situated at the remote end.
        c1 = mmc(Vᵈᶜ = 800, vDCbase = 800, Vₘ = transmissionVoltage,
                P_max = 1500, P_min = -1500, P = -pHVDC1, Q = qC1, Q_max = 500, Q_min = -500,
                occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
                # ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
                ccc = PI_control(Kₚ = 0.067, Kᵢ = 30.8425),
                pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
                q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [qC1]),
                dc = PI_control(Kₚ = 5, Kᵢ = 15)
                )
        # MMC2 controls P&Q. It is connected to bus 7.
        c2 = mmc(Vᵈᶜ = 800, vDCbase = 800, Vₘ = transmissionVoltage,
                P_max = 1000, P_min = -1000, P = pHVDC1, Q = qC2, Q_max = 1000, Q_min = -1000,
                occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
                # ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
                ccc = PI_control(Kₚ = 0.067, Kᵢ = 30.8425),
                pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
                p = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [pHVDC1]),
                # vac_supp = PI_control(Kₚ = 20, ω_f = 100, ref = [1.0]),
                # q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [qC2])
                vac = PI_control(Kₚ = 0, Kᵢ = 100, ref = [1.0])
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

        # # Converter transformer represented internally within the MMC using a series RL impedance

        # HVDC link 2
        # MMC3 controls the DC voltage, and is situated at the remote end.
        # # PU gains
        c3 = mmc(Vᵈᶜ = 800, vDCbase = 800, Vₘ = transmissionVoltage,
                P_max = 1500, P_min = -1500, P = -pHVDC2, Q = qC3, Q_max = 500, Q_min = -500,
                occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
                # ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
                ccc = PI_control(Kₚ = 0.067, Kᵢ = 30.8425),
                pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
                q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [qC3]),
                dc = PI_control(Kₚ = 5, Kᵢ = 15)
                )
        # MMC4 controls P&Q. It is connected to bus 5.
        c4 = mmc(Vᵈᶜ = 800, vDCbase = 800, Vₘ = transmissionVoltage,
                P_max = 1000, P_min = -1000, P = pHVDC2, Q = qC4, Q_max = 1000, Q_min = -1000,
                occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
                # ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
                ccc = PI_control(Kₚ = 0.067, Kᵢ = 30.8425),
                # pll = PI_control(Kₚ = 1.4, Kᵢ = 314.1593), # 50 Hz BW
                pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664), # 10 Hz BW
                p = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [pHVDC2]),
                vac_supp = PI_control(Kₚ = 20, ω_f = 100, ref = [1.0]),
                q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [qC4])
                # vac = PI_control(Kₚ = 0, Kᵢ = 100, ref = [1.0])
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

        l5[2.1] == gndd
        l5[2.2] == gndq

        # tl64[2.1] == tl96[1.1] == Bus6d
        # tl64[2.2] == tl96[1.2] == Bus6q

        # tl54[2.1] == tl75[1.1] == c4[2.1] == Bus5d
        # tl54[2.2] == tl75[1.2] == c4[2.2] == Bus5q

        # # With HVDC #2, original
        tl54[2.1] == l5[1.1] == tl75[1.1] == c4[2.1] == Bus5d
        tl54[2.2] == l5[1.2] == tl75[1.2] == c4[2.2] == Bus5q
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
        tl75[2.1] == tl78[1.1] == c2[2.1] == Bus7d
        tl75[2.2] == tl78[1.2] == c2[2.2] == Bus7q
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
# MMC = net.elements[:c1]
# plot_data(MMC, omega_range = (0, 4, 1000), scale = :log)

# @time imp_ac1, omega_ac1 = determine_impedance(net, elim_elements=[:g3], input_pins=Any[:Bus3d,:Bus3q], 
# output_pins=Any[:gndd,:gndq], omega_range = (0,3,1000))

# writedlm("./files/imp_dq_MMC_SG.csv",  imp_ac1, ',')
# writedlm("./files/w_dq_MMC_SG.csv",  omega_ac1, ',')

# Stability analysis at HVDC1 terminals - no offshore cable

@time Z_MMC_AC_1, omega = determine_impedance(net, elim_elements=[:tl78,:tl75], input_pins=Any[:Bus7d,:Bus7q], output_pins=Any[:gndd,:gndq], omega_range = (0,2,1000))
@time Z_BUS_AC_1, omega = determine_impedance(net, elim_elements=[:c2], input_pins=Any[:Bus7d,:Bus7q], output_pins=Any[:gndd,:gndq], omega_range = (0,2,1000))

# Stability analysis at HVDC1 terminals - offshore cable
# @time Z_MMC_AC_1, omega = determine_impedance(net, elim_elements=[:offshore_cable], input_pins=Any[:BusOffshored,:BusOffshoreq], output_pins=Any[:gndd,:gndq], omega_range = (-2,2,1000))
# @time Z_BUS_AC_1, omega = determine_impedance(net, elim_elements=[:c2], input_pins=Any[:BusOffshored,:BusOffshoreq], output_pins=Any[:gndd,:gndq], omega_range = (-2,2,1000))

L_AC_1 = Z_BUS_AC_1 ./ Z_MMC_AC_1

nyquist_P2P_AC_1 = HVDCstability.nyquistplot(L_AC_1, omega, zoom = "yes", SM = "PM", 
# title = "Nyquist plots at bus 7 (HVDC 1 terminals), HVDC1 in vAC droop, HVDC2 in PQ")
# title = "Nyquist plots at bus 7 (HVDC 1 terminals), both links in vAC droop")
# savefig("./files/baseCase_HVDC1_HVDC2_droop_bus7.pdf")
title = "Nyquist plots at bus 7 (HVDC 1 terminals), both links in vAC control")
# savefig("./files/baseCase_HVDC1_HVDC2_PI_bus7.pdf")

Ycon = inv.(Z_MMC_AC_1)
Ybus = inv.(Z_BUS_AC_1)
Zcl_bus = inv.(Ybus + Ycon)

fmin = 1e-2
fmax = 1e3

EVD(Zcl_bus, omega, fmin, fmax)

# Stability analysis at HVDC2 terminals

@time Z_MMC_AC_2, omega = determine_impedance(net, elim_elements=[:tl75,:tl54,:l5], input_pins=Any[:Bus5d,:Bus5q], output_pins=Any[:gndd,:gndq], omega_range = (0,2,1000))
@time Z_BUS_AC_2, omega = determine_impedance(net, elim_elements=[:c4], input_pins=Any[:Bus5d,:Bus5q], output_pins=Any[:gndd,:gndq], omega_range = (0,2,1000))

L_AC_2 = Z_BUS_AC_2 ./ Z_MMC_AC_2

nyquist_P2P_AC_2 = HVDCstability.nyquistplot(L_AC_2, omega, zoom = "yes", SM = "PM", 
# title = "Nyquist plots at bus 5 (HVDC 2 terminals), HVDC1 in vAC droop, HVDC2 in PQ")
# title = "Nyquist plots at bus 5 (HVDC 2 terminals), both links in vAC droop")
# savefig("./files/baseCase_HVDC1_HVDC2_droop_bus5.pdf")
title = "Nyquist plots at bus 5 (HVDC 2 terminals), both links in vAC control")
# savefig("./files/baseCase_HVDC1_HVDC2_PI_bus5.pdf")

# Ycon = inv.(Z_MMC_AC_2)
# Ybus = inv.(Z_BUS_AC_2)
# Zcl_bus = inv.(Ybus + Ycon)

# fmin = 1e-2
# fmax = 1e3

# EVD(Zcl_bus, omega, fmin, fmax)


# @time imp_ac2, omega_ac2 = determine_impedance(net, elim_elements=[:g3,:c1,:c2,:dc_line,:g4], input_pins=Any[:Bus3d,:Bus3q], 
# output_pins=Any[:Bus7d,:Bus7q], omega_range = (-2,4,2000))

# writedlm("imp_Z12.csv",  imp_ac2, ',')
# writedlm("w_Z12.csv",  omega_ac2, ',')

# @time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:g3], input_pins=Any[:Bus3d,:Bus3q], 
# output_pins=Any[:gndd,:gndq], omega_range = (-2,4,2000))

# writedlm("imp_dq_MMC_SG.csv",  imp_ac, ',')
# writedlm("w_dq_MMC_SG.csv",  omega_ac, ',')

# @time imp_c1, omega_c1 = check_stability(net, net.elements[:c1], direction = :ac, omega_range = (0,4,1000))
# @time imp_sg1, omega_sg1 = check_stability(net, net.elements[:sg1], direction = :ac, omega_range = (0,4,1000))

# p = bode(imp_sg1, omega = omega_sg1)


# SI gains
        # c2 = mmc(Vᵈᶜ = 640, Vₘ = transmissionVoltage,
        #         P_max = 1000, P_min = -1000, P = pHVDC1, Q = qC2, Q_max = 1000, Q_min = -1000,
        #         occ = PI_control(Kₚ = 111.0624, Kᵢ = 7.5487e+04),
        #         ccc = PI_control(Kₚ = 42.9123, Kᵢ = 1.9739e+04),
        #         pll = PI_control(Kₚ = 2.8351e-04, Kᵢ = 0.0127),
        #         p = PI_control(Kₚ = 2.1487e-07, Kᵢ = 6.7503e-05),
        #         q = PI_control(Kₚ = 2.1487e-07, Kᵢ = 6.7503e-05)
        #         )

        # c1 = mmc(Vᵈᶜ = 640, Vₘ = transmissionVoltage,
        #         P_max = 1500, P_min = -1500, P = -pHVDC1, Q = qC1, Q_max = 500, Q_min = -500,
        #         occ = PI_control(Kₚ = 111.0624, Kᵢ = 7.5487e+04),
        #         ccc = PI_control(Kₚ = 42.9123, Kᵢ = 1.9739e+04),
        #         pll = PI_control(Kₚ = 2.8351e-04, Kᵢ = 0.0127),
        #         q = PI_control(Kₚ = 2.1487e-07, Kᵢ = 6.7503e-05),
        #         dc = PI_control(Kₚ = 0.0168, Kᵢ = 0.0504)
        #         )

        # SI gains
        # c4 = mmc(Vᵈᶜ = 640, Vₘ = transmissionVoltage,
        #         P_max = 1000, P_min = -1000, P = pHVDC2, Q = qC4, Q_max = 1000, Q_min = -1000,
        #         occ = PI_control(Kₚ = 111.0624, Kᵢ = 7.5487e+04),
        #         ccc = PI_control(Kₚ = 42.9123, Kᵢ = 1.9739e+04),
        #         pll = PI_control(Kₚ = 2.8351e-04, Kᵢ = 0.0127),
        #         p = PI_control(Kₚ = 2.1487e-07, Kᵢ = 6.7503e-05),
        #         q = PI_control(Kₚ = 2.1487e-07, Kᵢ = 6.7503e-05)
        #         )

        # c3 = mmc(Vᵈᶜ = 640, Vₘ = transmissionVoltage,
        #         P_max = 1500, P_min = -1500, P = -pHVDC2, Q = qC3, Q_max = 500, Q_min = -500,
        #         occ = PI_control(Kₚ = 111.0624, Kᵢ = 7.5487e+04),
        #         ccc = PI_control(Kₚ = 42.9123, Kᵢ = 1.9739e+04),
        #         pll = PI_control(Kₚ = 2.8351e-04, Kᵢ = 0.0127),
        #         q = PI_control(Kₚ = 2.1487e-07, Kᵢ = 6.7503e-05),
        #         dc = PI_control(Kₚ = 0.0168, Kᵢ = 0.0504)
        #         )

        # Old definition with a MV ideal source and transformer
        # g3 = ac_source(V = 13.8/sqrt(3), P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)
        # t3 = transformer(n = 13.8/380 , Lₚ = 3.5523e-5/2, Rₚ = 0, Rₛ = 0, Lₛ = 0.0269/2, Cₜ = 7e-9, Cₛ = 12e-9, pins = 3, transformation = true)
        # t3 = impedance(z = 0.0269s, pins = 3, transformation = true)

                # l5 = impedance(z = 960 + s, pins = 3, transformation = true) 
        # l6 = impedance(z = 960 + s, pins = 3, transformation = true) 
        # l8 = impedance(z = 960 + s, pins = 3, transformation = true)
end