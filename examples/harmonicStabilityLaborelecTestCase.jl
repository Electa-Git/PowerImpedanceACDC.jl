using DelimitedFiles,SymEngine
s = symbols("s")
ACvoltage = 333/sqrt(3)
DCvoltage = 640
activePower = 1000
reactivePower1 = 0
reactivePower2 = 0
@time net = @network begin
    
        
        g2 = ac_source(V = ACvoltage, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, P = -activePower, Q = reactivePower2, pins = 3, transformation = true)
        g1 = ac_source(V = ACvoltage, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, P = activePower, Q = reactivePower1, pins = 3, transformation = true)

        # OHL1
        tl1 = overhead_line(length = 50e3,
                conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                                Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
                groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
                earth_parameters = (1,1,100), transformation = true)
        
        tl2 = overhead_line(length = 50e3,
                conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                                Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
                groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
                earth_parameters = (1,1,100), transformation = true)
        
        # t1 = transformer(n = 333/333 , Lₚ = 0.06/2, Rₚ = 0, Rₛ = 0, Lₛ = 0.06/2, Cₜ = 0, Cₛ = 0, pins = 3, transformation = true)
        # t2 = transformer(n = 333/333 , Lₚ = 0.06/2, Rₚ = 0, Rₛ = 0, Lₛ = 0.06/2, Cₜ = 0, Cₛ = 0, pins = 3, transformation = true)
        
        # SI gains using old code
        # c2 = mmc(Vᵈᶜ = DCvoltage, Vₘ = ACvoltage, Rᵣ = 0,
        #         P_max = 1000, P_min = -1000, P = activePower, Q = reactivePower2, Q_max = 1000, Q_min = -1000,
        #         # occ = PI_control(Kₚ = 111.08, Kᵢ = 7.55e+04), # 150 Hz
        #         occ = PI_control(Kₚ = 148.46, Kᵢ = 1.3422e+05), # 200 Hz
        #         ccc = PI_control(Kₚ = 42.9123, Kᵢ = 1.9739e+04),
        #         pll = PI_control(Kₚ = 3.2353e-04, Kᵢ = 0.0145),
        #         p = PI_control(Kₚ = 2.45e-07, Kᵢ = 7.7e-05),
        #         q = PI_control(Kₚ = 2.45e-07, Kᵢ = 7.7e-05),
        #         timeDelay = 116e-6, padeOrderNum = 3, padeOrderDen =3
        #         )

        # c1 = mmc(Vᵈᶜ = DCvoltage, Vₘ = ACvoltage, Rᵣ = 0,
        #         P_max = 1000, P_min = -1000, P = -activePower, Q = reactivePower1, Q_max = 1000, Q_min = -1000,
        #         # occ = PI_control(Kₚ = 111.08, Kᵢ = 7.55e+04),
        #         occ = PI_control(Kₚ = 148.46, Kᵢ = 1.3422e+05), # 200 Hz
        #         ccc = PI_control(Kₚ = 42.9123, Kᵢ = 1.9739e+04),
        #         pll = PI_control(Kₚ = 3.2353e-04, Kᵢ = 0.0145),
        #         dc = PI_control(Kₚ = 0.0192, Kᵢ = 0.0575),
        #         q = PI_control(Kₚ = 2.45e-07, Kᵢ = 7.7e-05)
        #         )
        # PU gains
        c2 = mmc(Vᵈᶜ = DCvoltage, Vₘ = ACvoltage, Rᵣ = 0,
                P_max = 1000, P_min = -1000, P = activePower, Q = reactivePower2, Q_max = 1000, Q_min = -1000,
                occ = PI_control(Kₚ = 0.8359, Kᵢ = 569.24), # 150 Hz
                # occ = PI_control(Kₚ = 1.1178, Kᵢ = 1.012e3), # 200 Hz
                ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
                pll = PI_control(Kₚ = 0.28, Kᵢ =12.5664),
                p = PI_control(Kₚ = 0.1, Kᵢ = 31.4159),
                q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159),
                timeDelay = 116e-6, padeOrderNum = 3, padeOrderDen =3
                )

        c1 = mmc(Vᵈᶜ = DCvoltage, Vₘ = ACvoltage, Rᵣ = 0,
                P_max = 1000, P_min = -1000, P = -activePower, Q = reactivePower1, Q_max = 1000, Q_min = -1000,
                occ = PI_control(Kₚ = 0.8359, Kᵢ = 569.24),
                # occ = PI_control(Kₚ = 1.1178, Kᵢ = 1.012e3), # 200 Hz
                ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
                pll = PI_control(Kₚ = 0.28, Kᵢ =12.5664),
                dc = PI_control(Kₚ = 5, Kᵢ = 15),
                q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
                )
        # Paper Aleksandra with a PLL
        # c2 = mmc(Vᵈᶜ = DCvoltage, Vₘ = ACvoltage, Rᵣ = 0,
        #         P_max = 1000, P_min = -1000, P = activePower, Q = reactivePower2, Q_max = 1000, Q_min = -1000,
        #         occ = PI_control(Kₚ = 1.0634, Kᵢ = 766.5),
        #         ccc = PI_control(Kₚ = 0.0487, Kᵢ = 10.9863),
        #         # pll = PI_control(Kₚ = 0.28, Kᵢ =12.5664),
        #         p = PI_control(Kₚ = 0.0816, Kᵢ = 40.8),
        #         q = PI_control(Kₚ = 0.0816, Kᵢ = 40.8)
        #         )

        # c1 = mmc(Vᵈᶜ = DCvoltage, Vₘ = ACvoltage, Rᵣ = 0,
        #         P_max = 1000, P_min = -1000, P = -activePower, Q = reactivePower1, Q_max = 1000, Q_min = -1000,
        #         occ = PI_control(Kₚ = 1.0634, Kᵢ = 766.5),
        #         ccc = PI_control(Kₚ = 0.0487, Kᵢ = 10.9863),
        #         # pll = PI_control(Kₚ = 0.28, Kᵢ =12.5664),
        #         dc = PI_control(Kₚ = 2.61, Kᵢ = 522/2),
        #         q = PI_control(Kₚ = 0.0816, Kᵢ = 40.8)
        #         )
        dc_line = cable(length = 140e3, positions = [(-0.5,1), (0.5,1)],
            C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
            C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
            C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
            I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
            I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
            I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)


        g2[1.1] ⟷ tl2[1.1] ⟷ Bus1d
        g2[1.2] ⟷ tl2[1.2] ⟷ Bus1q
 
        g1[1.1] ⟷ tl1[1.1] ⟷ Bus2d
        g1[1.2] ⟷ tl1[1.2] ⟷ Bus2q

        # tl2[2.1] ⟷ t2[1.1] ⟷ Bus3d
        # tl2[2.2] ⟷ t2[1.2] ⟷ Bus3q

        # tl1[2.1] ⟷ t1[1.1] ⟷ Bus4d
        # tl1[2.2] ⟷ t1[1.2] ⟷ Bus4q

        tl2[2.1] ⟷ c2[2.1] ⟷ Bus5d
        tl2[2.2] ⟷ c2[2.2] ⟷ Bus5q
 
        tl1[2.1] ⟷ c1[2.1] ⟷ Bus6d
        tl1[2.2] ⟷ c1[2.2] ⟷ Bus6q

        g2[2.1] ⟷ gndd 
        g2[2.2] ⟷ gndq 

        g1[2.1] ⟷ gndd 
        g1[2.2] ⟷ gndq 

        # g2[2.1] ⟷ g2[2.2] ⟷ gnd
        # g1[2.1] ⟷ g1[2.2] ⟷ gnd
        
        c2[1.1] ⟷ dc_line[1.1] ⟷ BusDC1
        c1[1.1] ⟷ dc_line[2.1] ⟷ BusDC2


end


# @time imp_c2, omega_c2 = check_stability(net, net.elements[:c2], direction = :ac, omega_range = (0,4,1000))

# p = bode(imp_c2, omega = omega_c2)

@time Z_MMC_AC, omega = determine_impedance(net, elim_elements=[:tl2], input_pins=Any[:Bus5d,:Bus5q], output_pins=Any[:gndd,:gndq], omega_range = (3,4,1000))
@time Z_BUS_AC, omega = determine_impedance(net, elim_elements=[:c2], input_pins=Any[:Bus5d,:Bus5q], output_pins=Any[:gndd,:gndq], omega_range = (3,4,1000))

L_AC = Z_BUS_AC ./ Z_MMC_AC

@time nyquist_P2P_AC = nyquistplot(L_AC, omega, zoom = "yes", SM = "PM")
display(nyquist_P2P_AC)

# @time Z_MMC_DC, omega = determine_impedance(net, elim_elements=[:dc_line], input_pins=Any[:BusDC2], output_pins=Any[:gndd], omega_range = (-1,4,1000))
# @time Z_BUS_DC, omega = determine_impedance(net, elim_elements=[:c1], input_pins=Any[:BusDC2], output_pins=Any[:gndd], omega_range = (-1,4,1000))

# L_DC = Z_BUS_DC ./ Z_MMC_DC

# @time nyquist_P2P_DC = nyquistplot(L_DC, zoom = "yes")
# display(nyquist_P2P_DC)