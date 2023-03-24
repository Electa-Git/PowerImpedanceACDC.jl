using DelimitedFiles,SymEngine
s = symbols("s")
@time net = @network begin

        # Values used in the SG validation
    
        # sg1 = synchronousmachine(V = 1.0004 * 380/sqrt(3),  Vᵃᶜ_base = 380.0, P = 900, P_max = 1000)
        # g2 = ac_source(V = 380/sqrt(3), pins = 3, transformation = true)

        # sg1[2.1] == gndd
        # sg1[2.2] == gndq

        # sg1[1.1] == g2[1.1] == Bus1d
        # sg1[1.2] == g2[1.2] == Bus1q

        # g2[2.1] == gndd
        # g2[2.2] == gndq

    
        g2 = ac_source(V = 380/sqrt(3), P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)

    
        # Mercator - Van Eyck
        tl96 = overhead_line(length = 90e3,
                conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                                Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
                groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
                earth_parameters = (1,1,100), transformation = true) 
                
        # tl96_otherside = overhead_line(length = 90e3,
        #         conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
        #                         Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        #         groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        #         earth_parameters = (1,1,100), transformation = true) 

        c1 = mmc(Vᵈᶜ = 640, Vₘ = 380/sqrt(3),
                P_max = 1000, P_min = -1000, P = 400, Q = -400, Q_max = 1000, Q_min = -1000,
                occ = PI_control(Kₚ = 111.0624, Kᵢ = 7.5487e+04),
                ccc = PI_control(Kₚ = 42.9123, Kᵢ = 1.9739e+04),
                pll = PI_control(Kₚ = 2.8351e-04, Kᵢ = 0.0127),
                p = PI_control(Kₚ = 2.1487e-07, Kᵢ = 6.7503e-05),
                q = PI_control(Kₚ = 2.1487e-07, Kᵢ = 6.7503e-05)
                )

        # c2 = mmc(Vᵈᶜ = 640, Vₘ = 380/sqrt(3),
        #         P_max = 1500, P_min = -1500, P = -400, Q = 0, Q_max = 500, Q_min = -500,
        #         # P_max = -50, P_min = -1500, P = -400, Q = 0, Q_max = 500, Q_min = -500,
        #         occ = PI_control(Kₚ = 111.0624, Kᵢ = 7.5487e+04),
        #         ccc = PI_control(Kₚ = 42.9123, Kᵢ = 1.9739e+04),
        #         pll = PI_control(Kₚ = 2.8351e-04, Kᵢ = 0.0127),
        #         q = PI_control(Kₚ = 2.1487e-07, Kᵢ = 6.7503e-05),
        #         dc = PI_control(Kₚ = 0.0168, Kᵢ = 0.0504)
        #         )

        # dc_line = cable(length = 100e3, positions = [(-0.5,1), (0.5,1)],
        #     C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
        #     C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
        #     C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
        #     I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
        #     I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
        #     I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)

        # g4 = ac_source(V = 380/sqrt(3), P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)

        dc = dc_source(V = 640, P_min = -2000, P_max = 2000, P = 500)

        dc[1.1] == gnd

        dc[2.1] == c1[1.1]

        # g4[1.1] == tl96_otherside[2.1]
        # g4[1.2] == tl96_otherside[2.2]
        
        # tl96_otherside[1.1] == c2[2.1]
        # tl96_otherside[1.2] == c2[2.2]
     
        # c1[1.1] == dc_line[1.1]
        # c2[1.1] == dc_line[2.1]

        c1[2.1] == tl96[1.1] == Bus4d
        c1[2.2] == tl96[1.2] == Bus4q

        tl96[2.1] == g2[1.1] == Bus1d
        tl96[2.2] == g2[1.2] == Bus1q

        # g4[2.1] == g4[2.2] == gnd
        # g3[2.1] == g3[2.2] == gnd
        # g4[2.1] == gndd
        # g4[2.2] == gndq
        g2[2.1] == gndd
        g2[2.2] == gndq


end

@time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:g2], input_pins=Any[:Bus1d,:Bus1q], 
output_pins=Any[:gndd,:gndq], omega_range = (-2,4,2000))

# p = bode(imp_ac, omega = omega_ac)

# For dq impedances: [dd qd dq qq]
writedlm("imp_dq_MMC_SG.csv",  imp_ac, ',')
writedlm("w_dq_MMC_SG.csv",  omega_ac, ',')

# @time imp_c1, omega_c1 = check_stability(net, net.elements[:c1], direction = :ac, omega_range = (0,4,1000))
# @time imp_sg1, omega_sg1 = check_stability(net, net.elements[:sg1], direction = :ac, omega_range = (0,4,1000))

# p = bode(imp_sg1, omega = omega_sg1)
