using DelimitedFiles,SymEngine, HVDCstability
s = symbols("s")
voltage = 380/sqrt(3)
pHVDC = 100
@time net = @network begin

        voltageBase = voltage
        g2 = ac_source(V = voltage, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)
        
        # PU gains
        c2 = mmc(Vᵈᶜ = 640, Vₘ = voltage,
                P_max = 1000, P_min = -1000, P = pHVDC, Q = 0, Q_max = 1000, Q_min = -1000,
                occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
                ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
                pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
                p = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [pHVDC]),
                q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [0.0])
                )

        c1 = mmc(Vᵈᶜ = 640, Vₘ = voltage,
                P_max = 1500, P_min = -1500, P = -pHVDC, Q = 0, Q_max = 500, Q_min = -500,
                occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
                ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
                pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
                q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [0.0]),
                dc = PI_control(Kₚ = 5, Kᵢ = 15, ref = [1.0])
                )

        c1[1.1] == c2[1.1] == BusDC

        g4 = ac_source(V = voltage, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)

        c1[2.1] == g4[1.1] == Bus5d
        c1[2.2] == g4[1.2] == Bus5q

        # Excluding TL96
        c2[2.1] == g2[1.1] == Bus1d
        c2[2.2] == g2[1.2] == Bus1q

        g4[2.1] == gndd
        g4[2.2] == gndq
        g2[2.1] == gndd
        g2[2.2] == gndq


end

@time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:g2], input_pins=Any[:Bus1d,:Bus1q], 
output_pins=Any[:gndd,:gndq], omega_range = (-2,4,2000))
# @time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:g2], input_pins=Any[:Bus1d], 
# output_pins=Any[:Bus1q], omega_range = (-2,4,2000))
# p = bode(imp_ac, omega = omega_ac)

# For dq impedances: [dd qd dq qq]
writedlm("imp_dq_MMC_SG.csv",  imp_ac, ',')
writedlm("w_dq_MMC_SG.csv",  omega_ac, ',')

# @time imp_c1, omega_c1 = check_stability(net, net.elements[:c1], direction = :ac, omega_range = (0,4,1000))
# @time imp_sg1, omega_sg1 = check_stability(net, net.elements[:sg1], direction = :ac, omega_range = (0,4,1000))

# p = bode(imp_sg1, omega = omega_sg1)

writedlm("A_julia.csv",  net.elements[:c2].element_value.A, ',')
writedlm("B_julia.csv",  net.elements[:c2].element_value.B, ',')
writedlm("C_julia.csv",  net.elements[:c2].element_value.C, ',')
writedlm("D_julia.csv",  net.elements[:c2].element_value.D, ',')