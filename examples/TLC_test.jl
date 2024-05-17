using DelimitedFiles,SymEngine, HVDCstability
s = symbols("s")
transmissionVoltage = 220 / sqrt(3)
pHVDC1 = -250
qC1 = -50
qC2 = 0


# The P and Q defined here are what is injected into the network. 
# The setpoint of the reactive power controller is minus the value set here. This is adjusted internally, no action here needed.
@time net = @network begin
    
        voltageBase = transmissionVoltage
        g1 = ac_source(V = transmissionVoltage, P = pHVDC1, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)

        c1 = tlc(Vᵈᶜ = 3.592584956081994e+02, Vₘ = transmissionVoltage, Lᵣ = 0.024649917586073, Rᵣ = 0.07744,  
                P_max = 1500, P_min = -1500, P = -pHVDC1, Q = qC1, Q_max = 500, Q_min = -500,
                occ = PI_control(Kₚ = 0.254647908947033, Kᵢ = 0.8),
                pll = PI_control(Kₚ = 0.795774715459477, Kᵢ = 31.830988618379067, ω_f = 2*pi*80),
                v_meas_filt = PI_control(ω_f = 1e4),
                i_meas_filt = PI_control(ω_f = 1e4),
                vac_supp = PI_control(ω_f = 1/0.5, Kₚ =10),
                f_supp = PI_control(ω_f = 1/0.5, Kₚ =10),
                p = PI_control(Kₚ = 0.01, Kᵢ = 10),
                q = PI_control(Kₚ = 0.01, Kᵢ = 10)
                )
        c2 = mmc(Vᵈᶜ = 3.592584956081994e+02, vDCbase = 3.592584956081994e+02, Vₘ = transmissionVoltage,
                P_max = 1000, P_min = -1000, P = pHVDC1, Q = qC2, Q_max = 1000, Q_min = -1000,
                occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
                ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
                pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
                dc = PI_control(Kₚ = 5, Kᵢ = 15, ref = [1.0]),
                q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159, ref = [qC2])
                )

        g2 = ac_source(V = transmissionVoltage, P = pHVDC1, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)

        # TL at the remote end
        tl1 = overhead_line(length = 50e3,
                conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                                Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
                groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
                earth_parameters = (1,1,100), transformation = true)

        c2[2.1] ⟷ tl1[2.1]
        c2[2.2] ⟷ tl1[2.2]

        g2[1.1] ⟷ tl1[1.1]
        g2[1.2] ⟷ tl1[1.2]

        g2[2.1] ⟷ gndd
        g2[2.2] ⟷ gndq
        
        c1[1.1] ⟷ c2[1.1]

        g1[1.1] == c1[2.1] == Bus7d
        g1[1.2] == c1[2.2] == Bus7q

        g1[2.1] == gndd
        g1[2.2] == gndq


end

writedlm("./files/A_TLC.csv",  net.elements[:c1].element_value.A, ',')
writedlm("./files/B_TLC.csv",  net.elements[:c1].element_value.B, ',')
writedlm("./files/C_TLC.csv",  net.elements[:c1].element_value.C, ',')
writedlm("./files/D_TLC.csv",  net.elements[:c1].element_value.D, ',')



# @time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:g1], input_pins=Any[:Bus7d,:Bus7q], 
# output_pins=Any[:gndd,:gndq], omega_range = (-2,4,2000))

# writedlm("./files/imp_P2P.csv",  imp_ac, ',')
# writedlm("./files/w_P2P.csv",  omega_ac, ',')
