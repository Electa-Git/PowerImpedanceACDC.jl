using DelimitedFiles,SymEngine
s = symbols("s")
@time net = @network begin
    
        # The SG branch is added twice here!
        sg1 = synchronousmachine(V = 1.0005 * 380/sqrt(3), Vᵃᶜ_base = 380.0, P = 900)
        g2 = ac_source(V = 380/sqrt(3), pins = 3, transformation = true)

        sg1[2.1] == gndd
        sg1[2.2] == gndq

        sg1[1.1] == g2[1.1] == Bus1d
        sg1[1.2] == g2[1.2] == Bus1q

        g2[2.1] == gndd
        g2[2.2] == gndq

        # g1 = ac_source(V = 380 * sqrt(2/3), P_min = -2000, P_max = 2000, Q_max = 2000, Q_min = -2000, pins = 3, P = 300, transformation = true)
        # g2 = ac_source(V = 380 * sqrt(2/3), P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)

        # br = impedance(z = 3.7544 + 0.29876565917211*s, pins = 3, transformation = true) 

        # g1[2.1] == gndd
        # g1[2.2] == gndq

        # g1[1.1] == br[1.1] == Bus1d
        # g1[1.2] == br[1.2] == Bus1q

        # br[2.1] == g2[1.1] == Bus2d
        # br[2.2] == g2[1.2] == Bus2q

        # g2[2.1] == gndd
        # g2[2.2] == gndq

        

end


@time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:g2], input_pins=Any[:Bus1d,:Bus1q], 
output_pins=Any[:gndd,:gndq], omega_range = (-2,4,2000))

# p = bode(imp_ac, omega = omega_ac)

writedlm("imp_dq_MMC_SG.csv",  imp_ac, ',')
writedlm("w_dq_MMC_SG.csv",  omega_ac, ',')

# @time imp_c1, omega_c1 = check_stability(net, net.elements[:c1], direction = :ac, omega_range = (0,4,1000))
# @time imp_sg1, omega_sg1 = check_stability(net, net.elements[:sg1], direction = :ac, omega_range = (0,4,1000))

# p = bode(imp_sg1, omega = omega_sg1)
