using SymEngine,DelimitedFiles
s = symbols("s")
net = @network begin
    vs = ac_source(V = 380e3, pins = 3)
    cb = crossbonded_cable(section = cable(length = 834, positions = [(0,1.9), (0.5,1.9), (1,1.9)],
            C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8, μᵣ = 1),
            # I1 = Insulator(rᵢ = 31.75e-3, a = 31.75e-3 + 0.002, b = 60.85e-3-0.0013, rₒ = 60.85e-3, ϵᵣ = 2.26),
            I1 = Insulator(rᵢ = 31.75e-3, a = 31.75e-3, b = 60.85e-3, rₒ = 60.85e-3, ϵᵣ = 2.58908112),
            C2 = Conductor(rᵢ = 60.85e-3, rₒ = 61.05e-3, ρ = 1.72e-8, μᵣ = 1),
            I2 = Insulator(rᵢ = 61.05e-3, rₒ = 65.95e-3, ϵᵣ = 2.26),
            earth_parameters = (1,1,25)),
             nₛ = 3, #minor sections
             mₛ = 2, #major sections
             Zᶜᵇ = 1e-6s,  # cross-bonding impedance, zero corresponds to a direct connection
             Zₛ = 1e-3  # sheath grounding impedance
             ) 
    # cb = cable(length = 100e3, positions = [(0,1.9), (0.5,1.9), (1,1.9)],
    #          C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8, μᵣ = 1),
    #          I1 = Insulator(rᵢ = 31.75e-3, a = 31.75e-3 + 0.002, b = 60.85e-3-0.0013, rₒ = 60.85e-3, ϵᵣ = 2.26),
    #          C2 = Conductor(rᵢ = 60.85e-3, rₒ = 61.05e-3, ρ = 1.72e-8, μᵣ = 1),
    #          I2 = Insulator(rᵢ = 61.05e-3, rₒ = 65.95e-3, ϵᵣ = 2.26),
    #          earth_parameters = (1,1,25))


    r = impedance(pins = 3, z = 1e-3)

    vs[1.1] ⟷ cb[1.1] ⟷ Node1
    vs[1.2] ⟷ cb[1.2] ⟷ Node2
    vs[1.3] ⟷ cb[1.3] ⟷ Node3
    cb[2.1] ⟷ r[1.1]
    cb[2.2] ⟷ r[1.2]
    cb[2.3] ⟷ r[1.3]
    vs[2.1] ⟷ gnd1
    vs[2.2] ⟷ gnd2
    vs[2.3] ⟷ gnd3
    r[2.1] ⟷ gnd1
    r[2.2] ⟷ gnd2
    r[2.3] ⟷ gnd3
end

imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1,:Node2, :Node3],
                                        output_pins= Any[:gnd1,:gnd2, :gnd3], omega_range = (0, 4, 1000))
writedlm("./files/imp_abc_CB.csv",  imp, ',')
writedlm("./files/w_abc_CB.csv",  omega, ',')
