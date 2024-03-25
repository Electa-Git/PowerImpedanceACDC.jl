using SymEngine,DelimitedFiles
# For debugging 
# include("../../src/HVDCstability.jl")
# using .HVDCstability
s = symbols("s")
net = @network begin
    vs = ac_source(V = 380e3, pins = 3, transformation = false)
    # Cable used by Willem in his TPWRD paper
    cb = crossbonded_cable(section = cable(length = 834, positions = [(0,1.9), (0.5,1.9), (1,1.9)],
            C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8, μᵣ = 1),
            I1 = Insulator(rᵢ = 31.75e-3, a = 0.2e-3, b = 0.13e-3, rₒ = 60.58e-3, ϵᵣ = 2.26),
            C2 = Conductor(rᵢ = 60.58e-3, rₒ = 60.78e-3, ρ = 1.72e-8, μᵣ = 1),
            I2 = Insulator(rᵢ = 60.78e-3, rₒ = 65.68e-3, ϵᵣ = 2.26),
            earth_parameters = (1,1,25)),
             nₛ = 3, #minor sections
             mₛ = 2, #major sections
             Zᶜᵇ = 1e-6s,  # cross-bonding impedance, zero corresponds to a direct connection
             Zₛ = 1,  # sheath grounding impedance
             transformation = false) 
    # Test cable
    # cb = crossbonded_cable(section = cable(length = 10e3, positions = [(0,1.9), (0.5,1.9), (1,1.9)],
    #         C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8, μᵣ = 1),
    #         # I1 = Insulator(rᵢ = 31.75e-3, rₒ = 60.85e-3, ϵᵣ = 2.26),
    #         I1 = Insulator(rᵢ = 31.75e-3, a = 0.002, b = 0.0013, rₒ = 60.85e-3, ϵᵣ = 2.26),
    #         # I1 = Insulator(rᵢ = 31.75e-3, a = 31.75e-3, b = 60.85e-3, rₒ = 60.85e-3, ϵᵣ = 2.58908112),
    #         C2 = Conductor(rᵢ = 60.85e-3, rₒ = 61.05e-3, ρ = 1.72e-8, μᵣ = 1),
    #         I2 = Insulator(rᵢ = 61.05e-3, rₒ = 65.95e-3, ϵᵣ = 2.26),
    #         earth_parameters = (1,1,25)),
    #          nₛ = 3, #minor sections
    #          mₛ = 1, #major sections
    #          Zᶜᵇ = 1e-3+1e-1s,  # cross-bonding impedance, zero corresponds to a direct connection
    #          Zₛ = 1e-6s,  # sheath grounding impedance
    #          transformation = false) 
    # cb = cable(length = 30e3, positions = [(0,1.9), (0.5,1.9), (1,1.9)],
    #          C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8, μᵣ = 1),
    #          I1 = Insulator(rᵢ = 31.75e-3, a = 0.002, b = 0.0013, rₒ = 60.85e-3, ϵᵣ = 2.26),
    #          C2 = Conductor(rᵢ = 60.85e-3, rₒ = 61.05e-3, ρ = 1.72e-8, μᵣ = 1),
    #          I2 = Insulator(rᵢ = 61.05e-3, rₒ = 65.95e-3, ϵᵣ = 2.26),
    #          earth_parameters = (1,1,25))

    r = impedance(pins = 3, z = 1, transformation = false)

    vs[1.1] ⟷ cb[1.1] ⟷ Node1
    vs[1.2] ⟷ cb[1.2] ⟷ Node2
    vs[1.3] ⟷ cb[1.3] ⟷ Node3
    
    vs[2.1] ⟷ gnd1
    vs[2.2] ⟷ gnd2
    vs[2.3] ⟷ gnd3
    cb[2.1] ⟷ r[1.1] ⟷ Node4
    cb[2.2] ⟷ r[1.2] ⟷ Node5
    cb[2.3] ⟷ r[1.3] ⟷ Node6
    # r[2.1] ⟷ gnd1
    # r[2.2] ⟷ gnd2
    # r[2.3] ⟷ gnd3
    r[2.1] ⟷ gnd1
    r[2.2] ⟷ gnd2
    r[2.3] ⟷ gnd3
    # cb[2.1] ⟷ gnd1
    # cb[2.2] ⟷ gnd2
    # cb[2.3] ⟷ gnd3
    # vs[1.1] ⟷ cb[1.1] ⟷ Node1
    # vs[1.2] ⟷ cb[1.2] ⟷ Node2
    # cb[2.1] ⟷ r[1.1] ⟷ Node4
    # cb[2.2] ⟷ r[1.2] ⟷ Node5
    # vs[2.1] ⟷ gnd1
    # vs[2.2] ⟷ gnd2
    # r[2.1] ⟷ gnd1
    # r[2.2] ⟷ gnd2
end

imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1,:Node2, :Node3],
                                        output_pins= Any[:gnd1,:gnd2, :gnd3], omega_range = (0, 4, 1500))
writedlm("./files/imp_abc_CB.csv",  imp, ',')
writedlm("./files/w_abc_CB.csv",  omega, ',')

# imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1,:Node2, ],
#                                         output_pins= Any[:gnd1,:gnd2], omega_range = (0, 4, 1000))
# writedlm("./files/imp_dq_CB.csv",  imp, ',')
# writedlm("./files/w_dq_CB.csv",  omega, ',')


# min_ω = 0
# max_ω = 4
# n_ω = 1000
# omegas= 2*pi* 10 .^range(min_ω, max_ω, length= n_ω) 

# Y_BUS =  HVDCstability.make_y_matrix(net,elim_elements = [:vs],input_pins = Any[:Node1,:Node2,:Node3,:Node4,:Node5,:Node6], 
# omega_range = (min_ω, max_ω, n_ω))
# no_elim = [i for i in 1:3]
# Y_BUS_kron=[]
# for i in eachindex(Y_BUS)
#     push!(Y_BUS_kron, HVDCstability.kron(Y_BUS[i],no_elim))
# end

# ABCD_cable =[]
# Y_cable =[]
# no_elim = [1,3,5]
# for i in eachindex(omegas)
#     ABCD_cable = HVDCstability.kron_abcd(HVDCstability.eval_abcd(net.elements[:cb].element_value,omegas[i]*1im), 0, no_elim)
#     # display(ABCD_cable)
#     push!(Y_cable,HVDCstability.abcd_to_y(ABCD_cable))
# end