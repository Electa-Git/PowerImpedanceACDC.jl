# include("../../src/HVDCstability.jl")
# using HVDCstability
using DelimitedFiles
using SymEngine
## cross-bonding cable test
s = symbols(:s)
#ω=2*π*50
#TODO: Not updated after a change in the way semiconducting layers are modeled. Insulator data not accurate.

net = @network begin
    ac = ac_source(V = 380, pins = 3)
    cb = crossbonded_cable(section = cable(length = 834, positions = [(-0.5,1.9), (0,1.9), (0.5,1.9)],
             C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8, μᵣ = 1),
             C2 = Conductor(rᵢ = 60.58e-3, rₒ = 60.78e-3, ρ = 1.72e-8, μᵣ = 1),
             I1 = Insulator(rᵢ = 31.75e-3, a = 31.95e-3, b = 60.45e-3, rₒ = 60.58e-3, ϵᵣ = 2.26),
             I2 = Insulator(rᵢ = 60.78e-3, rₒ = 65.68e-3, ϵᵣ = 2.26)),
             nₛ = 3, #minor sections
             mₛ = 1, #major sections
             Zₛ = 1e-6,
             Zᶜᵇ = 1e-3*s)
    #r = impedance(z = 1e-3+s*1e-30, pins = 3) #why is not possible to have z= Real?
    ac[1.1] ⟷ cb[1.1] ⟷ Node1
    ac[1.2] ⟷ cb[1.2] ⟷ Node2
    ac[1.3] ⟷ cb[1.3] ⟷ Node3
    ac[2.1] ⟷ gnd1
    ac[2.2] ⟷ gnd2
    ac[2.3] ⟷ gnd3
    # ac[2.1] ⟷ ac[2.2] ⟷ ac[2.3] ⟷ gnd
    #cb[2.1] ⟷ r[1.1]
    #cb[2.2] ⟷ r[1.2]
    #cb[2.3] ⟷ r[1.3]
    cb[2.1] ⟷ gnd1
    cb[2.2] ⟷ gnd2
    cb[2.3] ⟷ gnd3
end
imp, omega = determine_impedance(net, elim_elements = [:ac], input_pins = Any[:Node1, :Node2, :Node3],
                            output_pins = Any[:gnd1, :gnd2, :gnd3], omega_range = (1,5,1000))
#save_data(net.elements[:Cable],"cb_trial", omega_range = (5, 5, 1), scale=:loglin)
bode(imp, omega = omega, axis_type = :loglin)
# writedlm( "imp_cable_CB_semi_1_MS_3_ns.csv", imp,',')
# writedlm( "omega_cable_CB_semi_1_MS_3_ns.csv", omega, ',')
