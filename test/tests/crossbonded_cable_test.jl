using SymEngine
s = symbols("s")
#TODO: Not updated after a change in the way semiconducting layers are modeled. Insulator data not accurate.

net = @network begin
    vs = ac_source(V = 380, pins = 3)
    cb = crossbonded_cable(section = cable(length = 866.7, positions = [(0,1.9), (0.5,1.9), (1,1.9)],
             C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8, μᵣ = 1),
             C2 = Conductor(rᵢ = 60.85e-3, rₒ = 61.05e-3, ρ = 1.72e-8, μᵣ = 1),
             I1 = Insulator(rᵢ = 31.75e-3, a = 33.75e-3, b = 59.55e-3, rₒ = 60.85e-3, ϵᵣ = 2.26),
             I2 = Insulator(rᵢ = 61.05e-3, rₒ = 65.95e-3, ϵᵣ = 2.26)),
             nₛ = 3, #minor sections
             mₛ = 1, 
             Zᶜᵇ = 1e-9 + 1e-9*s) #major sections
    r = impedance(pins = 3, z = 1e-3 + 1e-9*s)

    vs[1.1] ⟷ cb[1.1] ⟷ Node1
    vs[1.2] ⟷ cb[1.2] ⟷ Node2
    vs[1.3] ⟷ cb[1.3] ⟷ Node3
    cb[2.1] ⟷ r[1.1]
    cb[2.2] ⟷ r[1.2]
    cb[2.3] ⟷ r[1.3]
    vs[2.1] ⟷ vs[2.2] ⟷ vs[2.3] ⟷ gnd
    #cb[2.1] ⟷ cb[2.2] ⟷ cb[2.3] ⟷ gnd
    r[2.1] ⟷ gnd1
    r[2.2] ⟷ gnd2
    r[2.3] ⟷ gnd3
end

imp, omega = determine_impedance(net, elim_elements = [:as], input_pins = Any[:Node1,:Node2, :Node3],
                                        output_pins= Any[:gnd1,:gnd2, :gnd3], omega_range = (0.1, 5, 1000))
bode(imp, omega = omega, axis_type = :loglog)
#import as did last time in matlab
#see on the simulator tutorial if it is possible to report other axis data format
