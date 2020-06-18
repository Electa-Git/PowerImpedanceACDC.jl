using SymEngine
s = symbols("s")
net = @network begin
    vs = ac_source(V = 380, pins = 3)
    cb = crossbonded_cable(C1 = cable(length = 866.7, positions = [(0,1.9), (0.5,1.9), (1,1.9)],
             C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8, μᵣ = 1),
             C2 = Conductor(rᵢ = 60.85e-3, rₒ = 61.05e-3, ρ = 1.72e-8, μᵣ = 1),
             I1 = Insulator(rᵢ = 31.75e-3, a = 33.75e-3, b = 59.55e-3, rₒ = 60.85e-3, ϵᵣ = 2.26),
             I2 = Insulator(rᵢ = 61.05e-3, rₒ = 65.95e-3, ϵᵣ = 2.26)),
             nₛ = 3, #minor sections
             mₛ = 10) #major sections
    #z = impedance(pins = 3, z = 25e-3*s)

    vs[1.1] ⟷ cb[1.1] ⟷ Node1
    vs[1.2] ⟷ cb[1.2] ⟷ Node2
    vs[1.3] ⟷ cb[1.3] ⟷ Node3
    #t[2.1] ⟷ z[1.1]
    #t[2.2] ⟷ z[1.2]
    #t[2.3] ⟷ z[1.3]
    vs[2.1] ⟷ vs[2.2] ⟷ vs[2.3] ⟷ gndvs
    cb[2.1] ⟷ cb[2.2] ⟷ cb[2.3] ⟷ gnd
end

imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1,:Node2],
                                        output_pins= Any[:gnd,:gnd], omega_range = (0, 4.196119877, 2500))
 bode(imp, omega = omega)
#import as did last time in matlab
#see on the simulator tutorial if it is possible to report other axis data format
