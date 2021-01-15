using SymEngine

# semiconducting layers test
net = @network begin
    #ac = "3-phase ac voltage source with V amplitude 380kV
    ac = ac_source(V = 380, pins = 3)
    #c = cable with length 26km
    cb = crossbonded_cable(section = cable(length = 886.7, positions = [(-0.5,1.9), (0,1.9), (0.5,1.9)],
              C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8), #inner radius[m], resistivity [Ωm]
              C2 = Conductor(rᵢ = 57.55e-3, rₒ = 57.75e-3, ρ = 1.72e-8), # inner and outer radius [m], resistivity [Ωm]
              I1 = Insulator(rᵢ = 31.75e-3, rₒ = 57.55e-3, ϵᵣ = 2.26), #inner radius, inner semiconductor outer radius, outer semiconductor inner radius, outer radius, relative permittivity
              I2 = Insulator(rᵢ = 57.75e-3, rₒ = 62.65e-3, ϵᵣ = 2.26)),
              nₛ = 1, #minor sections
              mₛ = 1) #major sections
    ac[1.1] ⟷ cb[1.1] ⟷ Node1
    ac[1.2] ⟷ cb[1.2] ⟷ Node2
    ac[1.3] ⟷ cb[1.3] ⟷ Node3
    ac[2.1] ⟷ ac[2.2] ⟷ ac[2.3] ⟷ gnd
    cb[2.1] ⟷ gnd1
    cb[2.2] ⟷ gnd2
    cb[2.3] ⟷ gnd3
end

imp, omega = determine_impedance(net, elim_elements = [:ac], input_pins = Any[:Node1,:Node2,:Node3],
                            output_pins = Any[:gnd1,:gnd2,:gnd3], omega_range = (0.1,5,500))
bode(imp, omega = omega, axis_type = :loglog)
