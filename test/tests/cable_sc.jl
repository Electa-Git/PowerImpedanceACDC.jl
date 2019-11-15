net = @network begin
    vs = dc_source(voltage = 380e3)
    c = cable(length = 100e3, positions = [(0,2)], earth_parameters = (1,1,1),
              type = :aerial,
              C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8),
              # SC = Conductor(rᵢ = 60.85e-3, rₒ = 61.05e-3, ρ = 1.72e-8),
              C2 = Conductor(rᵢ = 60.85e-3, rₒ = 65.95e-3, ρ = 1.72e-8),
              I1 = Insulator(rᵢ = 31.75e-3, a = 33.75e-3, b = 59.55e-3, rₒ = 60.85e-3, ϵᵣ = 2.26),
              I2 = Insulator(rᵢ = 65.95e-3, rₒ = 68e-3, ϵᵣ = 2.26))
    vs[1.1] ⟷ c[1.1] ⟷ Node1
    vs[2.1] ⟷ c[2.1] ⟷  gnd
end
imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1], output_pins = Any[:gnd],
        omega_range = (-1,6,10000))
bode(imp, omega = omega)
