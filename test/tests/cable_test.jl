net = @network begin
    vs = dc_source(voltage = 5)
    c = cable(length = 100e3, positions = [(0,1)], earth_parameters = (1,1,1),
              C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8), C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
              C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
              I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
              I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
              I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3))
    vs[1.1] ⟷ c[1.1] ⟷ Node1
    vs[2.1] ⟷ c[2.1] ⟷  gnd
end
imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1], output_pins = Any[:gnd],
        omega_range = (-1,6,10000))
bode(imp, omega = omega)
