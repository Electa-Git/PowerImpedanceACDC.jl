net = @network begin
    c = cable(length = 10e3, positions = [(0,1)], earth_parameters = (1,1,25),
              type = :aerial,
              C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8, μᵣ = 1),
              C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
              C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
              I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
              I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
              I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3))
    vs[1.1] ⟷ c[1.1] ⟷ Node1
    vs[2.1] ⟷ c[2.1] ⟷  gnd
