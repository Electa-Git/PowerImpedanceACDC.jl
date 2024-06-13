# For debugging
# include("../src/HVDCstability.jl")
# using .HVDCstability
# using DelimitedFiles,SymEngine
# For normal operation
using DelimitedFiles,SymEngine,HVDCstability
s = symbols("s")
transmissionVoltage = 333 / sqrt(3)

@time net = @network begin

        g1 = ac_source(V = transmissionVoltage, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)

        z_ground = impedance(z = 1e-3, pins = 3, transformation = true) 

         tl1 = transformer(n = 380/380 , Lₚ = 0.0269/2, Rₚ = 0, Rₛ = 0, Lₛ = 0.0269/2, Cₜ = 7e-9, Cₛ = 12e-9, pins = 3, transformation = true)
         
         
        #  cable(length = 60e3, positions = [(0,1.9), (0.5,1.9), (1,1.9)],
        #      C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8, μᵣ = 1),
        #      I1 = Insulator(rᵢ = 31.75e-3, a = 0.002, b = 0.0013, rₒ = 60.85e-3, ϵᵣ = 2.26),
        #      C2 = Conductor(rᵢ = 60.85e-3, rₒ = 61.05e-3, ρ = 1.72e-8, μᵣ = 1),
        #      I2 = Insulator(rᵢ = 61.05e-3, rₒ = 65.95e-3, ϵᵣ = 2.26),
        #      earth_parameters = (1,1,25), transformation = true)



#= tl1 = overhead_line(length = 90e3,
                conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                                Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
                groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
                earth_parameters = (1,1,100), transformation = true) =#
 

        z_ground[2.1] ⟷ tl1[2.1]
        z_ground[2.2] ⟷ tl1[2.2]

        g1[1.1] ⟷ tl1[1.1] ⟷ Bus1d
        g1[1.2] ⟷ tl1[1.2] ⟷ Bus1q

        g1[2.1] ⟷ gndd
        g1[2.2] ⟷ gndq

        z_ground[1.1] == gndd
        z_ground[1.2] == gndq


end

@time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:g1], input_pins=Any[:Bus1d,:Bus1q], 
output_pins=Any[:gndd,:gndq], omega_range = (-2,4,20000))

writedlm("./files/imp_TL_validation.csv",  imp_ac, ',')
writedlm("./files/w_TL_validation.csv",  omega_ac, ',')
