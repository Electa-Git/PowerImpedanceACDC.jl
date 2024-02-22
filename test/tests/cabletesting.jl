# For debugging 
# include("../../src/HVDCstability.jl")
# using .HVDCstability

# For better match with sweep: In MMC.jl, VGd = Vm and VGq = 0

using DelimitedFiles
using SymEngine
using Plots
using LinearAlgebra

s = symbols("s")

#Reactive power compensation


net = @network begin

    # Zg1 = impedance(z = 1e-3 , pins = 3, transformation = false)
    Zg1 = impedance(z = 1e-3 , pins = 1, transformation = false)

    # CableDC12 = cable(length = 100e3, positions = [(-0.5,1), (0.5,1)],
    #     C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
    #     I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
    #     C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
    #     I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
    #     C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),     
    #     I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3),
    #     earth_parameters = (1,1,100), transformation = true)  
    
    # PhD thesis Thomas Roose
    CableDC12 = cable(length = 100e3, positions = [(-0.5,1), (0.5,1)],
        C1 = Conductor(rₒ = 25.13e-3, ρ = 1.72e-8),
        I1 = Insulator(rᵢ = 25.13e-3, rₒ = 48.38e-3, ϵᵣ = 2.5),
        C2 = Conductor(rᵢ = 48.38e-3, rₒ = 50.13e-3, ρ = 2.14e-7),
        I2 = Insulator(rᵢ = 50.13e-3, rₒ = 54.13e-3, ϵᵣ = 2.3),
        C3 = Conductor(rᵢ = 54.13e-3, rₒ = 59.73e-3, ρ = 1.38e-7, μᵣ = 10),     
        I3 = Insulator(rᵢ = 59.73e-3, rₒ = 64.73e-3, ϵᵣ = 2.3),
        earth_parameters = (1,1,100), transformation = true)    


    # CableAC23 = cable(length = 60e3, positions = [(0,1.33562), (-0.0575,1.4322), (0.0575,1.4322)],
    #     C1 = Conductor(rₒ = 25.45e-3, ρ = 2.63e-8, μᵣ = 1),
    #     I1 = Insulator(rᵢ = 25.45e-3, rₒ = 50.65e-3, a = 25.45e-3 + 2e-3, b = 50.65e-3 - 1.5e-3, ϵᵣ = 2.3),
    #     C2 = Conductor(rᵢ = 50.65e-3, rₒ = 52.35e-3, ρ = 22e-8, μᵣ = 1),
    #     I2 = Insulator(rᵢ = 52.35e-3, rₒ = 55.75e-3, ϵᵣ = 100),
    #     earth_parameters = (1,1,1), transformation = false)

    # CableAC23 = cable(length = 60e3, positions = [(0,1.317150236), (-0.07175,1.441424882), (0.07175,1.441424882)],
    #     C1 = Conductor(rₒ = 25.45e-3, ρ = 2.63e-8, μᵣ = 1),
    #     I1 = Insulator(rᵢ = 25.45e-3, rₒ = 50.65e-3, a = 25.45e-3 + 2e-3, b = 50.65e-3 - 1.5e-3, ϵᵣ = 2.3),
    #     C2 = Conductor(rᵢ = 50.65e-3, rₒ = 52.35e-3, ρ = 22e-8, μᵣ = 1),
    #     I2 = Insulator(rᵢ = 52.35e-3, rₒ = 55.75e-3, ϵᵣ = 2.3,),
    #     C3 = Conductor(rᵢ = 55.75e-3, rₒ = 66.55e-3, ρ = 18e-8, μᵣ = 10),
    #     I3 = Insulator(rᵢ = 66.55e-3, rₒ = 71.75e-3, ϵᵣ = 2.3), 
    #     earth_parameters = (1,1,1), transformation = false)

    # g2 = ac_source(V = Vm, P = -Powf, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = false)

    g2 = dc_source(V = Vm, P = -Powf, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 1, transformation = false)


    # g2[1.1] ⟷ CableAC23[1.1] == TestNodeD
    # g2[1.2] ⟷ CableAC23[1.2] == TestNodeQ

    # g2[2.1] ⟷ gndd
    # g2[2.2] ⟷ gndq

    # CableAC23[2.1] ⟷ Zg1[1.1] == gndd
    # CableAC23[2.2] ⟷ Zg1[1.2] == gndq

    # Zg1[2.1] ⟷ gndd
    # Zg1[2.2] ⟷ gndq

    # g2[1.1] ⟷ CableAC23[1.1] == TestNodeA
    # g2[1.2] ⟷ CableAC23[1.2] == TestNodeB
    # g2[1.3] ⟷ CableAC23[1.3] == TestNodeC

    # g2[2.1] ⟷ gnda
    # g2[2.2] ⟷ gndb
    # g2[2.3] ⟷ gndc

    # CableAC23[2.1] ⟷ Zg1[1.1] 
    # CableAC23[2.2] ⟷ Zg1[1.2] 
    # CableAC23[2.3] ⟷ Zg1[1.3] 

    # Zg1[2.1] ⟷ gnda
    # Zg1[2.2] ⟷ gndb
    # Zg1[2.3] ⟷ gndc

    g2[1.1] ⟷ CableDC12[1.1] == TestNodeA

    g2[2.1] ⟷ gnda

    CableDC12[2.1] ⟷ Zg1[1.1] 


    Zg1[2.1] ⟷ gnda
    # Zg1[2.2] ⟷ gndb
    # Zg1[2.3] ⟷ gndc


    
end

# @time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:g2], input_pins=Any[:TestNodeD,:TestNodeQ], 
# output_pins=Any[:gndd,:gndq], omega_range = (-2,4,2000))
# @time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:g2], input_pins=Any[:TestNodeA,:TestNodeB,:TestNodeC], 
# output_pins=Any[:gnda,:gndb,:gndc], omega_range = (-2,4,2000))
@time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:g2], input_pins=Any[:TestNodeA], 
output_pins=Any[:gnda], omega_range = (-2,4,2000))


writedlm("./files/imp_cable.csv",  imp_ac, ',')
writedlm("./files/w_cable.csv",  omega_ac, ',')