# For debugging 
# include("../../src/HVDCstability.jl")
# using .HVDCstability

# For better match with sweep: In MMC.jl, VGd = Vm and VGq = 0

using DelimitedFiles
using SymEngine
using Plots
using LinearAlgebra

s = symbols("s")

Powf = 700
Qowf = -150
Pmmc1 = (1/3)*Powf
Pmmc2 = -(2/3)*Powf
Qref = 0
Vm = 220 / sqrt(3) #Vln,rms
Vdc = 640
Ztrafo_base = 220^2/1000
Ltrafo = 0.18 * Ztrafo_base /2/pi/50
Rtrafo = 0.005 * Ztrafo_base

#Reactive power compensation
Qcomp = 50 #MVar
Lcomp = (3*Vm^2)/(2*pi*50*Qcomp)

multiTerminalAnalysis = true
singleTerminalAnalysis = false

net = @network begin

    Zg3 = impedance(z = 0.1 + 0.005*s, pins = 3, transformation = false)

    G3 = ac_source(pins = 3, V = Vm, transformation = false)

    G2 = ac_source(pins = 3, V = Vm, transformation = false)
    
    CableAC23 = cable(length = 60e3, positions = [(0,1.33562), (-0.0575,1.4322), (0.0575,1.4322)],
        C1 = Conductor(rₒ = 25.45e-3, ρ = 2.63e-8, μᵣ = 1),
        I1 = Insulator(rᵢ = 25.45e-3, rₒ = 50.65e-3, a = 2e-3, b = 1.5e-3, ϵᵣ = 2.3),
        C2 = Conductor(rᵢ = 50.65e-3, rₒ = 52.35e-3, ρ = 22e-8, μᵣ = 1),
        I2 = Insulator(rᵢ = 52.35e-3, rₒ = 55.75e-3, ϵᵣ = 100),
        earth_parameters = (1,1,1), transformation = false)

    CableAC2g2 = cable(length = 10e3, positions = [(0,1.33562), (-0.0575,1.4322), (0.0575,1.4322)],
        C1 = Conductor(rₒ = 25.45e-3, ρ = 2.63e-8, μᵣ = 1),
        I1 = Insulator(rᵢ = 25.45e-3, rₒ = 50.65e-3, a = 2e-3, b = 1.5e-3, ϵᵣ = 2.3),
        C2 = Conductor(rᵢ = 50.65e-3, rₒ = 52.35e-3, ρ = 22e-8, μᵣ = 1),
        I2 = Insulator(rᵢ = 52.35e-3, rₒ = 55.75e-3, ϵᵣ = 100),
        earth_parameters = (1,1,1), transformation = false)

    OHLAC3g3 = overhead_line(length = 50e3,
        conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                        Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        earth_parameters = (1,1,100), transformation = false)

    G3[2.1] == gndA
    G3[2.2] == gndB
    G3[2.3] == gndC

    G2[1.1] == gndA
    G2[1.2] == gndB
    G2[1.3] == gndC

    

    CableAC23[2.1] ⟷ CableAC2g2[2.1]  ⟷ NodeA2
    CableAC23[2.2] ⟷ CableAC2g2[2.2]  ⟷ NodeB2
    CableAC23[2.3] ⟷ CableAC2g2[2.3]  ⟷ NodeC2

    CableAC23[1.1] ⟷ OHLAC3g3[1.1]  ⟷ NodeA3
    CableAC23[1.2] ⟷ OHLAC3g3[1.2]  ⟷ NodeB3
    CableAC23[1.3] ⟷ OHLAC3g3[1.3]  ⟷ NodeC3

    G2[2.1] ⟷ CableAC2g2[1.1] ⟷ NodeAg2
    G2[2.2] ⟷ CableAC2g2[1.2] ⟷ NodeBg2
    G2[2.3] ⟷ CableAC2g2[1.3] ⟷ NodeCg2

    # With external grid impedance
    G3[1.1] ⟷ Zg3[1.1] ⟷ NodeA4
    G3[1.2] ⟷ Zg3[1.2] ⟷ NodeB4
    G3[1.3] ⟷ Zg3[1.3] ⟷ NodeC4

    Zg3[2.1] ⟷ OHLAC3g3[2.1] ⟷ NodeAg3
    Zg3[2.2] ⟷ OHLAC3g3[2.2] ⟷ NodeBg3
    Zg3[2.3] ⟷ OHLAC3g3[2.3] ⟷ NodeCg3

    
end

    omega_min = 0
    omega_max = 4
    omegaₙ = 1000
   

    # First get the full matrix, then apply Kron reduction
    Y_BUS_AC2_new =  HVDCstability.make_y_matrix(net,elim_elements = [:Zg3],input_pins = Any[:NodeAg3,:NodeBg3,:NodeCg3,:NodeAg2,:NodeBg2,:NodeCg2,
    :NodeA3,:NodeB3,:NodeC3,:NodeA2,:NodeB2,:NodeC2], omega_range = (omega_min, omega_max, omegaₙ))
    no_elim = [i for i in 1:6]
    Y_BUS_AC2_new_kron=[]
    for i in eachindex(Y_BUS_AC2_new)
        push!(Y_BUS_AC2_new_kron, HVDCstability.kron(Y_BUS_AC2_new[i],no_elim))
    end

    Y_BUS_AC2 = []
    
    omega= 2*pi* 10 .^range(omega_min, omega_max, length= omegaₙ) 
    for i in eachindex(omega)

        ABCD23 = eval_abcd(net.elements[:CableAC23].element_value, 1im*omega[i])
        ABCD2g2 = eval_abcd(net.elements[:CableAC2g2].element_value, 1im*omega[i])
        ABCD3g3 = eval_abcd(net.elements[:OHLAC3g3].element_value, 1im*omega[i])

        ABCDeq = Matrix{Complex}(ABCD3g3 * (ABCD23 * ABCD2g2))
        
        push!(Y_BUS_AC2, HVDCstability.abcd_to_y(ABCDeq))
    end

    
