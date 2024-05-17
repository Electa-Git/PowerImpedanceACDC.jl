# Test file to have fun :D
# OWF AC connected
# Single node stability analysis


#= include("../../src/HVDCstability.jl")
using .HVDCstability 
 =#

using DelimitedFiles
using SymEngine
using Plots
using LinearAlgebra



s = symbols("s")
Pmmc = 700
Qmmc = 0
Vm = 380 / sqrt(3) #Vln,rms
Vdc = 640

Ztrafo_base = 380^2/1000
Ltrafo = 0.18 * Ztrafo_base /2/pi/50
Rtrafo = 0.005 * Ztrafo_base


grid=@network begin
    


    voltageBase = Vm # Apparently needed for correct per-unit calculation of the power flow

    G3 = ac_source(pins = 3, V = Vm, transformation = true)

    Zg3 = impedance(z = 0.28 + 0.009*s, pins = 3, transformation = true)

    OHLg3MMC1 = overhead_line(length = 50e3,
        conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                        Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        earth_parameters = (1,1,100), transformation = true)
    
    MMC1 = mmc(Vᵈᶜ = Vdc, vDCbase = 640, Sbase = 1000, vACbase_LL_RMS = 333, turnsRatio = 333/380, Vₘ = Vm, Lᵣ = Ltrafo, Rᵣ = Rtrafo,
        P_max = 1500, P_min = -1500, P = Pmmc, Q = Qmmc, Q_max = 1000, Q_min = -1000,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        p = PI_control(Kₚ = 0.1, Kᵢ = 31.4159),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159), 
        #p = VSE(H= 1,K_d=0,K_ω=0,n_f=0,ω_f=0)
        ) 



    CableDC12 = cable(length = 60e3, positions = [(-0.5,1), (0.5,1)],
        C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
        I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
        C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
        I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
        C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),     
        I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)       



    MMC2 = mmc(Vᵈᶜ = Vdc, vDCbase = 640, Sbase = 1000, vACbase_LL_RMS = 333, turnsRatio = 333/380, Vₘ = Vm, Lᵣ = Ltrafo, Rᵣ = Rtrafo,
        P_max = 1500, P_min = -1500, P = -Pmmc, Q = Qmmc, Q_max = 1000, Q_min = -1000,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        dc = PI_control(Kₚ = 5, Kᵢ = 15),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )    
   
    OHLg2MMC2 = overhead_line(length = 50e3,
        conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                        Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        earth_parameters = (1,1,100), transformation = true)        


    Zg2 = impedance(z = 0.28 + 0.009*s, pins = 3, transformation = true)

    G2 = ac_source(pins = 3, V = Vm, transformation = true)

        
    G3[2.1] ⟷ gndD
    G3[2.2] ⟷ gndQ
   
    G2[2.1] ⟷ gndD
    G2[2.2] ⟷ gndQ

    G3[1.1] ⟷ Zg3[1.1] ⟷ NodeG3d
    G3[1.2] ⟷ Zg3[1.2] ⟷ NodeG3q

    Zg3[2.1] ⟷ OHLg3MMC1[2.1] ⟷ Nodeg3d
    Zg3[2.2] ⟷ OHLg3MMC1[2.2] ⟷ Nodeg3q

    OHLg3MMC1[1.1] == MMC1[2.1] ⟷ NodeMMC1d
    OHLg3MMC1[1.2] == MMC1[2.2] ⟷ NodeMMC1q
    
    MMC1[1.1]== CableDC12[1.1] ⟷ NodeDC1
    MMC2[1.1]== CableDC12[2.1] ⟷ NodeDC2

    OHLg2MMC2[1.1] == MMC2[2.1] ⟷ NodeMMC2d
    OHLg2MMC2[1.2] == MMC2[2.2] ⟷ NodeMMC2q


    Zg2[2.1] ⟷ OHLg2MMC2[2.1] ⟷ Nodeg2d
    Zg2[2.2] ⟷ OHLg2MMC2[2.2] ⟷ Nodeg2q

    G2[1.1] ⟷ Zg2[1.1] ⟷ NodeG2d
    G2[1.2] ⟷ Zg2[1.2] ⟷ NodeG2q



end

omega_Y_MMC1 = collect(range(2*pi*0.01, stop=2*pi*1000, step=1))



Y_MMC1 = []



for i in 1:length(omega_Y_MMC1)
    Y1 = eval_abcd(grid.elements[:MMC1].element_value, 1im*omega_Y_MMC1[i]) 
    push!(Y_MMC1, [transpose(Y1[1,:]); transpose(-Y1[2,:]); transpose(-Y1[3,:])]) # Keep sign of Ydc, swapping sign of Yacs
    # push!(Y_dc2, Y2[1,1]) 
end


writedlm("./files/Y_MMC1.csv",  Y_MMC1, ',')
writedlm("./files/omega.csv",  omega_Y_MMC1, ',')
