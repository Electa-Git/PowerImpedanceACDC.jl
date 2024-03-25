# For debugging 
# include("../../src/HVDCstability.jl")
# using .HVDCstability

# For better match with sweep: In MMC.jl, VGd = Vm and VGq = 0

using DelimitedFiles
using SymEngine
using Plots
using LinearAlgebra

s = symbols("s")

Ptransfer = 1000
Q1 = 0
Q2 = 0
Q3 = 0
Q4 = 0

Vm_offshore = 273 / sqrt(3) #Vln,rms
Vm_onshore = 380 / sqrt(3) #Vln,rms
Vdcu = 525
Vdcl = -525
Vdc = Vdcu - Vdcl
Ztrafo_base_offshore = 273^2/1250
Ltrafo_offshore = 0.06 * Ztrafo_base_offshore /2/pi/50
Rtrafo_offshore = 0.005 * Ztrafo_base_offshore

Ztrafo_base_onshore = 273^2/1100
Ltrafo_onshore = 0.06 * Ztrafo_base_onshore /2/pi/50
Rtrafo_onshore = 0.005 * Ztrafo_base_onshore


net = @network begin

    Zg1 = impedance(z = 0.1 + 0.005*s, pins = 3, transformation = true)
    Zg2 = impedance(z = 0.1 + 0.005*s, pins = 3, transformation = true)

    G1 = ac_source(pins = 3, V = Vm_onshore, transformation = true)
    G2 = ac_source(pins = 3, V = Vm_offshore, transformation = true)

    # Onshore converters
    MMC1 = mmc_bi(Vᵈᶜᵘ = Vdcu, Vᵈᶜˡ = Vdcl, vDCbase = 525, Sbase = 1000, vACbase_LL_RMS = 273, turnsRatio = 273/380, Vₘ = Vm_onshore, Lᵣ = Ltrafo_onshore, Rᵣ = Rtrafo_onshore,
        P = Ptransfer/2, Q = Q1, 
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        dc = PI_control(Kₚ = 5, Kᵢ = 15),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )
    
    MMC2 = mmc_bi(Vᵈᶜᵘ = Vdcu, Vᵈᶜˡ = Vdcl, vDCbase = 525, Sbase = 1000, vACbase_LL_RMS = 273, turnsRatio = 273/380, Vₘ = Vm_onshore, Lᵣ = Ltrafo_onshore, Rᵣ = Rtrafo_onshore,
        P = Ptransfer/2, Q = Q2, 
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        dc = PI_control(Kₚ = 5, Kᵢ = 15),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )

    # Offshore converters
    MMC3 = mmc_bi(Vᵈᶜᵘ = Vdcu, Vᵈᶜˡ = Vdcl, vDCbase = 525, Sbase = 1000, vACbase_LL_RMS = 273, turnsRatio = 273/220, Vₘ = Vm_offshore, Lᵣ = Ltrafo_offshore, Rᵣ = Rtrafo_offshore,
        P = Ptransfer/2, Q = Q3, 
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        p = PI_control(Kₚ = 0.1, Kᵢ = 31.4159),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )

    MMC4 = mmc_bi(Vᵈᶜᵘ = Vdcu, Vᵈᶜˡ = Vdcl, vDCbase = 525, Sbase = 1000, vACbase_LL_RMS = 273, turnsRatio = 273/220, Vₘ = Vm_offshore, Lᵣ = Ltrafo_offshore, Rᵣ = Rtrafo_offshore,
        P = Ptransfer/2, Q = Q4, 
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        p = PI_control(Kₚ = 0.1, Kᵢ = 31.4159),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )


    CableDC_upper = cable(length = 60e3, positions = [(0,1)],
        C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
        I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
        C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
        I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
        C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),     
        I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = false)      

    CableDC_lower = cable(length = 60e3, positions = [(0,1)],
        C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
        I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
        C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
        I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
        C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),     
        I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = false)
        
    CableDC_MR = cable(length = 60e3, positions = [(0,1)],
        C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
        I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
        C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
        I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
        C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),     
        I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = false) 
    
    G1[2.1] ⟷ gndD
    G1[2.2] ⟷ gndQ
    G2[2.1] ⟷ gndD
    G2[2.2] ⟷ gndQ

    G1[1.1] ⟷ Zg1[1.1] ⟷ OnshoreD1
    G1[1.2] ⟷ Zg1[1.2] ⟷ OnshoreQ1

    Zg1[2.1] ⟷ MMC1[2.1] ⟷ MMC2[2.1] ⟷ OnshoreD2
    Zg1[2.2] ⟷ MMC1[2.2] ⟷ MMC2[2.2] ⟷ OnshoreQ2

    MMC1[1.1] ⟷ CableDC_upper[1.1] ⟷ OnshoreDC1
    MMC1[1.2] ⟷ MMC2[1.1] ⟷ CableDC_MR[1.1] ⟷ OnshoreDC2
    MMC2[1.2] ⟷ CableDC_lower[1.1] ⟷ OnshoreDC3

    MMC3[1.1] ⟷ CableDC_upper[2.1] ⟷ OffshoreDC1
    MMC3[1.2] ⟷ MMC4[1.1] ⟷ CableDC_MR[2.1] ⟷ OffshoreDC2
    MMC4[1.2] ⟷ CableDC_lower[2.1] ⟷ OffshoreDC3

    G2[1.1] ⟷ Zg2[1.1] ⟷ OffshoreD1
    G2[1.2] ⟷ Zg2[1.2] ⟷ OffshoreQ1

    Zg2[2.1] ⟷ MMC3[2.1] ⟷ MMC4[2.1] ⟷ OffshoreD2
    Zg2[2.2] ⟷ MMC3[2.2] ⟷ MMC4[2.2] ⟷ OffshoreQ2
    
end

##### ----- STABILITY ANALYSIS WITH DETERMINE_IMPEDANCE ----- #####

# DC-side
# Z_MMC_DC1, omega = determine_impedance(net, elim_elements = [:CableDC12], input_pins = Any[:NodeDC1], output_pins = Any[:gndD], omega_range = (-3, 4, 1000))
# Z_BUS_DC1, omega = determine_impedance(net, elim_elements = [:MMC1], input_pins = Any[:NodeDC1], output_pins = Any[:gndD], omega_range = (-3, 4, 1000))
# L_DC1 = Z_BUS_DC1 ./ Z_MMC_DC1

# nyquist_L_DC1 = nyquistplot(L_DC1, omega, zoom = "no", SM = "VM") 

# AC-side
Z_MMC_AC1, omega = determine_impedance(net, elim_elements = [:Zg1], input_pins = Any[:OnshoreD2, :OnshoreQ2], output_pins = Any[:gndD, :gndQ], omega_range = (-3, 4, 1000))
Z_BUS_AC1, omega = determine_impedance(net, elim_elements = [:MMC1,:MMC2], input_pins = Any[:OnshoreD2, :OnshoreQ2], output_pins = Any[:gndD, :gndQ], omega_range = (-3, 4, 1000))
L_AC1 = Z_BUS_AC1 ./ Z_MMC_AC1

nyquist_L_AC1 = nyquistplot(L_AC1, omega, zoom = "yes", SM = "PM") 



