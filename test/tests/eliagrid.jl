# include("../../src/HVDCstability.jl")
# using .HVDCstability
using SymEngine

#Notations:
#Node notation: Node xy -> x number of the node from Deliv. 1 scheme, y phase. e.g. Node 31 -> Node 3, Phase 1
#Cable notation: Cable cdan -> c cable, d departing substation, a arriving substation, n number of the cable (if more than 1 needed)
#e.g. cMG1 -> c Cable, m MAERL380 substation, g GEZEL380 substation, 1 cable number 1.
#Notation connections ->  for an ac source: ac[s.p] where s side and p phase. e.g. ac[1.2] ac source: side=1 phase=2.

# cross-bonded cables from MAERL380 to GEZEL380-> generic notation: cMGn where n: cable number
s = symbols(:s)
#TODO: Not updated after a change in the way semiconducting layers are modeled. Insulator data not accurate.

net = @network begin
    # Equivalent Grid FR_EQ380 380kV
    acFR_EQ380 = ac_source(V = 380, pins = 3)
    zeqFR = impedance(z= s*0.05, pins=3)

    # Equivalent Grid MERCA380
    acMERCA380 = ac_source(V = 380, pins = 3)
    zeqMERCA = impedance(z = s*0.05, pins = 3)

    # Equivalent Reactor at GEZEL380 Q=234.67Mvar
    reactorGE = impedance(z = s*1.9587, pins = 3)

    # Equivalent Reactor at STEVN380  Q=117.33Mvar
    reactorSTEVN = impedance(z = s*3.917439, pins = 3)

    #Autotransformer20 STEVN220 1
    atST1 = autotransformer(S=600, Xₕₗ = 0.15, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #Autotransformer21 STEVN220 2
    atST2 = autotransformer(S=600, Xₕₗ = 0.15, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #Autotransformer22 STEVN220 3
    atST3 = autotransformer(S=600, Xₕₗ = 0.15, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #Autotransformer23 STEVN220 4
    atST4 = autotransformer(S=600, Xₕₗ = 0.15, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #Autotransformer24 STEVN150 1
    atST5 = autotransformer(S=555, Xₕₗ = 0.23, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #Autotransformer25 STEVN150 2
    atST6 = autotransformer(S=555, Xₕₗ = 0.23, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #Autotransformer16 GEZEL220 1
    atGE1 = autotransformer(S=600, Xₕₗ = 0.15, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #Autotransformer17 GEZEL200 2
    atGE2 = autotransformer(S=600, Xₕₗ = 0.15, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #Autotransformer18 GEZEL220 3
    atGE3 = autotransformer(S=600, Xₕₗ = 0.15, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #Autotransformer19 GEZEL200 4
    atGE4 = autotransformer(S=600, Xₕₗ = 0.15, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #autotransformer26 IZGEM150 1
    atIZ1 = autotransformer(S=555, Xₕₗ = 0.23, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #Autotransformer27 IZGEM150 2
    atIZ2 = autotransformer(S=555, Xₕₗ = 0.23, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #Autotransformer28 AVLGM150 1
    atAV1 = autotransformer(S=555, Xₕₗ = 0.23, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #Autotransformer29 AVLGM150 2
    atAV2 = autotransformer(S=555, Xₕₗ = 0.23, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #Autotransformer30 EEKLN150 1
    atEE1 = autotransformer(S=555, Xₕₗ = 0.23, Xₕₜ = 0.02, Xₜₗ = 0.02, pins = 3)

    #Load AVLGM150
    loadAV1 = impedance(z = 150+1/10.61e-6/s, pins = 3) #Node14

    #Load IZGEM150
    loadIZ1 = impedance(z = 150+1/10.61e-6/s, pins = 3) #Node13

    #Load STEVN150
    loadST1 = impedance(z = 150+1/10.61e-6/s, pins = 3) #Node12

    #Load EEKLN150
    loadEE1 = impedance(z = 300+1/5.3e-6/s, pins = 3) #Node15

    #Load GEZEL220
    loadGE2 = impedance(z = 193.6+1/8.22e-6/s, pins = 3) #Node10

    #Load STEVN220
    loadST2 = impedance(z = 193.6+1/8.22e-6/s, pins = 3) #Node11

    # Undergound cable n.1 between MAERL380 and GEZEL380
    cMG1 = crossbonded_cable(C1 = cable(length = 834, positions = [(0,1.9), (0.5,1.9), (1,1.9)],
              C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8, μᵣ = 1),
              C2 = Conductor(rᵢ = 60.85e-3, rₒ = 61.05e-3, ρ = 1.72e-8, μᵣ = 1),
              I1 = Insulator(rᵢ = 31.75e-3, a = 33.75e-3, b = 59.55e-3, rₒ = 60.85e-3, ϵᵣ = 2.26),
              I2 = Insulator(rᵢ = 61.05e-3, rₒ = 65.95e-3, ϵᵣ = 2.26)),
              nₛ = 3, #minor sections
              mₛ = 4) #major sections

    # Undergound cable n.2 between MAERL380 and GEZEL380
    cMG2 = crossbonded_cable(C1 = cable(length = 834, positions = [(2.32,1.9), (2.82,1.9), (3.32,1.9)],
             C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8, μᵣ = 1),
             C2 = Conductor(rᵢ = 60.85e-3, rₒ = 61.05e-3, ρ = 1.72e-8, μᵣ = 1),
             I1 = Insulator(rᵢ = 31.75e-3, a = 33.75e-3, b = 59.55e-3, rₒ = 60.85e-3, ϵᵣ = 2.26),
             I2 = Insulator(rᵢ = 61.05e-3, rₒ = 65.95e-3, ϵᵣ = 2.26)),
             nₛ = 3, #minor sections
             mₛ = 4) #major sections

    # Undergound cable n.3 between MAERL380 and GEZEL380
    cMG3 = crossbonded_cable(C1 = cable(length = 834, positions = [(4.64,1.9), (5.14,1.9), (5.64,1.9)],
             C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8, μᵣ = 1),
             C2 = Conductor(rᵢ = 60.85e-3, rₒ = 61.05e-3, ρ = 1.72e-8, μᵣ = 1),
             I1 = Insulator(rᵢ = 31.75e-3, a = 33.75e-3, b = 59.55e-3, rₒ = 60.85e-3, ϵᵣ = 2.26),
             I2 = Insulator(rᵢ = 61.05e-3, rₒ = 65.95e-3, ϵᵣ = 2.26)),
             nₛ = 3, #minor sections
             mₛ = 4) #major sections

    # Undergound cable n.4 between MAERL380 and GEZEL380
    cMG4 = crossbonded_cable(C1 = cable(length = 834, positions = [(6.96,1.9), (7.46,1.9), (7.96,1.9)],
             C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8, μᵣ = 1),
             C2 = Conductor(rᵢ = 60.85e-3, rₒ = 61.05e-3, ρ = 1.72e-8, μᵣ = 1),
             I1 = Insulator(rᵢ = 31.75e-3, a = 33.75e-3, b = 59.55e-3, rₒ = 60.85e-3, ϵᵣ = 2.26),
             I2 = Insulator(rᵢ = 61.05e-3, rₒ = 65.95e-3, ϵᵣ = 2.26)),
             nₛ = 3, #minor sections
             mₛ = 4) #major sections

    # OHL between EEKLN380 and MAERL380
    ohlEM = overhead_line(length = 17e3,
             conductors = Conductors(organization = :vertical,
                            nᵇ = 3, Rᵈᶜ = 0.05, rᶜ = 0.016,  yᵇᶜ = 30, Δxᵇᶜ = 14, Δyᵇᶜ= 9, dˢᵃᵍ = 10),
             groundwires = Groundwires(nᵍ = 1, Rᵍᵈᶜ = 0.1, rᵍ = 0.011, Δyᵍ = 27, dᵍˢᵃᵍ = 10),
             earth_parameters = (1,1,25))

    # OHL between HORTA380 and MAERL380
    ohlHM = overhead_line(length = 29e3,
             conductors = Conductors(organization = :vertical,
                            nᵇ = 3, Rᵈᶜ = 0.05, rᶜ = 0.016,  yᵇᶜ = 30, Δxᵇᶜ = 14, Δyᵇᶜ= 9, dˢᵃᵍ = 10),
             groundwires = Groundwires(nᵍ = 1, Rᵍᵈᶜ = 0.1, rᵍ = 0.011, Δyᵍ = 27, dᵍˢᵃᵍ = 10),
             earth_parameters = (1,1,25))

    # OHL between HORTA380 and EEKLN380
    ohlHE = overhead_line(length = 12e3,
             conductors = Conductors(organization = :vertical,
                            nᵇ = 3, Rᵈᶜ = 0.05, rᶜ = 0.016,  yᵇᶜ = 30, Δxᵇᶜ = 14, Δyᵇᶜ= 9, dˢᵃᵍ = 10),
             groundwires = Groundwires(nᵍ = 1, Rᵍᵈᶜ = 0.1, rᵍ = 0.011, Δyᵍ = 27, dᵍˢᵃᵍ = 10),
             earth_parameters = (1,1,25))

    # OHLx2 between MERCA380 and HORTA380
    ohlMH = overhead_line(length = 71e3,
             conductors = Conductors(organization = :vertical,
                            nᵇ = 6, Rᵈᶜ = 0.05, rᶜ = 0.016,  yᵇᶜ = 30, Δxᵇᶜ = 14, Δyᵇᶜ= 9, dˢᵃᵍ = 10),
             groundwires = Groundwires(nᵍ = 1, Rᵍᵈᶜ = 0.1, rᵍ = 0.011, Δyᵍ = 27, dᵍˢᵃᵍ = 10),
             earth_parameters = (1,1,25))

    #OHLx2 between AVLGM380 and HORTA380
    ohlAH = overhead_line(length = 54e3,
             conductors = Conductors(organization = :vertical,
                            nᵇ = 6, Rᵈᶜ = 0.05, rᶜ = 0.016,  yᵇᶜ = 30, Δxᵇᶜ = 14, Δyᵇᶜ= 9, dˢᵃᵍ = 10),
             groundwires = Groundwires(nᵍ = 1, Rᵍᵈᶜ = 0.1, rᵍ = 0.011, Δyᵍ = 27, dᵍˢᵃᵍ = 10),
             earth_parameters = (1,1,25))

    #OHLx2 between FR_EQ380 and AVLGM380
    ohlFA = overhead_line(length = 54e3,
             conductors = Conductors(organization = :vertical,
                            nᵇ = 6, Rᵈᶜ = 0.05, rᶜ = 0.016,  yᵇᶜ = 30, Δxᵇᶜ = 14, Δyᵇᶜ= 9, dˢᵃᵍ = 10),
             groundwires = Groundwires(nᵍ = 1, Rᵍᵈᶜ = 0.1, rᵍ = 0.011, Δyᵍ = 27, dᵍˢᵃᵍ = 10),
             earth_parameters = (1,1,25))

    #OHLx2 between AVLGM380 and IZGEM380  Line that changes its composition in respect to the case considered. Case 1: OHLx2
    ohlAI = overhead_line(length = 23e3,
             conductors = Conductors(organization = :vertical,
                            nᵇ = 6, Rᵈᶜ = 0.05, rᶜ = 0.016,  yᵇᶜ = 30, Δxᵇᶜ = 14, Δyᵇᶜ= 9, dˢᵃᵍ = 10),
             groundwires = Groundwires(nᵍ = 1, Rᵍᵈᶜ = 0.1, rᵍ = 0.011, Δyᵍ = 27, dᵍˢᵃᵍ = 10),
             earth_parameters = (1,1,25))

    #OHLx2 between IGZEM380 and GEZEL380  Line that changes its composition in respect to the case considered. Case 1: OHLx2
    ohlIG = overhead_line(length = 52e3,
             conductors = Conductors(organization = :vertical,
                            nᵇ = 6, Rᵈᶜ = 0.05, rᶜ = 0.016,  yᵇᶜ = 30, Δxᵇᶜ = 14, Δyᵇᶜ= 9, dˢᵃᵍ = 10),
             groundwires = Groundwires(nᵍ = 1, Rᵍᵈᶜ = 0.1, rᵍ = 0.011, Δyᵍ = 27, dᵍˢᵃᵍ = 10),
             earth_parameters = (1,1,25))
    #OHLx2 between GEZEL380 and STEVN380
    ohlGS = overhead_line(length = 8e3,
             conductors = Conductors(organization = :vertical,
                            nᵇ = 6, Rᵈᶜ = 0.05, rᶜ = 0.016,  yᵇᶜ = 30, Δxᵇᶜ = 14, Δyᵇᶜ= 9, dˢᵃᵍ = 10),
    groundwires = Groundwires(nᵍ = 1, Rᵍᵈᶜ = 0.1, rᵍ = 0.011, Δyᵍ = 27, dᵍˢᵃᵍ = 10),
    earth_parameters = (1,1,25))

    # Circuit connections

    # Node9 FR_EQ380 + equivalent impedance of FR_EQ380
    ohlFA[1.1] ⟷ ohlFA[1.4] ⟷ zeqFR[2.1] ⟷ Node91
    ohlFA[1.2] ⟷ ohlFA[1.5] ⟷ zeqFR[2.2] ⟷ Node92
    ohlFA[1.3] ⟷ ohlFA[1.6] ⟷ zeqFR[2.3] ⟷ Node93

    # Connection of equivalent generator FR_EQ380 and corresponding impedance
    acFR_EQ380[1.1] ⟷ zeqFR[1.1] ⟷ eqFR1
    acFR_EQ380[1.2] ⟷ zeqFR[1.2] ⟷ eqFR2
    acFR_EQ380[1.3] ⟷ zeqFR[1.3] ⟷ eqFR3

    # Grounding of the equivalent generators of FR_EQ380
    acFR_EQ380[2.1] ⟷ acFR_EQ380[2.2] ⟷ acFR_EQ380[2.3] ⟷ gndFR_EQ380

    # Node7 AVLGM380
    ohlFA[2.1] ⟷ ohlFA[2.4] ⟷ ohlAH[1.1] ⟷ ohlAH[1.4] ⟷ ohlAI[1.1] ⟷ ohlAI[1.4] ⟷ atAV1[1.1] ⟷ atAV2[1.1] ⟷ Node71 #autotransformer28-29 [1.1]
    ohlFA[2.2] ⟷ ohlFA[2.5] ⟷ ohlAH[1.2] ⟷ ohlAH[1.5] ⟷ ohlAI[1.2] ⟷ ohlAI[1.5] ⟷ atAV1[1.2] ⟷ atAV2[1.2] ⟷ Node72
    ohlFA[2.3] ⟷ ohlFA[2.6] ⟷ ohlAH[1.3] ⟷ ohlAH[1.6] ⟷ ohlAI[1.3] ⟷ ohlAI[1.6] ⟷ atAV1[1.3] ⟷ atAV2[1.3] ⟷ Node73

    # Node8 MERCA380 + equivalent impedance of  MERCA380
    ohlMH[1.1] ⟷ ohlMH[1.4] ⟷ zeqMERCA[2.1] ⟷ Node81
    ohlMH[1.2] ⟷ ohlMH[1.5] ⟷ zeqMERCA[2.2] ⟷ Node82
    ohlMH[1.3] ⟷ ohlMH[1.6] ⟷ zeqMERCA[2.3] ⟷ Node83

    # Connection of equivalent generator acMERCA380 and corresponding impedance
    acMERCA380[1.1] ⟷ zeqMERCA[1.1] ⟷ eqMERCA1
    acMERCA380[1.2] ⟷ zeqMERCA[1.2] ⟷ eqMERCA2
    acMERCA380[1.3] ⟷ zeqMERCA[1.3] ⟷ eqMERCA3

    # Grounding of the equivalent generators of MERCA380
    acMERCA380[2.1] ⟷ acMERCA380[2.2] ⟷ acMERCA380[2.3] ⟷ gndMERCA380

    # Node6 HORTA380
    ohlMH[2.1] ⟷ ohlMH[2.4] ⟷ ohlAH[2.1] ⟷ ohlAH[2.4] ⟷ ohlHM[1.1] ⟷ ohlHE[1.1] ⟷ Node61
    ohlMH[2.2] ⟷ ohlMH[2.5] ⟷ ohlAH[2.2] ⟷ ohlAH[2.5] ⟷ ohlHM[1.2] ⟷ ohlHE[1.2] ⟷ Node62
    ohlMH[2.3] ⟷ ohlMH[2.6] ⟷ ohlAH[2.3] ⟷ ohlAH[2.6] ⟷ ohlHM[1.3] ⟷ ohlHE[1.3] ⟷ Node63

    # Node5 EEKLN380
    ohlHE[2.1] ⟷ ohlEM[1.1] ⟷ atEE1[1.1] ⟷ Node51 #autotransformer30 [1.1]
    ohlHE[2.2] ⟷ ohlEM[1.2] ⟷ atEE1[1.2] ⟷ Node52
    ohlHE[2.3] ⟷ ohlEM[1.3] ⟷ atEE1[1.3] ⟷ Node53

    # Node3 MAERL380
    ohlHM[2.1] ⟷ ohlEM[2.1] ⟷ cMG1[1.1] ⟷ cMG2[1.1] ⟷ cMG3[1.1] ⟷ cMG4[1.1] ⟷ Node31
    ohlHM[2.2] ⟷ ohlEM[2.2] ⟷ cMG1[1.2] ⟷ cMG2[1.2] ⟷ cMG3[1.2] ⟷ cMG4[1.2] ⟷ Node32
    ohlHM[2.3] ⟷ ohlEM[2.3] ⟷ cMG1[1.3] ⟷ cMG2[1.3] ⟷ cMG3[1.3] ⟷ cMG4[1.3] ⟷ Node33

    # Node4 IZGEM380
    ohlAI[2.1] ⟷ ohlAI[2.4] ⟷ ohlIG[1.1] ⟷ ohlIG[1.4] ⟷ atIZ1[1.1] ⟷ atIZ2[1.1] ⟷ Node41 #autotransformer26-27 [1.1]
    ohlAI[2.2] ⟷ ohlAI[2.5] ⟷ ohlIG[1.2] ⟷ ohlIG[1.5] ⟷ atIZ1[1.2] ⟷ atIZ2[1.2] ⟷ Node42
    ohlAI[2.3] ⟷ ohlAI[2.6] ⟷ ohlIG[1.3] ⟷ ohlIG[1.6] ⟷ atIZ1[1.3] ⟷ atIZ2[1.3] ⟷ Node43

    # Node1 GEZEL380
    cMG1[2.1] ⟷ cMG2[2.1] ⟷ cMG3[2.1] ⟷ cMG4[2.1] ⟷ ohlIG[2.1] ⟷ ohlIG[2.4] ⟷ ohlGS[1.1] ⟷ ohlGS[1.4] ⟷ reactorGE[1.1] ⟷ atGE1[1.1] ⟷ atGE2[1.1] ⟷ atGE3[1.1] ⟷ atGE4[1.1] ⟷ Node11 #autotransformer16-19 [1.1]
    cMG1[2.2] ⟷ cMG2[2.2] ⟷ cMG3[2.2] ⟷ cMG4[2.2] ⟷ ohlIG[2.2] ⟷ ohlIG[2.5] ⟷ ohlGS[1.2] ⟷ ohlGS[1.5] ⟷ reactorGE[1.2] ⟷ atGE1[1.2] ⟷ atGE2[1.2] ⟷ atGE3[1.2] ⟷ atGE4[1.2] ⟷ Node12
    cMG1[2.3] ⟷ cMG2[2.3] ⟷ cMG3[2.3] ⟷ cMG4[2.3] ⟷ ohlIG[2.3] ⟷ ohlIG[2.6] ⟷ ohlGS[1.3] ⟷ ohlGS[1.6] ⟷ reactorGE[1.3] ⟷ atGE1[1.3] ⟷ atGE2[1.3] ⟷ atGE3[1.3] ⟷ atGE4[1.3] ⟷ Node13
    # Grounding for the equivalent reactor at GEZEL380
    reactorGE[2.1] ⟷ reactorGE[2.2] ⟷ reactorGE[2.3] ⟷ gndGE

    # Node2 STEVN380 + momentaneous grounding Node notation is excluded otherwise gives errors with gnd
    ohlGS[2.1] ⟷ ohlGS[2.4] ⟷ reactorSTEVN[1.1] ⟷ atST1[1.1] ⟷ atST2[1.1] ⟷ atST3[1.1] ⟷ atST4[1.1] ⟷ atST5[1.1] ⟷ atST6[1.1] ⟷ Node21 #autotransformer20-25 [1.1]
    ohlGS[2.2] ⟷ ohlGS[2.5] ⟷ reactorSTEVN[1.2] ⟷ atST1[1.2] ⟷ atST2[1.2] ⟷ atST3[1.2] ⟷ atST4[1.2] ⟷ atST5[1.2] ⟷ atST6[1.2] ⟷ Node22
    ohlGS[2.3] ⟷ ohlGS[2.6] ⟷ reactorSTEVN[1.3] ⟷ atST1[1.3] ⟷ atST2[1.3] ⟷ atST3[1.3] ⟷ atST4[1.3] ⟷ atST5[1.3] ⟷ atST6[1.3] ⟷ Node23
    # Grounding for the equivalent reactor at GEZEL380
    reactorSTEVN[2.1] ⟷ reactorSTEVN[2.2] ⟷ reactorSTEVN[2.3] ⟷ gndSTEVN

    # Node14 AVLGM150
    loadAV1[1.1] ⟷ atAV1[2.1] ⟷ atAV2[2.1] ⟷ Node141 #autotransformer28-29 [2.1]
    loadAV1[1.2] ⟷ atAV1[2.2] ⟷ atAV2[2.2] ⟷ Node142
    loadAV1[1.3] ⟷ atAV1[2.3] ⟷ atAV2[2.3] ⟷ Node143
    #Grounding for the eqyuvalent load connected to AVLGM150
    loadAV1[2.1] ⟷ loadAV1[2.2] ⟷ loadAV1[2.3] ⟷ gndLAV1

    # Node10 GEZEL220
    loadGE2[1.1] ⟷ atGE1[2.1] ⟷ atGE2[2.1] ⟷ atGE3[2.1] ⟷ atGE4[2.1] ⟷ Node101 #autotransformer16-19[2.1]
    loadGE2[1.2] ⟷ atGE1[2.2] ⟷ atGE2[2.2] ⟷ atGE3[2.2] ⟷ atGE4[2.2] ⟷ Node102
    loadGE2[1.3] ⟷ atGE1[2.3] ⟷ atGE2[2.3] ⟷ atGE3[2.3] ⟷ atGE4[2.3] ⟷ Node103
    #Grounding for the equivalent load connected to GEZEL220
    loadGE2[2.1] ⟷ loadGE2[2.2] ⟷ loadGE2[2.3] ⟷ gndLGE2

    # Node13 IZGEM150
    loadIZ1[1.1] ⟷ atIZ1[2.1] ⟷ atIZ2[2.1] ⟷ Node131 #autotransformer26-27 [2.1]
    loadIZ1[1.2] ⟷ atIZ1[2.2] ⟷ atIZ2[2.2] ⟷ Node132
    loadIZ1[1.3] ⟷ atIZ1[2.3] ⟷ atIZ2[2.3] ⟷ Node133
    #Grounding for the equivalent load connected to IZGEM150
    loadIZ1[2.1] ⟷ loadIZ1[2.2] ⟷ loadIZ1[2.3] ⟷ gndLIZ1

    # Node11 STEVN220
    loadST2[1.1] ⟷ atST1[2.1] ⟷ atST2[2.1] ⟷ atST3[2.1] ⟷ atST4[2.1] ⟷ Node111 #autotransformer20-23 [2.1]
    loadST2[1.2] ⟷ atST1[2.2] ⟷ atST2[2.2] ⟷ atST3[2.2] ⟷ atST4[2.2] ⟷ Node112
    loadST2[1.3] ⟷ atST1[2.3] ⟷ atST2[2.3] ⟷ atST3[2.3] ⟷ atST4[2.3] ⟷ Node113
    #Grounding for the equivalent load connected to STEVN220
    loadST2[2.1] ⟷ loadST2[2.2] ⟷ loadST2[2.3] ⟷ gndLST2

    # Node12 STEVN150
    loadST1[1.1] ⟷ atST5[2.1] ⟷ atST6[2.1] ⟷ Node121 #autotransformer24-25 [2.1]
    loadST1[1.2] ⟷ atST5[2.2] ⟷ atST6[2.2] ⟷ Node123
    loadST1[1.3] ⟷ atST5[2.3] ⟷ atST6[2.3] ⟷ Node122
    #Grounding for the equivalent load connected to STEVN150
    loadST1[2.1] ⟷ loadST1[2.2] ⟷ loadST1[2.3] ⟷ gndLST1

    # Node15 EEKLN150
    loadEE1[1.1] ⟷ atEE1[2.1] ⟷ Node151 #autotransformer30 [2.1]
    loadEE1[1.2] ⟷ atEE1[2.2] ⟷ Node152 #autotransformer30
    loadEE1[1.3] ⟷ atEE1[2.3] ⟷ Node153 #autotransformer30
    #Grounding for the equivalent load connected to EEKLN150
    loadEE1[2.1] ⟷ loadEE1[2.2] ⟷ loadEE1[2.3] ⟷ gndLEE1


end
#imp, omega = determine_impedance(net, elim_elements = [:acFR_EQ380], input_pins = Any[:Node31, :Node32, :Node33],
                            #output_pins = Any[:gnd1, :gnd2, :gnd3], omega_range = (-3,5,1000))
#bode(imp, omega = omega, axis_type = :loglog)
