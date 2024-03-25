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

    Zg3 = impedance(z = 0.1 + 0.005*s, pins = 3, transformation = true)

    Q2 = impedance(z = Lcomp*s, pins = 3, transformation = true)
    Q3 = impedance(z = Lcomp*s, pins = 3, transformation = true)

    G3 = ac_source(pins = 3, V = Vm, transformation = true)

    OWF = tlc(Vᵈᶜ = Vdc, Vₘ = Vm, Lᵣ = 0.012, Rᵣ = 0.03872, 
        Sbase = 1000, vACbase_LL_RMS = 220, 
        P = Powf, Q = Qowf,
        occ = PI_control(Kₚ = 0.254647908947033, Kᵢ = 0.8),
        # pll = PI_control(Kₚ = 0.8754, Kᵢ = 38.5155, ω_f = 2*pi*80), #110 rad/s
        # pll = PI_control(Kₚ = 0.795774715459477, Kᵢ = 31.830988618379067, ω_f = 2*pi*80), # These PLL gains result in an instability in PSCAD. 
        # This instability can be detected (negative phase margin and a very low vector margin at around 59 Hz, which is close to the oscillation frequency in PSCAD) here if single terminal analysis is used.
        # As of 14/02/24 multi-terminal analysis cannot detect this instability.
        pll = PI_control(Kₚ = 0.397887357729738, Kᵢ = 7.957747154594767, ω_f = 2*pi*80), # These gains are fine.
        v_meas_filt = PI_control(ω_f = 1e4),
        i_meas_filt = PI_control(ω_f = 1e4),
        #vac_supp = PI_control(ω_f = 1/0.5, Kₚ =10),
        #f_supp = PI_control(ω_f = 1/0.5, Kₚ =10),
        p = PI_control(Kₚ = 0.01, Kᵢ = 10),
        q = PI_control(Kₚ = 0.01, Kᵢ = 10)
        )

    MMC2 = mmc(Vᵈᶜ = Vdc, vDCbase = 640, Sbase = 1000, vACbase_LL_RMS = 333, turnsRatio = 333/220, Vₘ = Vm, Lᵣ = Ltrafo, Rᵣ = Rtrafo,
        P_max = 1500, P_min = -1500, P = Pmmc2, Q = Qref, Q_max = 500, Q_min = -500,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        p = PI_control(Kₚ = 0.1, Kᵢ = 31.4159),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )  

    MMC3 = mmc(Vᵈᶜ = Vdc, vDCbase = 640, Sbase = 1000, vACbase_LL_RMS = 333, turnsRatio = 333/220, Vₘ = Vm, Lᵣ = Ltrafo, Rᵣ = Rtrafo,
        P_max = 1500, P_min = -1500, P = -Pmmc2, Q = Qref, Q_max = 500, Q_min = -500,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159),
        dc = PI_control(Kₚ = 5, Kᵢ = 15)
        )        


    CableDC23 = cable(length = 60e3, positions = [(-0.5,1), (0.5,1)],
        C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
        I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
        C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
        I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
        C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),     
        I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)   
        
    CableAC23 = cable(length = 60e3, positions = [(0,1.33562), (-0.0575,1.4322), (0.0575,1.4322)],
        C1 = Conductor(rₒ = 25.45e-3, ρ = 2.63e-8, μᵣ = 1),
        I1 = Insulator(rᵢ = 25.45e-3, rₒ = 50.65e-3, a = 2e-3, b = 1.5e-3, ϵᵣ = 2.3),
        C2 = Conductor(rᵢ = 50.65e-3, rₒ = 52.35e-3, ρ = 22e-8, μᵣ = 1),
        I2 = Insulator(rᵢ = 52.35e-3, rₒ = 55.75e-3, ϵᵣ = 100),
        earth_parameters = (1,1,1), transformation = true)

    CableAC2g2 = cable(length = 10e3, positions = [(0,1.33562), (-0.0575,1.4322), (0.0575,1.4322)],
        C1 = Conductor(rₒ = 25.45e-3, ρ = 2.63e-8, μᵣ = 1),
        I1 = Insulator(rᵢ = 25.45e-3, rₒ = 50.65e-3, a = 2e-3, b = 1.5e-3, ϵᵣ = 2.3),
        C2 = Conductor(rᵢ = 50.65e-3, rₒ = 52.35e-3, ρ = 22e-8, μᵣ = 1),
        I2 = Insulator(rᵢ = 52.35e-3, rₒ = 55.75e-3, ϵᵣ = 100),
        earth_parameters = (1,1,1), transformation = true)

    OHLAC3g3 = overhead_line(length = 50e3,
        conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                        Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        earth_parameters = (1,1,100), transformation = true)

    ######
    # Extra AC/DC converter to represent DC voltage and get power flow converged 
    c2 = mmc(Vᵈᶜ = Vdc, vDCbase = 640, Vₘ = Vm, 
        P_max = 1000, P_min = -1000, P = -Powf, Q = 0, Q_max = 1000, Q_min = -1000,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        dc = PI_control(Kₚ = 5, Kᵢ = 15),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )

    g2 = ac_source(V = Vm, P = -Powf, P_min = -2000, P_max = 2000, Q_max = 1000, Q_min = -1000, pins = 3, transformation = true)

    # TL at the remote end
    tl1 = overhead_line(length = 50e3,
        conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                        Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        earth_parameters = (1,1,100), transformation = true)
    #####

    c2[2.1] ⟷ tl1[2.1]
    c2[2.2] ⟷ tl1[2.2]

    g2[1.1] ⟷ tl1[1.1]
    g2[1.2] ⟷ tl1[1.2]

    g2[2.1] ⟷ gndD
    g2[2.2] ⟷ gndQ
    
    OWF[1.1] ⟷ c2[1.1]

    G3[2.1] ⟷ gndD
    G3[2.2] ⟷ gndQ

    Q2[2.1] ⟷ gndD
    Q2[2.2] ⟷ gndQ
    Q3[2.1] ⟷ gndD
    Q3[2.2] ⟷ gndQ

    MMC2[1.1] ⟷ CableDC23[1.1] ⟷ NodeDC2
    MMC3[1.1] ⟷ CableDC23[2.1] ⟷ NodeDC3

    MMC2[2.1] ⟷ CableAC23[1.1] ⟷ CableAC2g2[1.1] ⟷ Q2[1.1] ⟷ NodeD2
    MMC2[2.2] ⟷ CableAC23[1.2] ⟷ CableAC2g2[1.2] ⟷ Q2[1.2] ⟷ NodeQ2

    MMC3[2.1] ⟷ CableAC23[2.1] ⟷ OHLAC3g3[1.1] ⟷ Q3[1.1] ⟷ NodeD3
    MMC3[2.2] ⟷ CableAC23[2.2] ⟷ OHLAC3g3[1.2] ⟷ Q3[1.2] ⟷ NodeQ3

    OWF[2.1] ⟷ CableAC2g2[2.1] ⟷ NodeDg2
    OWF[2.2] ⟷ CableAC2g2[2.2] ⟷ NodeQg2

    # With external grid impedance
    G3[1.1] ⟷ Zg3[1.1] ⟷ NodeD4
    G3[1.2] ⟷ Zg3[1.2] ⟷ NodeQ4
    Zg3[2.1] ⟷ OHLAC3g3[2.1] ⟷ NodeDg3
    Zg3[2.2] ⟷ OHLAC3g3[2.2] ⟷ NodeQg3

    # Without external grid impedance
    # G3[1.1] ⟷ OHLAC3g3[2.1] ⟷ NodeDg3
    # G3[1.2] ⟷ OHLAC3g3[2.2] ⟷ NodeQg3
    
end

if singleTerminalAnalysis
    ##### ----- STABILITY ANALYSIS WITH DETERMINE_IMPEDANCE ----- #####

    # DC-side
    # Z_MMC_DC1, omega = determine_impedance(net, elim_elements = [:CableDC12], input_pins = Any[:NodeDC1], output_pins = Any[:gndD], omega_range = (-3, 4, 1000))
    # Z_BUS_DC1, omega = determine_impedance(net, elim_elements = [:MMC1], input_pins = Any[:NodeDC1], output_pins = Any[:gndD], omega_range = (-3, 4, 1000))
    # L_DC1 = Z_BUS_DC1 ./ Z_MMC_DC1

    # nyquist_L_DC1 = nyquistplot(L_DC1, omega, zoom = "no", SM = "VM") 

    # AC-side
    Z_MMC_AC1, omega = determine_impedance(net, elim_elements = [:CableAC2g2], input_pins = Any[:NodeDg2, :NodeQg2], output_pins = Any[:gndD, :gndQ], omega_range = (-3, 4, 1000))
    Z_BUS_AC1, omega = determine_impedance(net, elim_elements = [:OWF], input_pins = Any[:NodeDg2, :NodeQg2], output_pins = Any[:gndD, :gndQ], omega_range = (-3, 4, 1000))
    L_AC1 = Z_BUS_AC1 ./ Z_MMC_AC1

    nyquist_L_AC1 = nyquistplot(L_AC1, omega, zoom = "yes", SM = "PM") 
end

##### ----- ADMITTANCE MATRIX NODES ----- #####

# C. Zhang, M. Molinas, A. Rygg, and X. Cai, “Impedance-based analysis of interconnected power electronics systems: 
# Impedance network modeling and comparative studies of stability criteria,” 
# IEEE J. Emerg. Sel. Top. Power Electron., vol. 8, no. 3, pp. 2520–2533, Sep. 2020, doi: 10.1109/JESTPE.2019.2914560.

# General dq-reference frame 
# The dq frame is with d-axis alignment to alpha and q-axis leading the d-axis?
# Tdq(θ₁) = [cos(θ₁) sin(θ₁); -sin(θ₁) cos(θ₁)]

# Julia: The dq frame is with d-axis alignment to alpha and q-axis lagging the d-axis
# Tdq(θ₁) = [cos(θ₁) -sin(θ₁); sin(θ₁) cos(θ₁)]
# Tdq_inv(θ₁) = [cos(θ₁) sin(θ₁); -sin(θ₁) cos(θ₁)]
# Page 30 of tutorial, MMC.jl line 205-206

#Zdq_global = Tdq(-θ₁)*Zdq_local*Tdq(θ₁)

# 3x3 matrix converter 
# Tdq_c(θ₁) = [1 0 0; 0 cos(θ₁) -sin(θ₁); 0 sin(θ₁) cos(θ₁)]
# Zdq_c_global(θ₁) = Tdq_c(-θ₁)*Zdq_c_local*Tdq_c(θ₁)

# The HVDC link decouples the ac systems in terms of reference frames. 
# For each ac system, the ac impedances in that area should be transformed into 
# a reference frame provided by the corresponding ac/dc interface, e.g., the sending- or receiving-VSC

#=
# Use of two-level converter input data
FS_2LC = readdlm("Y_2LC.txt", ComplexF64) 
omegaₙ = size(FS_2LC,1)
n = size(FS_2LC,2)
dim = Int(sqrt(n-1))
omega = real(FS_2LC[:,1])

omega_min = log10(omega[1]/(2*pi))
omega_max = log10(omega[omegaₙ]/(2*pi))

Y_2LC = []
for i in 1:omegaₙ
    push!(Y_2LC, Tdq_inv*transpose(reshape(FS_2LC[i,2:n],dim,dim))*Tdq)
end
=#

@time if multiTerminalAnalysis
    omega_min = 0
    omega_max = 4
    omegaₙ = 1000

    
    # Aggregated OWF (represented using a single TLC)
    Zg2DQ, omega = determine_impedance(net, elim_elements = [:CableAC2g2], input_pins = Any[:NodeDg2, :NodeQg2], output_pins = Any[:gndD, :gndQ], omega_range = (omega_min, omega_max, omegaₙ)) 
    Yg2DQ = inv.(Zg2DQ)
    
    f = omega/(2*pi)

    # Yg2DQ_Data = readdlm("Y_TLC.txt", ',', ComplexF64)
    # #f = real(Yg2DQ_Data[:,1])
    # dim = size(Yg2DQ_Data[1,:])[1]-1
    # Yg2DQ_Data = Yg2DQ_Data[:,2:dim+1]
    # Yg2DQ = []
    # dim = Int(sqrt(dim))
    # for i in 1:length(f)
    #     push!(Yg2DQ, transpose(reshape(Yg2DQ_Data[i,:],dim,dim)))
    # end

    θₒ = -net.elements[:OWF].element_value.θ

    Tdq_co = [cos(θₒ) -sin(θₒ); sin(θₒ) cos(θₒ)]
    for i in 1:length(f)
        Yg2DQ[i] = inv(Tdq_co)*Yg2DQ[i]*Tdq_co
    end
    # External grid 2, representing BE grid (bottom left)
    Zg3DQ, omega = determine_impedance(net, elim_elements = [:OHLAC3g3], input_pins = Any[:NodeDg3, :NodeQg3], output_pins = Any[:gndD, :gndQ], omega_range = (omega_min, omega_max, omegaₙ)) 
    Yg3DQ = inv.(Zg3DQ)

    Y_MMC2 = []
    Y_MMC3 = []
    # Not sure why these two are needed, disabled for now.
    # Y_dc1 = []
    # Y_dc2 = []

    

    for i in 1:length(omega)
        Y2 = eval_abcd(net.elements[:MMC2].element_value, 1im*omega[i]) 
        push!(Y_MMC2, [transpose(Y2[1,:]); transpose(-Y2[2,:]); transpose(-Y2[3,:])]) 
        # push!(Y_dc2, Y2[1,1]) 
        Y3 = eval_abcd(net.elements[:MMC3].element_value, 1im*omega[i]) 
        push!(Y_MMC3, [transpose(Y3[1,:]); transpose(-Y3[2,:]); transpose(-Y3[3,:])]) 
    end

    θ₂ = -net.elements[:MMC2].element_value.θ
    θ₃ = -net.elements[:MMC3].element_value.θ

    Tdq_c2 = [1 0 0; 0 cos(θ₂) -sin(θ₂); 0 sin(θ₂) cos(θ₂)]
    Tdq_c3 = [1 0 0; 0 cos(θ₃) -sin(θ₃); 0 sin(θ₃) cos(θ₃)]
    for i in 1:length(f)
        Y_MMC2[i] = inv(Tdq_c2)*Y_MMC2[i]*Tdq_c2
        Y_MMC3[i] = inv(Tdq_c3)*Y_MMC3[i]*Tdq_c3
    end

    Yn = [] #Nodes
    for i in 1:length(omega)
        Ynᵢ = zeros(Complex{Float64},10,10)
        Ynᵢ[1:3,1:3] = Y_MMC2[i]
        Ynᵢ[4:6,4:6] = Y_MMC3[i]
        Ynᵢ[7:8,7:8] = Yg2DQ[i]
        Ynᵢ[9:10,9:10] = Yg3DQ[i]
        push!(Yn, Ynᵢ)
    end

    ##### ----- Y_BUS_DC ----- #####

    ABCD_CableDC23 = []
    for i in 1:length(omega)
        ABCD = eval_abcd(net.elements[:CableDC23].element_value, 1im*omega[i]) 
        n = Int(size(ABCD, 1)/2)
        (a, b, c, d) = (ABCD[1:n,1:n], ABCD[1:n,n+1:end], ABCD[n+1:end,1:n], ABCD[n+1:end, n+1:end])
        ABCD = [(a[1,1]+a[2,2]-a[1,2]-a[2,1])/2 (b[1,1]+b[2,2]-b[1,2]-b[2,1]); (c[1,1]+c[2,2]-c[1,2]-c[2,1])/4 (d[1,1]+d[2,2]-d[1,2]-d[2,1])/2] # Section 2.5 Tutorial
        push!(ABCD_CableDC23, ABCD) 
    end

    Y_BUS_DC_old = []

    for i in 1:length(omega)
        ABCD23 = ABCD_CableDC23[i]
        push!(Y_BUS_DC_old, HVDCstability.abcd_to_y(ABCD23))
    end

    Y_BUS_DC = HVDCstability.make_y_matrix(net,elim_elements = [:MMC2,:MMC3], input_pins = Any[:NodeDC2,:NodeDC3], omega_range = (omega_min, omega_max, omegaₙ))

    ##### ----- Y_BUS_AC2 ----- #####

    Y_BUS_AC2 =  HVDCstability.make_y_matrix(net,elim_elements = [:Zg2,:Zg3],input_pins = Any[:NodeD2,:NodeQ2,:NodeD3,:NodeQ3,:NodeDg2,:NodeQg2,:NodeDg3,:NodeQg3], omega_range = (omega_min, omega_max, omegaₙ))

    # AC/DC bus admittance matrix

    # I = [idc2 id2 iq2 idc3 id3 iq3 idg2 iqg2 idg3 iqg3]^T 
    # V = [vdc2 vd2 vq2 vdc3 vd3 vq3 vdg2 vqg2 vdg3 vqg3]^T
    # I = Y_BUS*V

    Ye = [] #Edges

    for i in 1:length(omega)

        Yeᵢ = zeros(Complex{Float64},10,10)
        
        #Y_BUS_DC
        Yeᵢ[1,1] = Y_BUS_DC[i][1,1]
        Yeᵢ[1,4] = Y_BUS_DC[i][1,2]
        Yeᵢ[4,1] = Y_BUS_DC[i][2,1]
        Yeᵢ[4,4] = Y_BUS_DC[i][2,2]

        #Y_BUS_AC2
        Yeᵢ[2:3,2:3] = Y_BUS_AC2[i][1:2,1:2]
        Yeᵢ[2:3,5:6] = Y_BUS_AC2[i][1:2,3:4]
        Yeᵢ[2:3,7:8] = Y_BUS_AC2[i][1:2,5:6]
        Yeᵢ[2:3,9:10] = Y_BUS_AC2[i][1:2,7:8]

        Yeᵢ[5:6,2:3] = Y_BUS_AC2[i][3:4,1:2]
        Yeᵢ[5:6,5:6] = Y_BUS_AC2[i][3:4,3:4]
        Yeᵢ[5:6,7:8] = Y_BUS_AC2[i][3:4,5:6]
        Yeᵢ[5:6,9:10] = Y_BUS_AC2[i][3:4,7:8]

        Yeᵢ[7:8,2:3] = Y_BUS_AC2[i][5:6,1:2]
        Yeᵢ[7:8,5:6] = Y_BUS_AC2[i][5:6,3:4]
        Yeᵢ[7:8,7:8] = Y_BUS_AC2[i][5:6,5:6]
        Yeᵢ[7:8,9:10] = Y_BUS_AC2[i][5:6,7:8]

        Yeᵢ[9:10,2:3] = Y_BUS_AC2[i][7:8,1:2]
        Yeᵢ[9:10,5:6] = Y_BUS_AC2[i][7:8,3:4]
        Yeᵢ[9:10,7:8] = Y_BUS_AC2[i][7:8,5:6]
        Yeᵢ[9:10,9:10] = Y_BUS_AC2[i][7:8,7:8]

        push!(Ye, Yeᵢ)

    end

    L = inv.(Ye) .* Yn
    nyquist_Energy_Hub = nyquistplot(L, omega, zoom = "yes", SM = "VM")
    # savefig("nyquist_Energy_Hub.png")

    Zcl_bus = inv.(Ye + Yn)

    fmin = 1
    fmax = 10000

    EVD_Energy_Hu1b = EVD(Zcl_bus, omega, fmin, fmax)
    # savefig("EVD_Energy_Hub.png")

    #=
    # Export to L_Matlab
    L_Matlab = []
    for i in 1:length(omega)
        Lᵢ = L[i]
        push!(L_Matlab, [transpose(Lᵢ[1,:]) transpose(Lᵢ[2,:]) transpose(Lᵢ[3,:]) transpose(Lᵢ[4,:]) transpose(Lᵢ[5,:]) transpose(Lᵢ[6,:]) transpose(Lᵢ[7,:]) transpose(Lᵢ[8,:]) transpose(Lᵢ[9,:]) transpose(Lᵢ[10,:]) transpose(Lᵢ[11,:]) transpose(Lᵢ[12,:]) transpose(Lᵢ[13,:]) transpose(Lᵢ[14,:]) transpose(Lᵢ[15,:]) omega[i]])
    end
    writedlm("L_Matlab.csv", L_Matlab, ',')
    =#

    Λ₂ = eigvals(net.elements[:MMC2].element_value.A)
    for i in 1:length(Λ₂)
        if real(Λ₂[i]) >= 0
            println("MMC2 RHP pole : ", Λ₂[i])
        end
    end

    Λ₃ = eigvals(net.elements[:MMC3].element_value.A)
    for i in 1:length(Λ₃)
        if real(Λ₃[i]) >= 0
            println("MMC3 RHP pole : ", Λ₃[i])
        end
    end

    Λₒ = eigvals(net.elements[:OWF].element_value.A)
    for i in 1:length(Λₒ)
        if real(Λₒ[i]) >= 0
            println("OWF RHP pole : ", Λₒ[i])
        end
    end
end
