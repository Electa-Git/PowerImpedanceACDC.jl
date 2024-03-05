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
Pmmc2 = -200
Qref = 0
Vm = 222 / sqrt(3) #Vln,rms
Vdc = 640
Ztrafo_base = 222^2/1000
Ltrafo = 0.18*Ztrafo_base/2/pi/50
Rtrafo = 0.005*Ztrafo_base

net = @network begin

    Zg3 = impedance(z = 0.1 + 0.005*s, pins = 3, transformation = true)

    OWF = synchronousmachine(V = Vm, Vᵃᶜ_base = 380, P_min = -2000, P_max = 2000, Q_max = 2000, Q_min = -2000, P = Powf)
    G3 = ac_source(pins = 3, V = Vm, transformation = true)

    MMC2 = mmc(Vᵈᶜ = Vdc, vDCbase = 640, vACbase_LL_RMS = 222, Vₘ = Vm, Lᵣ = Ltrafo, Rᵣ = Rtrafo,
        P_max = 1500, P_min = -1500, P = Pmmc2, Q = Qref, Q_max = 500, Q_min = -500,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        p = PI_control(Kₚ = 0.1, Kᵢ = 31.4159),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )  

    MMC3 = mmc(Vᵈᶜ = Vdc, vDCbase = 640, vACbase_LL_RMS = 222, Vₘ = Vm, Lᵣ = Ltrafo, Rᵣ = Rtrafo,
        P_max = 1500, P_min = -1500, P = Pmmc3, Q = Qref, Q_max = 500, Q_min = -500,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159),
        dc = PI_control(Kₚ = 5, Kᵢ = 15)
        )       

    CableDC23 = cable(length = 60e3, positions = [(-0.5,1), (0.5,1)],
        C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
        C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
        C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
        I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
        I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
        I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)   

    CableAC23 = cable(length = 60e3, positions = [(-0.5,1), (0,1), (0.5,1)],
        C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
        C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
        C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
        I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
        I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
        I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)   

    CableAC2g2 = cable(length = 10e3, positions = [(-0.5,1), (0,1), (0.5,1)],
        C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8),
        C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
        C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
        I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
        I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
        I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3), transformation = true)

    OHLAC3g3 = overhead_line(length = 50e3,
        conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
                        Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
        groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
        earth_parameters = (1,1,100), transformation = true)

    OWF[2.1] ⟷ gndD
    OWF[2.2] ⟷ gndQ
    G3[2.1] ⟷ gndD
    G3[2.2] ⟷ gndQ

    MMC2[1.1] ⟷ CableDC23[1.1] ⟷ NodeDC2
    MMC3[1.1] ⟷ CableDC23[2.1] ⟷ NodeDC3

    MMC2[2.1] ⟷ CableAC23[1.1] ⟷ CableAC2g2[1.1] ⟷ NodeD2
    MMC2[2.2] ⟷ CableAC23[1.2] ⟷ CableAC2g2[1.2] ⟷ NodeQ2

    MMC3[2.1] ⟷ CableAC23[2.1] ⟷ OHLAC3g3[1.1] ⟷ NodeD3
    MMC3[2.2] ⟷ CableAC23[2.2] ⟷ OHLAC3g3[1.2] ⟷ NodeQ3

    OWF[1.1] ⟷ CableAC2g2[2.1] ⟷ NodeDg2
    OWF[1.2] ⟷ CableAC2g2[2.2] ⟷ NodeQg2

    G3[1.1] ⟷ Zg3[1.1]
    G3[1.2] ⟷ Zg3[1.2]
    Zg3[2.1] ⟷ OHLAC3g3[2.1] ⟷ NodeDg3
    Zg3[2.2] ⟷ OHLAC3g3[2.2] ⟷ NodeQg3
    
end

# omega_min = 0
# omega_max = 4
# omegaₙ = 1000

# Zg2DQ, omega = determine_impedance(net, elim_elements = [:CableAC2g2], input_pins = Any[:NodeDg2, :NodeQg2], output_pins = Any[:gndD, :gndQ], omega_range = (omega_min, omega_max, omegaₙ)) 
# Yg2DQ = inv.(Zg2DQ) 
# f = omega/(2*pi)

# Zg3DQ, omega = determine_impedance(net, elim_elements = [:OHLAC3g3], input_pins = Any[:NodeDg3, :NodeQg3], output_pins = Any[:gndD, :gndQ], omega_range = (omega_min, omega_max, omegaₙ)) 
# Yg3DQ = inv.(Zg3DQ)
# f = omega/(2*pi)

# Y_MMC2 = []
# Y_dc2 = []
# for i in 1:length(omega)
#     Y = eval_abcd(net.elements[:MMC2].element_value, 1im*omega[i]) 
#     push!(Y_MMC2, [transpose(Y[1,:]); transpose(-Y[2,:]); transpose(-Y[3,:])]) 
#     push!(Y_dc2, Y[1,1]) 
# end

# bode_Ydc2 = bodeplot(Y_dc2, omega, legend = "Y_dc2")
# display(bode_Ydc2)

# Y_MMC3 = []
# for i in 1:length(omega)
#     Y = eval_abcd(net.elements[:MMC3].element_value, 1im*omega[i]) 
#     push!(Y_MMC3, [transpose(Y[1,:]); transpose(-Y[2,:]); transpose(-Y[3,:])]) 
# end

# omega = 2*pi*f

# θ₂ = net.elements[:MMC2].element_value.θ
# θ₃ = net.elements[:MMC3].element_value.θ

# Tdq_c2 = [1 0 0; 0 cos(θ₂) -sin(θ₂); 0 sin(θ₂) cos(θ₂)]
# for i in 1:length(f)
#     Y_MMC2[i] = inv(Tdq_c2)*Y_MMC2[i]*Tdq_c2
# end

# Tdq_c3 = [1 0 0; 0 cos(θ₃) -sin(θ₃); 0 sin(θ₃) cos(θ₃)]
# for i in 1:length(f)
#     Y_MMC3[i] = inv(Tdq_c3)*Y_MMC3[i]*Tdq_c3
# end

# Yn = [] #Nodes
# for i in 1:length(omega)
#     Ynᵢ = zeros(Complex{Float64},10,10)
#     Ynᵢ[1:3,1:3] = Y_MMC2[i]
#     Ynᵢ[4:6,4:6] = Y_MMC3[i]
#     Ynᵢ[7:8,7:8] = Yg2DQ[i]
#     Ynᵢ[9:10,9:10] = Yg3DQ[i]
#     push!(Yn, Ynᵢ)
# end

# ##### ----- Y_BUS_DC ----- #####

# ABCD_CableDC23 = []
# for i in 1:length(omega)
#     ABCD = eval_abcd(net.elements[:CableDC23].element_value, 1im*omega[i]) 
#     n = Int(size(ABCD, 1)/2)
#     (a, b, c, d) = (ABCD[1:n,1:n], ABCD[1:n,n+1:end], ABCD[n+1:end,1:n], ABCD[n+1:end, n+1:end])
#     ABCD = [(a[1,1]+a[2,2]-a[1,2]-a[2,1])/2 (b[1,1]+b[2,2]-b[1,2]-b[2,1]); (c[1,1]+c[2,2]-c[1,2]-c[2,1])/4 (d[1,1]+d[2,2]-d[1,2]-d[2,1])/2] # Section 2.5 Tutorial
#     push!(ABCD_CableDC23, ABCD) 
# end

# Y_BUS_DC = []

# for i in 1:length(omega)
#     ABCD23 = ABCD_CableDC23[i]
#     (a23, b23, c23, d23) = (ABCD23[1,1], ABCD23[1,2], ABCD23[2,1], ABCD23[2,2])
#     push!(Y_BUS_DC, [a23/b23 -1/b23; -1/b23 a23/b23])
# end

# ##### ----- Y_BUS_AC2 ----- #####

# ABCD_CableAC23_DQ = []

# for i in 1:length(omega)
#     ω₀ = 100*π
#     ABCD₁ = eval_abcd(net.elements[:CableAC23].element_value, 1im*omega[i] + 1im*ω₀)
#     ABCD₂ = eval_abcd(net.elements[:CableAC23].element_value, 1im*omega[i] - 1im*ω₀)

#     n = Int(size(ABCD₁, 1)/2)
#     (a₁, b₁, c₁, d₁) = (ABCD₁[1:n,1:n], ABCD₁[1:n,n+1:end], ABCD₁[n+1:end,1:n], ABCD₁[n+1:end, n+1:end])
#     (a₂, b₂, c₂, d₂) = (ABCD₂[1:n,1:n], ABCD₂[1:n,n+1:end], ABCD₂[n+1:end,1:n], ABCD₂[n+1:end, n+1:end])

#     # Implemented Park transformation: d-axis alignment, q-axis lagging, peak-invariant
#     T = [1 -1im; -1im -1]
#     CK = (2/3)*[1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2]
#     CKinv = [1 0; -1/2 sqrt(3)/2; -1/2 -sqrt(3)/2]

#     a_dq = T * CK * a₂ * CKinv * conj(T) + conj(T) * CK * a₁ * CKinv * T 
#     b_dq = T * CK * b₂ * CKinv * conj(T) + conj(T) * CK * b₁ * CKinv * T
#     c_dq = T * CK * c₂ * CKinv * conj(T) + conj(T) * CK * c₁ * CKinv * T
#     d_dq = T * CK * d₂ * CKinv * conj(T) + conj(T) * CK * d₁ * CKinv * T

#     push!(ABCD_CableAC23_DQ, [a_dq b_dq; c_dq d_dq])

# end

# ABCD_CableAC2g2_DQ = []

# for i in 1:length(omega)
#     ω₀ = 100*π
#     ABCD₁ = eval_abcd(net.elements[:CableAC2g2].element_value, 1im*omega[i] + 1im*ω₀)
#     ABCD₂ = eval_abcd(net.elements[:CableAC2g2].element_value, 1im*omega[i] - 1im*ω₀)

#     n = Int(size(ABCD₁, 1)/2)
#     (a₁, b₁, c₁, d₁) = (ABCD₁[1:n,1:n], ABCD₁[1:n,n+1:end], ABCD₁[n+1:end,1:n], ABCD₁[n+1:end, n+1:end])
#     (a₂, b₂, c₂, d₂) = (ABCD₂[1:n,1:n], ABCD₂[1:n,n+1:end], ABCD₂[n+1:end,1:n], ABCD₂[n+1:end, n+1:end])

#     # Implemented Park transformation: d-axis alignment, q-axis lagging, peak-invariant
#     T = [1 -1im; -1im -1]
#     CK = (2/3)*[1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2]
#     CKinv = [1 0; -1/2 sqrt(3)/2; -1/2 -sqrt(3)/2]

#     a_dq = T * CK * a₂ * CKinv * conj(T) + conj(T) * CK * a₁ * CKinv * T 
#     b_dq = T * CK * b₂ * CKinv * conj(T) + conj(T) * CK * b₁ * CKinv * T
#     c_dq = T * CK * c₂ * CKinv * conj(T) + conj(T) * CK * c₁ * CKinv * T
#     d_dq = T * CK * d₂ * CKinv * conj(T) + conj(T) * CK * d₁ * CKinv * T

#     push!(ABCD_CableAC2g2_DQ, [a_dq b_dq; c_dq d_dq])

# end

# ABCD_OHLAC3g3_DQ = []

# for i in 1:length(omega)
#     ω₀ = 100*π
#     ABCD₁ = eval_abcd(net.elements[:OHLAC3g3].element_value, 1im*omega[i] + 1im*ω₀)
#     ABCD₂ = eval_abcd(net.elements[:OHLAC3g3].element_value, 1im*omega[i] - 1im*ω₀)

#     n = Int(size(ABCD₁, 1)/2)
#     (a₁, b₁, c₁, d₁) = (ABCD₁[1:n,1:n], ABCD₁[1:n,n+1:end], ABCD₁[n+1:end,1:n], ABCD₁[n+1:end, n+1:end])
#     (a₂, b₂, c₂, d₂) = (ABCD₂[1:n,1:n], ABCD₂[1:n,n+1:end], ABCD₂[n+1:end,1:n], ABCD₂[n+1:end, n+1:end])

#     # Implemented Park transformation: d-axis alignment, q-axis lagging, peak-invariant
#     T = [1 -1im; -1im -1]
#     CK = (2/3)*[1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2]
#     CKinv = [1 0; -1/2 sqrt(3)/2; -1/2 -sqrt(3)/2]

#     a_dq = T * CK * a₂ * CKinv * conj(T) + conj(T) * CK * a₁ * CKinv * T 
#     b_dq = T * CK * b₂ * CKinv * conj(T) + conj(T) * CK * b₁ * CKinv * T
#     c_dq = T * CK * c₂ * CKinv * conj(T) + conj(T) * CK * c₁ * CKinv * T
#     d_dq = T * CK * d₂ * CKinv * conj(T) + conj(T) * CK * d₁ * CKinv * T

#     push!(ABCD_OHLAC3g3_DQ, [a_dq b_dq; c_dq d_dq])

# end

# Y_BUS_AC2 = []

# for i in 1:length(omega)
#     ABCD23 = ABCD_CableAC23_DQ[i]
#     n = Int(size(ABCD23, 1)/2)
#     (a23, b23, c23, d23) = (ABCD23[1:n,1:n], ABCD23[1:n,n+1:end], ABCD23[n+1:end,1:n], ABCD23[n+1:end, n+1:end])

#     ABCD2g2 = ABCD_CableAC2g2_DQ[i]
#     n = Int(size(ABCD2g2, 1)/2)
#     (a2g2, b2g2, c2g2, d2g2) = (ABCD2g2[1:n,1:n], ABCD2g2[1:n,n+1:end], ABCD2g2[n+1:end,1:n], ABCD2g2[n+1:end, n+1:end])

#     ABCD3g3 = ABCD_OHLAC3g3_DQ[i]
#     n = Int(size(ABCD3g3, 1)/2)
#     (a3g3, b3g3, c3g3, d3g3) = (ABCD3g3[1:n,1:n], ABCD3g3[1:n,n+1:end], ABCD3g3[n+1:end,1:n], ABCD3g3[n+1:end, n+1:end])

#     M0 = zeros(Complex{Float64},2,2)

#     push!(Y_BUS_AC2, [inv(b23)*a23+inv(b2g2)*a2g2 -inv(b23) -inv(b2g2) M0; -inv(b23) inv(b23)*a23+inv(b3g3)*a3g3 M0 -inv(b3g3); -inv(b2g2) M0 inv(b2g2)*a2g2 M0; M0 -inv(b3g3) M0 inv(b3g3)*a3g3])
# end

# # AC/DC bus admittance matrix

# # I = [idc2 id2 iq2 idc3 id3 iq3 idg2 iqg2 idg3 iqg3]^T 
# # V = [vdc2 vd2 vq2 vdc3 vd3 vq3 vdg2 vqg2 vdg3 vqg3]^T
# # I = Y_BUS*V

# Ye = [] #Edges

# for i in 1:length(omega)

#     Yeᵢ = zeros(Complex{Float64},10,10)
    
#     #Y_BUS_DC
#     Yeᵢ[1,1] = Y_BUS_DC[i][1,1]
#     Yeᵢ[1,4] = Y_BUS_DC[i][1,2]
#     Yeᵢ[4,1] = Y_BUS_DC[i][2,1]
#     Yeᵢ[4,4] = Y_BUS_DC[i][2,2]

#     #Y_BUS_AC2
#     Yeᵢ[2:3,2:3] = Y_BUS_AC2[i][1:2,1:2]
#     Yeᵢ[2:3,5:6] = Y_BUS_AC2[i][1:2,3:4]
#     Yeᵢ[2:3,7:8] = Y_BUS_AC2[i][1:2,5:6]
#     Yeᵢ[2:3,9:10] = Y_BUS_AC2[i][1:2,7:8]

#     Yeᵢ[5:6,2:3] = Y_BUS_AC2[i][3:4,1:2]
#     Yeᵢ[5:6,5:6] = Y_BUS_AC2[i][3:4,3:4]
#     Yeᵢ[5:6,7:8] = Y_BUS_AC2[i][3:4,5:6]
#     Yeᵢ[5:6,9:10] = Y_BUS_AC2[i][3:4,7:8]

#     Yeᵢ[7:8,2:3] = Y_BUS_AC2[i][5:6,1:2]
#     Yeᵢ[7:8,5:6] = Y_BUS_AC2[i][5:6,3:4]
#     Yeᵢ[7:8,7:8] = Y_BUS_AC2[i][5:6,5:6]
#     Yeᵢ[7:8,9:10] = Y_BUS_AC2[i][5:6,7:8]

#     Yeᵢ[9:10,2:3] = Y_BUS_AC2[i][7:8,1:2]
#     Yeᵢ[9:10,5:6] = Y_BUS_AC2[i][7:8,3:4]
#     Yeᵢ[9:10,7:8] = Y_BUS_AC2[i][7:8,5:6]
#     Yeᵢ[9:10,9:10] = Y_BUS_AC2[i][7:8,7:8]

#     push!(Ye, Yeᵢ)

# end

# L = inv.(Ye) .* Yn
# nyquist_Energy_Hub = nyquistplot(L, omega, zoom = "no", SM = "VM")
# savefig("nyquist_Energy_Hub.png")

# Zcl_bus = inv.(Ye + Yn)

# fmin = 1
# fmax = 10000

# EVD_Energy_Hu1b = EVD(Zcl_bus, omega, fmin, fmax)
# savefig("EVD_Energy_Hub.png")

# Λ₂ = eigvals(net.elements[:MMC2].element_value.A)
# for i in 1:length(Λ₂)
#     if real(Λ₂[i]) >= 0
#         println("MMC2 RHP pole : ", Λ₂[i])
#     end
# end

# Λ₃ = eigvals(net.elements[:MMC3].element_value.A)
# for i in 1:length(Λ₃)
#     if real(Λ₃[i]) >= 0
#         println("MMC3 RHP pole : ", Λ₃[i])
#     end
# end

# Λₒ = eigvals(net.elements[:OWF].element_value.A)
# for i in 1:length(Λₒ)
#     if real(Λₒ[i]) >= 0
#         println("OWF RHP pole : ", Λₒ[i])
#     end
# end
