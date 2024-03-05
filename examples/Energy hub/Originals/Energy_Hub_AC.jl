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
Vm = 380 / sqrt(3) #Vln,rms

net = @network begin

    Zg3 = impedance(z = 0.1 + 0.005*s, pins = 3, transformation = true)

    OWF = synchronousmachine(V = Vm, Vᵃᶜ_base = 380, P_min = -2000, P_max = 2000, Q_max = 2000, Q_min = -2000, P = Powf)
    G3 = ac_source(pins = 3, V = Vm, transformation = true)

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

    CableAC23[1.1] ⟷ CableAC2g2[1.1] ⟷ NodeD2
    CableAC23[1.2] ⟷ CableAC2g2[1.2] ⟷ NodeQ2

    CableAC23[2.1] ⟷ OHLAC3g3[1.1] ⟷ NodeD3
    CableAC23[2.2] ⟷ OHLAC3g3[1.2] ⟷ NodeQ3

    OWF[1.1] ⟷ CableAC2g2[2.1] ⟷ NodeDg2
    OWF[1.2] ⟷ CableAC2g2[2.2] ⟷ NodeQg2

    G3[1.1] ⟷ Zg3[1.1]
    G3[1.2] ⟷ Zg3[1.2]
    Zg3[2.1] ⟷ OHLAC3g3[2.1] ⟷ NodeDg3
    Zg3[2.2] ⟷ OHLAC3g3[2.2] ⟷ NodeQg3
    
end

omega_min = 0
omega_max = 4
omegaₙ = 1000

Zg2DQ, omega = determine_impedance(net, elim_elements = [:CableAC2g2], input_pins = Any[:NodeDg2, :NodeQg2], output_pins = Any[:gndD, :gndQ], omega_range = (omega_min, omega_max, omegaₙ)) 
Yg2DQ = inv.(Zg2DQ) 
f = omega/(2*pi)

Zg3DQ, omega = determine_impedance(net, elim_elements = [:OHLAC3g3], input_pins = Any[:NodeDg3, :NodeQg3], output_pins = Any[:gndD, :gndQ], omega_range = (omega_min, omega_max, omegaₙ)) 
Yg3DQ = inv.(Zg3DQ)
f = omega/(2*pi)

Yn = [] #Nodes
for i in 1:length(omega)
    Ynᵢ = zeros(Complex{Float64},4,4)
    Ynᵢ[1:2,1:2] = Yg2DQ[i]
    Ynᵢ[3:4,3:4] = Yg3DQ[i]
    push!(Yn, Ynᵢ)
end

##### ----- Y_BUS_AC2 ----- #####

ABCD_CableAC23_DQ = []

for i in 1:length(omega)
    ω₀ = 100*π
    ABCD₁ = eval_abcd(net.elements[:CableAC23].element_value, 1im*omega[i] + 1im*ω₀)
    ABCD₂ = eval_abcd(net.elements[:CableAC23].element_value, 1im*omega[i] - 1im*ω₀)

    n = Int(size(ABCD₁, 1)/2)
    (a₁, b₁, c₁, d₁) = (ABCD₁[1:n,1:n], ABCD₁[1:n,n+1:end], ABCD₁[n+1:end,1:n], ABCD₁[n+1:end, n+1:end])
    (a₂, b₂, c₂, d₂) = (ABCD₂[1:n,1:n], ABCD₂[1:n,n+1:end], ABCD₂[n+1:end,1:n], ABCD₂[n+1:end, n+1:end])

    # Implemented Park transformation: d-axis alignment, q-axis lagging, peak-invariant
    T = [1 -1im; -1im -1]
    CK = (2/3)*[1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2]
    CKinv = [1 0; -1/2 sqrt(3)/2; -1/2 -sqrt(3)/2]

    a_dq = T * CK * a₂ * CKinv * conj(T) + conj(T) * CK * a₁ * CKinv * T 
    b_dq = T * CK * b₂ * CKinv * conj(T) + conj(T) * CK * b₁ * CKinv * T
    c_dq = T * CK * c₂ * CKinv * conj(T) + conj(T) * CK * c₁ * CKinv * T
    d_dq = T * CK * d₂ * CKinv * conj(T) + conj(T) * CK * d₁ * CKinv * T

    push!(ABCD_CableAC23_DQ, [a_dq b_dq; c_dq d_dq])

end

ABCD_CableAC2g2_DQ = []

for i in 1:length(omega)
    ω₀ = 100*π
    ABCD₁ = eval_abcd(net.elements[:CableAC2g2].element_value, 1im*omega[i] + 1im*ω₀)
    ABCD₂ = eval_abcd(net.elements[:CableAC2g2].element_value, 1im*omega[i] - 1im*ω₀)

    n = Int(size(ABCD₁, 1)/2)
    (a₁, b₁, c₁, d₁) = (ABCD₁[1:n,1:n], ABCD₁[1:n,n+1:end], ABCD₁[n+1:end,1:n], ABCD₁[n+1:end, n+1:end])
    (a₂, b₂, c₂, d₂) = (ABCD₂[1:n,1:n], ABCD₂[1:n,n+1:end], ABCD₂[n+1:end,1:n], ABCD₂[n+1:end, n+1:end])

    # Implemented Park transformation: d-axis alignment, q-axis lagging, peak-invariant
    T = [1 -1im; -1im -1]
    CK = (2/3)*[1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2]
    CKinv = [1 0; -1/2 sqrt(3)/2; -1/2 -sqrt(3)/2]

    a_dq = T * CK * a₂ * CKinv * conj(T) + conj(T) * CK * a₁ * CKinv * T 
    b_dq = T * CK * b₂ * CKinv * conj(T) + conj(T) * CK * b₁ * CKinv * T
    c_dq = T * CK * c₂ * CKinv * conj(T) + conj(T) * CK * c₁ * CKinv * T
    d_dq = T * CK * d₂ * CKinv * conj(T) + conj(T) * CK * d₁ * CKinv * T

    push!(ABCD_CableAC2g2_DQ, [a_dq b_dq; c_dq d_dq])

end

ABCD_OHLAC3g3_DQ = []

for i in 1:length(omega)
    ω₀ = 100*π
    ABCD₁ = eval_abcd(net.elements[:OHLAC3g3].element_value, 1im*omega[i] + 1im*ω₀)
    ABCD₂ = eval_abcd(net.elements[:OHLAC3g3].element_value, 1im*omega[i] - 1im*ω₀)

    n = Int(size(ABCD₁, 1)/2)
    (a₁, b₁, c₁, d₁) = (ABCD₁[1:n,1:n], ABCD₁[1:n,n+1:end], ABCD₁[n+1:end,1:n], ABCD₁[n+1:end, n+1:end])
    (a₂, b₂, c₂, d₂) = (ABCD₂[1:n,1:n], ABCD₂[1:n,n+1:end], ABCD₂[n+1:end,1:n], ABCD₂[n+1:end, n+1:end])

    # Implemented Park transformation: d-axis alignment, q-axis lagging, peak-invariant
    T = [1 -1im; -1im -1]
    CK = (2/3)*[1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2]
    CKinv = [1 0; -1/2 sqrt(3)/2; -1/2 -sqrt(3)/2]

    a_dq = T * CK * a₂ * CKinv * conj(T) + conj(T) * CK * a₁ * CKinv * T 
    b_dq = T * CK * b₂ * CKinv * conj(T) + conj(T) * CK * b₁ * CKinv * T
    c_dq = T * CK * c₂ * CKinv * conj(T) + conj(T) * CK * c₁ * CKinv * T
    d_dq = T * CK * d₂ * CKinv * conj(T) + conj(T) * CK * d₁ * CKinv * T

    push!(ABCD_OHLAC3g3_DQ, [a_dq b_dq; c_dq d_dq])

end

Y_BUS_AC2 = []

for i in 1:length(omega)
    ABCD23 = ABCD_CableAC23_DQ[i]
    ABCD2g2 = ABCD_CableAC2g2_DQ[i]
    ABCD3g3 = ABCD_OHLAC3g3_DQ[i]

    #Section 2.1 tutorial
    ABCD2g3g = ABCD2g2*ABCD23*ABCD3g3
    n = Int(size(ABCD2g3g, 1)/2)
    (a2g3g, b2g3g, c2g3g, d2g3g) = (ABCD2g3g[1:n,1:n], ABCD2g3g[1:n,n+1:end], ABCD2g3g[n+1:end,1:n], ABCD2g3g[n+1:end, n+1:end])
    push!(Y_BUS_AC2, [inv(b2g3g)*a2g3g -inv(b2g3g); -inv(b2g3g) inv(b2g3g)*a2g3g])
end

# AC/DC bus admittance matrix

# I = [idg2 iqg2 idg3 iqg3]^T 
# V = [vdg2 vqg2 vdg3 vqg3]^T
# I = Y_BUS*V

Ye = Y_BUS_AC2

L = inv.(Ye) .* Yn
nyquist_Energy_Hub = nyquistplot(L, omega, zoom = "no", SM = "VM")
savefig("nyquist_Energy_Hub.png")

Zcl_bus = inv.(Ye + Yn)

fmin = 1
fmax = 10000

EVD_Energy_Hu1b = EVD(Zcl_bus, omega, fmin, fmax)
savefig("EVD_Energy_Hub.png")

Λₒ = eigvals(net.elements[:OWF].element_value.A)
for i in 1:length(Λₒ)
    if real(Λₒ[i]) >= 0
        println("OWF RHP pole : ", Λₒ[i])
    end
end
