# include("C:/Users/osakinci/Documents/GitHub/HVDCstability.jl/src/HVDCstability.jl")
using HVDCstability, SymEngine, DelimitedFiles
s = symbols("s")
@time net = @network begin
    
    gen1 = ac_source(V = 220 * sqrt(2/3), pins = 3, transformation = true)
    
    # tl78 = overhead_line(length = 35e3,
    #             conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
    #                             Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
    #             groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
    #             earth_parameters = (1,1,100), transformation = true)

    # gen1[2.1] == gndd
    # gen1[2.2] == gndq
    # gen1[1.1] == tl78[1.1] == TestNode1d
    # gen1[1.2] == tl78[1.2] == TestNode1q
    # tl78[2.1] == gndd
    # tl78[2.2] == gndq

    sg = synchronousmachine(P = 49.9)

    gen1[2.1] == gndd
    gen1[2.2] == gndq
    gen1[1.1] == sg[1.1] == TestNode1d
    gen1[1.2] == sg[1.2] == TestNode1q
    sg[2.1] == gndd
    sg[2.2] == gndq


    # c1 = mmc(Vᵈᶜ = 640, Vₘ = 380*sqrt(2/3),
    #             P_max = 1000, P_min = -1000, P = 124, Q = -1, Q_max = 500, Q_min = -500, P_dc = -124,
    #             occ = PI_control(Kₚ = 135.0011, Kᵢ = 9.1603e+04),
    #             ccc = PI_control(Kₚ = 42.9123, Kᵢ = 1.9739e+04),
    #             pll = PI_control(Kₚ = 2.8351e-04, Kᵢ = 0.0127),
    #             power = PI_control(Kₚ = 2.1487e-07, Kᵢ = 6.7503e-05)
    #             )
    # dc = dc_source(V = 640)

    # dc[1.1] == gnd
        
    # dc[2.1] == c1[1.1]

    # gen1[2.1] == gnd1d
    # gen1[2.2] == gnd2q
    # c1[2.1] == gen1[1.1] == sg[1.1] == TestNode1
    # c1[2.2] == gen1[1.2] == sg[1.2] == TestNode2
    # sg[2.1] == gnd1d
    # sg[2.2] == gnd2q


end

@time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:gen1], input_pins=Any[:TestNode1d,:TestNode1q], 
output_pins=Any[:gndd,:gndq], omega_range = (-2,4,2000))

p = bode(imp_ac, omega = omega_ac)

writedlm("imp_dq_SG.csv",  imp_ac, ',')
writedlm("w_dq_SG.csv",  omega_ac, ',')
