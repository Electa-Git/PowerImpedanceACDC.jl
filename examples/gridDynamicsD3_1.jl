using DelimitedFiles,SymEngine
s = symbols("s")
@time net = @network begin

    g1 = ac_source(V = 380 * sqrt(2/3), pins = 3, transformation = false)

    # OHLx2 between AVLGM380 and HORTA380
    ohlAH = overhead_line(length = 39.7e3,
    conductors = Conductors(organization = :vertical,
                nᵇ = 3, Rᵈᶜ = 0.05, rᶜ = 0.016,  yᵇᶜ = 30, Δxᵇᶜ = 14, Δyᵇᶜ= 9, dˢᵃᵍ = 10),
    groundwires = Groundwires(nᵍ = 1, Rᵍᵈᶜ = 0.1, rᵍ = 0.011, Δyᵍ = 27, dᵍˢᵃᵍ = 10),
    earth_parameters = (1,1,25), transformation = false)

    load = impedance(z = s, pins = 3, transformation = false)

    # ohlAH = overhead_line(length = 39.7e3,
    #         conductors = Conductors(organization = :flat, nᵇ = 3, nˢᵇ = 1, Rᵈᶜ = 0.063, rᶜ = 0.015,  yᵇᶜ = 30,
    #                         Δyᵇᶜ = 0, Δxᵇᶜ = 10,  Δ̃xᵇᶜ = 0, dˢᵇ = 0,  dˢᵃᵍ = 10),
    #         groundwires = Groundwires(nᵍ = 2, Rᵍᵈᶜ = 0.92, rᵍ = 0.0062, Δxᵍ = 6.5, Δyᵍ = 7.5, dᵍˢᵃᵍ   = 10),
    #         earth_parameters = (1,1,100), transformation = false)

    g1[2.1] == gnda
    g1[2.2] == gndb
    g1[2.3] == gndc

    g1[1.1] == ohlAH[1.1] == Bus1a
    g1[1.2] == ohlAH[1.2] == Bus1b
    g1[1.3] == ohlAH[1.3] == Bus1c

    # ohlAH[2.1] == gnda
    # ohlAH[2.2] == gndb
    # ohlAH[2.3] == gndc

    ohlAH[2.1] == load[1.1]
    ohlAH[2.2] == load[1.2]
    ohlAH[2.3] == load[1.3]

    load[2.1] == gnda
    load[2.2] == gndb
    load[2.3] == gndc

    # g1[1.1] == ohlAH[1.1] == ohlAH[1.4] == Bus1a
    # g1[1.2] == ohlAH[1.2] == ohlAH[1.5] == Bus1b
    # g1[1.3] == ohlAH[1.3] == ohlAH[1.6] == Bus1c

    # ohlAH[2.1] == ohlAH[2.4] == gnda
    # ohlAH[2.2] == ohlAH[2.5] == gndb
    # ohlAH[2.3] == ohlAH[2.6] == gndc

   

end

@time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:g1], input_pins=Any[:Bus1a,:Bus1b,:Bus1c], 
output_pins=Any[:gnda,:gndb,:gndc], omega_range = (0,4,2000))
p = bode(imp_ac, omega = omega_ac)

writedlm("imp_abc_D31.csv",  imp_ac, ',')
writedlm("w_D31.csv",  omega_ac, ',')

# @time net_dq = @network begin


#     g1 = ac_source(V = 380 * sqrt(2/3), pins = 3, transformation = true)

#     #OHLx2 between AVLGM380 and HORTA380
#     ohlAH = overhead_line(length = 39.7e3,
#     conductors = Conductors(organization = :vertical,
#                 nᵇ = 3, Rᵈᶜ = 0.05, rᶜ = 0.016,  yᵇᶜ = 30, Δxᵇᶜ = 14, Δyᵇᶜ= 9, dˢᵃᵍ = 10),
#     groundwires = Groundwires(nᵍ = 1, Rᵍᵈᶜ = 0.1, rᵍ = 0.011, Δyᵍ = 27, dᵍˢᵃᵍ = 10),
#     earth_parameters = (1,1,25), transformation = true)

#     g1[2.1] == gndd
#     g1[2.2] == gndq

#     g1[1.1] == ohlAH[1.1] == Bus1d
#     g1[1.2] == ohlAH[1.2] == Bus1q

#     ohlAH[2.1] ==  gndd
#     ohlAH[2.2] ==  gndq

# end

# @time imp_ac_dq, omega_ac_dq = determine_impedance(net_dq, elim_elements=[:g1], input_pins=Any[:Bus1d,:Bus1q], 
# output_pins=Any[:gndd,:gndq], omega_range = (0,4,2000))

# writedlm("imp_dq_D31.csv",  imp_ac_dq, ',')
# writedlm("w_dq_D31.csv",  omega_ac_dq, ',')



