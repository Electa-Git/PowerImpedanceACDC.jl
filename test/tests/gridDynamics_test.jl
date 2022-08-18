using SymEngine
s = symbols("s")
@time net = @network begin
    
    gen1 = ac_source(pins = 3, transformation = true)

    z = impedance(z=3.16e-4*s+0.00995, pins = 3, transformation = true)

    # gen1[1.1] == gen1[1.2] == gen1[1.3] ==gnd
    # gen1[2.1] == z[1.1] == TestNode1
    # gen1[2.2] == z[1.2] == TestNode2
    # gen1[2.3] == z[1.3] == TestNode3
    # z[2.1] == z[2.2] == z[2.3] == gnd

    gen1[2.1] == gnd
    gen1[2.2] == gnd
    gen1[1.1] == z[1.1] == TestNode1
    gen1[1.2] == z[1.2] == TestNode2
    z[2.1] == gnd1
    z[2.2] == gnd2
end

#@time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:gen1], input_pins=Any[:TestNode1,:TestNode2,:TestNode3], 
#output_pins=Any[:gnd,:gnd,:gnd], omega_range = (0,4,1000))

@time imp_ac, omega_ac = determine_impedance(net, elim_elements=[:gen1], input_pins=Any[:TestNode1,:TestNode2], 
output_pins=Any[:gnd1,:gnd2], omega_range = (0,4,1000))

# imp_ac_d = []
# imp_ac_q = []

# for n in eachindex(imp_ac)
#     push!(imp_ac_d,imp_ac[n,1][1])
#     push!(imp_ac_q,imp_ac[n,1][2])
# end
#p = bode(imp_ac, omega = omega_ac)

#save_plot(p)
