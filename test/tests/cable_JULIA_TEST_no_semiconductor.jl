# include("../../src/HVDCstability.jl")
# using .HVDCstability
using DelimitedFiles
using SymEngine

# semiconducting layers test
net = @network begin
    #ac = "3-phase ac voltage source with V amplitude 380kV
    ac = ac_source(V = 380, pins = 3)
    #c = cable with length 26km
    c = cable(length = 26e3, positions = [(-0.5,1.9), (0,1.9), (0.5,1.9)],
              C1 = Conductor(rₒ = 31.75e-3, ρ = 2.18e-8), #inner radius[m], resistivity [Ωm]
              C2 = Conductor(rᵢ = 57.55e-3, rₒ = 57.75e-3, ρ = 1.72e-8), # inner and outer radius [m], resistivity [Ωm]
              I1 = Insulator(rᵢ = 31.75e-3, rₒ = 57.55e-3, ϵᵣ = 2.26), #inner radius, inner semiconductor outer radius, outer semiconductor inner radius, outer radius, relative permittivity
              I2 = Insulator(rᵢ = 57.75e-3, rₒ = 62.65e-3, ϵᵣ = 2.26))
    ac[1.1] ⟷ c[1.1] ⟷ Node1
    ac[1.2] ⟷ c[1.2] ⟷ Node2
    ac[1.3] ⟷ c[1.3] ⟷ Node3
    ac[2.1] ⟷ ac[2.2] ⟷ ac[2.3] ⟷ gnd
    c[2.1] ⟷ gnd1
    c[2.2] ⟷ gnd2
    c[2.3] ⟷ gnd3

end



imp, omega = determine_impedance(net, elim_elements = [:ac], input_pins = Any[:Node1,:Node2,:Node3],
                            output_pins = Any[:gnd1,:gnd2,:gnd3], omega_range = (0.1,5,500))
bode(imp, omega = omega, axis_type = :loglin)
writedlm( "impauto.csv", imp,',')
writedlm( "omega.csv", omega, ',')
save_data( net.elements[:c], "cable_no_semi", omega_range= (0.1, 5, 500), scale= :log)
