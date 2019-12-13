""" data from Lei Wu's PhD thesis """

# org   = :Yg;    # three-phase organization (:Y, :Yg or :Δ) --- not necessary for single phase connection
pins  = 1;      # marks single or three phase
N     = 11;     # number of layers
V_nom = 50e3;   # [V]
Q_nom = 100e6;  # [VAr]
f_nom = 50;     # [Hz]
L     = V_nom^2 / (2*π*f_nom*Q_nom);  # inductance [H]
R     = 3.4e-3; # resistance [ohm]
C_CO  = [0.047 0.049 0.051 0.054 0.056 0.058 0.061 0.063 0.065 0.068 0.070];    # cross-over capacitance [F]
C_IL  = [31.1 32.6 34.1 35.6 37.1 38.6 40.1 41.6 43.0 44.5];                    # inter-layer capacitance [F]
C_1_E = 8.6;   # capacitance between layer 1 and the earthed screen [F]
Ro    = 50;    # [Ω], output measurement resistance

""" first possibility: use determine_impedance and invert each impedance value to obtain the admittance """
""" to check, does not seem to work """

# net = @network begin
#     shunt = shunt_reactor(pins=pins, N=N, L=L, R=R, C_CO=C_CO, C_IL=C_IL, C_1_E=C_1_E)
#     vs1 = ac_source(V = 50e3, P = 0, pins = 1) # dummy source
#     vs2 = ac_source(V = 50e3, P = 0, pins = 1) # dummy source
#
#     vs1[2.1] ⟷ gnd1
#     vs1[1.1] ⟷ shunt[1.1] ⟷ Node1
#     vs2[1.1] ⟷ shunt[2.1] ⟷ Node2
#     vs2[2.1] ⟷ gnd2
# end
# imp, omega = determine_impedance(net, input_pins = Any[:Node1], output_pins = Any[:Node2], elim_elements = [:vs1, :vs2], omega_range = (1, 3, 100))
# bode(imp, omega = omega)


""" other possibility: use eval_abcd and extract the admittance """

shunt = shunt_reactor(pins=pins, N=N, L=L, R=R, C_CO=C_CO, C_IL=C_IL, C_1_E=C_1_E)

omega_range = (1, 6, 10)
(min_ω, max_ω, n_ω) = omega_range
n = (max_ω - min_ω) / n_ω
omegas = [exp10(min_ω)*10^(i*n) for i in 1:Int(n_ω)]

Y = []
nb_omegas = length(omegas)
for k in 1:nb_omegas
    print("$k/$nb_omegas\n")
    abcd = eval_abcd(shunt.element_value, 1im*omegas[k])
    # ABCD = [ 1 0]
    #        [-Y 1]
    push!(Y, -1*abcd[pins+1:2pins,1:pins])
end

bode(Ro.*Y, omega = omegas)
