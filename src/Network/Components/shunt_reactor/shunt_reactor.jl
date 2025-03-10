export shunt_reactor

@with_kw mutable struct Shunt_reactor
    ABCD :: Union{Array{Any}, Array{Basic}}  = []
    organization :: Symbol                   = :Y
    pins  :: Int                             = 1
    N     :: Int                             = 1
    L     :: Union{Int, Float64, Array{Any}, Array{Int}, Array{Float64}} = 1e-3
    R     :: Union{Int, Float64, Array{Any}, Array{Int}, Array{Float64}} = 1e-3
    C_CO  :: Union{Int, Float64, Array{Any}, Array{Int}, Array{Float64}} = 10e-6
    C_IL  :: Union{Int, Float64, Array{Any}, Array{Int}, Array{Float64}} = []
    C_1_E :: Union{Int, Float64}             = 10e-6
end

"""
```julia
organization :: Symbol                                   = :Y     # three-phase organization (:Y, :Yg or :Δ)
pins  :: Int                                             = 1      # marks single or three phase
N     :: Int                                             = 1,     # number of layers
L     :: Union{Int, Float64, Array{Any, Int or Float64}} = 1e-3,  # inductance [H]
R     :: Union{Int, Float64, Array{Any, Int or Float64}} = 1e-3,  # resistance [ohm]
C_CO  :: Union{Int, Float64, Array{Any, Int or Float64}} = 10e-6, # cross-over capacitance [F]
C_IL  :: Union{Int, Float64, Array{Any, Int or Float64}} = 10e-6, # inter-layer capacitance [F]
C_1_E :: Union{Int, Float64}                             = 10e-6  # capacitance between layer 1 and the earthed screen [F]
transformation :: Bool                                   = false  # set to true to transform from abc to dq components
```
Notes:
- The first layer is the innermost layer, the Nth layer is the outermost layer;

L:
- scalar, it is the *total* inductance of all the layers in series;
- array, it must be of length N and each element is the inductance of a layer;
    `L = L_tot`

    `L = [L1, L2, ..., LN]`

R:
- scalar, it is the *total* inductance of all the layers in series;
- array, it must be of length N and each element is the resistance of a layer;

    `R = R_tot`

    `R = [R1, R2, ..., RN]`

C_C0:
- scalar, it is the *average* cross-over capacitance;
- array, it must be of length N and each element is the cross-over capacitance of a layer;

    `C_C0 = C_C0_avg`

    `C_C0 = [C1, C2, ..., CN]`

C_IL:
- scalar, it is the *average* inter-layer capacitance;
- array, it must be of length N-1 and each element is an inter-layer capacitance;

    `C_C0 = C_C0_avg`

    `C_C0 = [C2_1, C3_2, ..., CN_N-1]`

organization:
- if pins == 1, then this parameter is disregarded;
- if pins == 3, then:

    :Yg ⟶ wye connection with grounded neutral

    :Y  ⟶ wye connection with ungrounded neutral

    :Δ  ⟶ delta connection
"""
function shunt_reactor(;
    organization :: Symbol                   = :Y,
    pins  :: Int                             = 1,
    N     :: Int                             = 1,
    L     :: Union{Int, Float64, Array{Any}, Array{Int}, Array{Float64}} = 1e-3,
    R     :: Union{Int, Float64, Array{Any}, Array{Int}, Array{Float64}} = 1e-3,
    C_CO  :: Union{Int, Float64, Array{Any}, Array{Int}, Array{Float64}} = 10e-6,
    C_IL  :: Union{Int, Float64, Array{Any}, Array{Int}, Array{Float64}} = [],
    C_1_E :: Union{Int, Float64}             = 10e-6,
    transformation :: Bool = false)

    # Checking the input parameters
    if pins != 1 && pins != 3
        error("Parameter 'pins' must either be equal to 1 or 3. Here, pins = $pins.")
    end
    if N < 1
        error("Parameter 'N' must either be a stricly positive integer. Here, N = $N.")
    end
    len_L = length(L);
    if len_L == 1 && N > 1
        L = L/N * ones(N);
    elseif len_L != N
        error("Parameter 'L' must either be a scalar or have length 'N'. Here, 'L' has length $len_L but N = $N.")
    elseif len_L > 1
        L = vec(L);
    end
    len_R = length(R);
    if len_R == 1 && N > 1
        R = R/N * ones(N);
    elseif len_R != N
        error("Parameter 'R' must either be a scalar or have length 'N'. Here, 'R' has length $len_R but N = $N.")
    elseif len_R > 1
        R = vec(R);
    end
    len_C_CO = length(C_CO);
    if len_C_CO == 1 && N > 1
        C_CO = C_CO * ones(N);
    elseif len_C_CO != N
        error("Parameter 'C_CO' must either be a scalar or have length 'N'. Here, 'C_CO' has length $len_C_CO but N = $N.")
    elseif len_C_CO > 1
        C_CO = vec(C_CO);
    end
    len_C_IL = length(C_IL);
    if len_C_IL == 1 && N > 2
        C_IL = C_IL * ones(N-1);
    elseif len_C_IL != N-1
        error("Parameter 'C_IL' must either be a scalar or have length 'N-1'. Here, 'C_IL' has length $len_C_IL but N-1 = $(N-1).")
    elseif len_C_IL > 1
        C_IL = vec(C_IL);
    end

    sh       = Shunt_reactor()
    sh.organization = organization;
    sh.pins  = pins;
    sh.N     = N;
    sh.L     = L;
    sh.R     = R;
    sh.C_CO  = C_CO;
    sh.C_IL  = C_IL;
    sh.C_1_E = C_1_E;
    element  = Element(element_value = sh, input_pins = pins, output_pins = pins, transformation = transformation)
    return element
end

function eval_abcd(sh :: Shunt_reactor, s :: Complex)

    """ 0 -- Retrieving the data from the element """

    organization = sh.organization;
    pins  = sh.pins;
    N     = sh.N;
    L     = sh.L;
    R     = sh.R;
    C_CO  = sh.C_CO;
    C_IL  = sh.C_IL;
    C_1_E = sh.C_1_E;

    """ 1 -- build 2N-by-2N ABCD matrix, with all layers disconnected """

    Id = Matrix(1.0I, N, N); # identity matrix
    Ze = zeros(N,N); # zero matrix

    # Building the impedance matrix of (R-L)//C_CO elements:
    if N > 1
        Z_RLC = diagm(0 => -1 ./ (1 ./ (s*L + R) + s*C_CO));
    else # N = 1, only one layer, cannot use diagm()
        Z_RLC = -1 ./ (1 ./ (s*L + R) + s*C_CO);
    end

    # Building the admittance matrix of p-side and q-side IL capacitances:
    if N > 2
        Y_Cp = 1/2*s* diagm(0 => vcat(-C_1_E, -C_IL) + vcat(-C_IL, 0), 1 => C_IL, -1 => C_IL);
    elseif N == 2 # only one inter-layer capacitance, cannot use diagm()
        Y_Cp = 1/2*s* [-C_1_E-C_IL C_IL; C_IL -C_IL];
    else # N = 1, no inter-layer capacitance
        Y_Cp = 1/2*s* -C_1_E;
    end
    Y_Cq = Y_Cp; # q-side is the same as p-side

    # Building the ABCD matrices from the admittances and impedances
    ABCD_Cp  = [Id   Ze;
                Y_Cp Id];
    ABCD_Cq  = [Id   Ze;
                Y_Cq Id];
    ABCD_RLC = [Id Z_RLC;
                Ze   Id];
    ABCD = ABCD_Cq * ABCD_RLC * ABCD_Cp; # this is the 2N-by-2N ABCD matrix with all layers disconnected

    """ 2 -- transform ABCD into Y matrix before applying boundary conditions """

    A = ABCD[1:N,    1:N   ];
    B = ABCD[1:N,    N+1:2N];
    C = ABCD[N+1:2N, 1:N   ];
    D = ABCD[N+1:2N, N+1:2N];

    Y12 = inv(B); # B is easy to invert because it is diagonal, using inv() is ok.
    Y11 = -Y12*A;
    Y22 = D*Y12;
    Y21 = C-D*Y12*A;

    Y = [Y11 Y12; Y21 Y22];

    if N > 1 # if N == 1, Y is already in the reduced form
        """ 3 -- determining which layer connections have to be made"""

        if mod(N, 2) == 0
            p_connections = [i for i in 1:N];
            q_connections = [i for i in 2:N-1];
        else
            p_connections = [i for i in 1:N-1];
            q_connections = [i for i in 2:N];
        end
        p_shape = (2, Int(length(p_connections)/2));
        q_shape = (2, Int(length(q_connections)/2));
        p_map = transpose(reshape(p_connections, p_shape));
        q_map = transpose(reshape(q_connections, q_shape));
        map = vcat(p_map, N .+ q_map);

        """
        In the admittance representation Y, we have:
        input vector : [U1p U2p ... UNp U1q U2q ... UNq]
        output vector: [I1p I2p ... INp I1q I2q ... INq]

        Each row of 'map' gives the index of two pins that have same voltage and
        equal but opposite currents. Consequently, the pins are numbered as follows:
        [1p 2p ... Np 1q 2q ... Nq], so that the notation is valid for both voltage
        and current vectors.
        """

        """ 4 -- Applying boundary conditions for voltage equalities and zero currents sums """

        Y[:, map[:,1]] = Y[:, map[:,1]] + Y[:, map[:,2]];
        Y[map[:,1], :] = Y[map[:,1], :] + Y[map[:,2], :];

        # determining which rows and columns must be kept
        keep = [true for i in 1:2N];
        for i in 1:2N
            if i in map[:,2]
                keep[i] = false;
            end
        end
        Y = Y[keep, keep]; # has now size (N+1)-by-(N+1)

        """ 5 -- Exchanging intermediate voltages from the input vector with zero currents from the output vector """

        # determining which indices correspond to non-relevant voltages (the same indices correspond to zero currents)
        if mod(N, 2) == 0
            exchange = vcat([Int(i) for i in 1:(N/2)], [Int(i) for i in (N/2+2):N]); # indices N/2+1 and N+1 must be preserved
        else
            exchange = vcat([Int(i) for i in 1:((N-1)/2)], [Int(i) for i in ((N+1)/2+2):(N+1)]); # indices (N-1)/2+1 and (N+1)/2+1 must be preserved
        end

        # applying the operation for each index in exchange
        n = size(Y,1);
        for x in exchange
            Y_new = 0 .* Y;
            for r in 1:n
                for c in 1:n
                    if r == x && c == x
                        Y_new[r,c] = 1/Y[x,x];
                    elseif r == x && c != x
                        Y_new[r,c] = -Y[x,c]/Y[x,x];
                    elseif r != x && c == x
                        Y_new[r,c] = Y[r,x]/Y[x,x];
                    else # r ~= x && c ~= x
                        Y_new[r,c] = Y[r,c] - Y[x,c]*Y[r,x]/Y[x,x];
                    end
                end
            end
            Y = Y_new;
        end

        """ 6 -- Removing all zero currents by deleting columns and all intermediate voltages by deleting rows """

        keep = [true for i in 1:n];
        for i in 1:n
            if i in exchange
                keep[i] = false;
            end
        end
        Y = Y[keep, keep];

        """ 7 -- Rearranging into proper ABCD, Y and Z representations """

        if mod(N, 2) == 0
            # current Y representation    final Y representation
            # [Iq1] = [Y11 Y1N][Uq1] ---> [IqN] = [YNN YN1][UqN]
            # [IqN]   [YN1 YNN][UqN]      [Iq1]   [Y1N Y11][Uq1]

            Y_new = 0 .* Y;
            Y_new[1,1] = Y[2,2];
            Y_new[1,2] = Y[2,1];
            Y_new[2,1] = Y[1,2];
            Y_new[2,2] = Y[1,1];
            Y = Y_new;
        end
    end


    """ 8 -- Building the ABCD matrix for the different connection types """
    Y11 = Y[1,1];
    Y12 = Y[1,2];

    Id = Matrix(1.0I, pins, pins); # identity matrix
    Ze = zeros(pins, pins); # zero matrix

    if pins == 1
        Y_sub = -Y11;

    elseif pins == 3
        if organization == :Δ # delta connection
            Y_sub = [(Y12-Y11) -Y12 Y11;
                     Y11 (Y12-Y11) -Y12;
                     -Y12 Y11 (Y12-Y11)];

        elseif organization == :Yg # wye connection with grounded neutral
            Y_sub = [-Y11 0 0;
                     0 -Y11 0;
                     0 0 -Y11];

        elseif organization == :Y # wye connection with ungrounded neutral
            Y_sub = 1/3*Y11.*[-2 1 1;
                              1 -2 1;
                              1 1 -2];

        else
            error("Parameter 'organization' must either be Y, Yg or Δ. Here, organization = $organization.")
        end
    end

    ABCD_final = [Id    Ze;
                  Y_sub Id];

    """ 9 -- Returning the ABCD matrix """

    return ABCD_final
end

function eval_y(sh :: Shunt_reactor, s :: Complex)
    return abcd_to_y(eval_abcd(sh, s))
end

# POWER FLOW
function make_power_flow_ac!(shunt :: Shunt_reactor, dict :: Dict{String, Any},
         global_dict :: Dict{String, Any})

    key = length(dict["shunt"])
    abcd = eval_abcd(shunt, global_dict["omega"] * 1im)
    n = 3
    Y = (abcd[n+1:end,1:n])[1,1] * global_dict["Z"]
    dict["shunt"][string(key)]["gs"] = real(Y)
    dict["shunt"][string(key)]["bs"] = imag(Y)
end

function make_power_flow_dc!(shunt :: Shunt_reactor, dict :: Dict{String, Any},
                         global_dict :: Dict{String, Any})
    nothing
end
