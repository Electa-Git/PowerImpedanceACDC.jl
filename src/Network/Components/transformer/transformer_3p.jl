export transformer_ΔY
export transformer_YY


"""
    transformer_YY(exp)

Creates 3 phase YY transformer with the same parameters as a single phase transformer.

Pins: `1.1`, `1.2`, `1.3`, and `2.1`, `2.2`, `2.3`
"""
function transformer_YY(;Lₚ = 0, Rₚ = 0, Rₛ = 0, Lₛ = 0,
    Lₘ = 0, Rₘ = 0, Cₜ = 0, C_stray = 0, Lₗ = 0, n = 1)

    Zᵖ_winding = convert(Array{Basic},Diagonal([1 for i in 1:2]))
    Zˢ_winding = convert(Array{Basic},Diagonal([1 for i in 1:2]))
    Y_turn = convert(Array{Basic},Diagonal([1 for i in 1:2]))
    Z_load = convert(Array{Basic},Diagonal([1 for i in 1:2]))
    Y_iron = convert(Array{Basic},Diagonal([1 for i in 1:2]))
    Z_stray = convert(Array{Basic},Diagonal([1 for i in 1:2]))

    # definision of elements
    Zᵖ_winding[1,2] += s*Lₚ+Rₚ
    Zˢ_winding[1,2] += s*Lₛ+Rₛ
    Y_turn[2,1] += s*Cₜ

    # Y_iron
    if (Lₘ != 0)
        Y_iron[2,1] += 1/s/Lₘ
    end
    if (Rₘ != 0)
        Y_iron[2,1] += 1/Rₘ
    end
    if (C_stray != 0)
        Z_stray = [1 1/s/C_stray; 0 1]
    end

    N_tr = [n 0; 0 1/n]

    # calculate ABCD matrix
    Z_inner = Zᵖ_winding * Y_iron * N_tr * Zˢ_winding

    if (C_stray != 0)
        Z = connect_parallel!(Z_inner, Z_stray) * Z_load
    else
        Z = Z_inner * Z_load
    end

    Z₃ₚ = zeros(Basic, 6, 6)
    for i in 1:3
        Z₃ₚ[i,i] = Z[1,1]
        Z₃ₚ[i,3+i] = Z[1,2]
        Z₃ₚ[3+i,i] = Z[2,1]
        Z₃ₚ[3+i,3+i] = Z[2,2]
    end

    elem = Element(number_input_pins = 3, number_output_pins = 3,
            element_type= :transformer_YY, element_value = Z₃ₚ)
    abcd = ABCD_port_representation(elem)
    setfield!(elem.ABCD, :transfer_function, Z₃ₚ)
    elem
end

"""
    function transformer_ΔY(exp)
Creates 3 phase YY transformer with the same parameters as a single phase transformer.

Pins: `1.1`, `1.2`, `1.3`, and `2.1`, `2.2`, `2.3`
"""
function transformer_ΔY(;Lₚ = 0, Rₚ = 0, Rₛ = 0, Lₛ = 0,
    Lₘ = 0, Rₘ = 0, Cₜ = 0, C_stray = 0, Lₗ = 0, Rₗ = 0, n = 1)

    Zᵖ_winding = convert(Array{Basic},Diagonal([1 for i in 1:2]))
    Zˢ_winding = convert(Array{Basic},Diagonal([1 for i in 1:6]))
    Y_turn = convert(Array{Basic},Diagonal([1 for i in 1:6]))
    Z_load = convert(Array{Basic},Diagonal([1 for i in 1:6]))
    Y_iron = convert(Array{Basic},Diagonal([1 for i in 1:2]))
    Z_stray = convert(Array{Basic},Diagonal([1 for i in 1:6]))

    # definision of elements
    Zᵖ_winding[1,2] += s*Lₚ+Rₚ
    # Y_iron
    if (Lₘ != 0)
        Y_iron[2,1] += 1/s/Lₘ
    end
    if (Rₘ != 0)
        Y_iron[2,1] += 1/Rₘ
    end
    N_tr = [n 0; 0 1/n]

    for i in 1:3
        Zˢ_winding[i,i+3] += s*Lₛ + Rₛ
        Y_turn[i+3,i] += s*Cₜ
        Z_load[i,i+3] += s*Lₗ + Rₗ
        if (C_stray != 0)
            Z_stray[i,i+3] += 1/s/C_stray
        end
    end


    # calculate ABCD matrix
    (A,B,C,D) = Zᵖ_winding * Y_iron * N_tr
    Z_inner = zeros(Basic,6,6)
    current_conversion = [1/sqrt(3) 0 -1/sqrt(3);
                          -1/sqrt(3) 1/sqrt(3) 0;
                          0 -1/sqrt(3) 1/sqrt(3)] .* D
    for i in 1:3
        Z_inner[i,i] = A * sqrt(3)
        Z_inner[i, 3+i] = B / sqrt(3)
        Z_inner[3+i, i] = C * sqrt(3)
        for j in 1:3
            Z_inner[3+i, 3+j] = current_conversion[i, j]
        end
    end

    #println(Z_inner, Zˢ_winding)
    if (C_stray != 0)
        Z = connect_parallel!(Z_inner, Z_stray) * Z_load
    else
        Z = Z_inner * Z_load
    end

    elem = Element(number_input_pins = 3, number_output_pins = 3,
        element_type= :transformer_ΔY, element_value = Z)
    abcd = ABCD_port_representation(elem)
    setfield!(elem.ABCD, :transfer_function, Z)
    elem
end
