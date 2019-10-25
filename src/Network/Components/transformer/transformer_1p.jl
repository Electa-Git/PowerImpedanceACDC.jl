export transformer_1p

"""
    transformer_1p(exp)
Creates a dc or a single phase transfromer.
`Exp` is a set of parameters: `Lₚ`, `Rₚ`, `Rₛ`, `Lₛ`, `Lₘ`, `Rₘ`, `Cₜ`, `C_stray`, `Lₗ`, `Rₗ`, `n`
or it can be written in an impedance form as: `Zᵖ_winding`, `Zˢ_winding`, `Y_turn`,
`Z_load`, `Y_iron`, `Z_stray`.

Pins: `1.1`, `2.1`
"""
function transformer_1p(;Lₚ = 0, Rₚ = 0, Rₛ = 0, Lₛ = 0,
    Lₘ = 0, Rₘ = 0, Cₜ = 0, C_stray = 0, Lₗ = 0, Rₗ = 0, n = 1,
    Zᵖ_winding = convert(Array{Basic},Diagonal([1 for i in 1:2])),
    Zˢ_winding = convert(Array{Basic},Diagonal([1 for i in 1:2])),
    Y_turn = convert(Array{Basic},Diagonal([1 for i in 1:2])),
    Z_load = convert(Array{Basic},Diagonal([1 for i in 1:2])),
    Y_iron = convert(Array{Basic},Diagonal([1 for i in 1:2])),
    Z_stray = convert(Array{Basic},Diagonal([1 for i in 1:2])))

    # definision of elements
    Zᵖ_winding[1,2] += s*Lₚ + Rₚ
    Zˢ_winding[1,2] += s*Lₛ + Rₛ
    Y_turn[2,1] += s*Cₜ
    Z_load[1,2] += s*Lₗ + Rₗ

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
    elem = Element(number_input_pins = 1, number_output_pins = 1, element_value = Z)
    abcd = ABCD_port_representation(elem)

    setfield!(elem.ABCD, :transfer_function, Z)
    elem
end
