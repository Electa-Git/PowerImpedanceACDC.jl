export transformer

@with_kw mutable struct Transformer
    ABCD :: Array{Basic} = []          # ABCD value

    pins :: Int = 1                     # marks single or three phase
    organization :: Symbol = :YY        # three phase organization (:YY or :ΔY)

    ω :: Union{Int, Float64} = 2*π*50   # rated frequency in [Hz]
    V₁ᵒ :: Union{Int, Float64} = 0      # open circuit primary voltage [V]
    V₁ˢ :: Union{Int, Float64} = 0      # short circuit primary voltage [V]
    I₁ᵒ :: Union{Int, Float64} = 0      # open circuit primary current [V]
    I₁ˢ :: Union{Int, Float64} = 0      # short circuit primary current [V]
    P₁ᵒ :: Union{Int, Float64} = 0      # open circuit losses on primary side [W]
    P₁ˢ :: Union{Int, Float64} = 0      # short circuit losses on primary side [W]
    V₂ᵒ :: Union{Int, Float64} = 0      # open circuit secondary voltage [V]
    V₂ˢ :: Union{Int, Float64} = 0      # short circuit secondary voltage [V]

    n :: Union{Int, Float64} = 0        # turn ratio
    Lₚ :: Union{Int, Float64} = 0       # primary side inductance [H]
    Rₚ :: Union{Int, Float64} = 0       # primary side resistance [Ω]
    Rₛ :: Union{Int, Float64} = 0       # secondary side resistance [Ω]
    Lₛ :: Union{Int, Float64} = 0       # secondary side inductance [H]
    Lₘ :: Union{Int, Float64} = 0      # magnetising inductance [H]
    Rₘ :: Union{Int, Float64} = 0      # magnetising resistance [Ω]
    Cₜ :: Union{Int, Float64} = 0       # turn-to-turn capacitance [F]
    Cₛ :: Union{Int, Float64} = 0       # stray capacitance [F]

end

"""
    function transformer(;args...)
Creates a dc or a single phase transfromer, or a three-phase transformer in YY
or ΔY configuration.

```julia
@with_kw mutable struct Transformer
    value :: Array{Basic} = []          # ABCD value

    pins :: Int = 1                     # marks single or three phase
    organization :: Symbol = :YY        # three phase organization (:YY or :ΔY)

    ω :: Union{Int, Float64} = 2*π*50   # rated frequency in [Hz]
    V₁ᵒ :: Union{Int, Float64} = 0      # open circuit primary voltage [V]
    V₁ˢ :: Union{Int, Float64} = 0      # short circuit primary voltage [V]
    I₁ᵒ :: Union{Int, Float64} = 0      # open circuit primary current [V]
    I₁ˢ :: Union{Int, Float64} = 0      # short circuit primary current [V]
    P₁ᵒ :: Union{Int, Float64} = 0      # open circuit losses on primary side [W]
    P₁ˢ :: Union{Int, Float64} = 0      # short circuit losses on primary side [W]
    V₂ᵒ :: Union{Int, Float64} = 0      # open circuit secondary voltage [V]
    V₂ˢ :: Union{Int, Float64} = 0      # short circuit secondary voltage [V]

    n :: Union{Int, Float64} = 0        # turn ratio
    Lₚ :: Union{Int, Float64} = 0       # primary side inductance [H]
    Rₚ :: Union{Int, Float64} = 0       # primary side resistance [Ω]
    Rₛ :: Union{Int, Float64} = 0       # secondary side resistance [Ω]
    Lₛ :: Union{Int, Float64} = 0       # secondary side inductance [H]
    Lₘ :: Union{Int, Float64} = 0      # magnetising inductance [H]
    Rₘ :: Union{Int, Float64} = 0      # magnetising resistance [Ω]
    Cₜ :: Union{Int, Float64} = 0       # turn-to-turn capacitance [F]
    Cₛ :: Union{Int, Float64} = 0       # stray capacitance [F]
end
```

Pins: `1.1`, `2.1` for single phase transformer and `1.1`, `1.2`, `1.3`, `2.1`,
`2.2`, `2.3` for a three-phase transformer.
"""
function transformer(;args...)
    t = Transformer()
    for (key, val) in kwargs_pairs(args)
        if in(key, propertynames(t))
            setfield!(t, key, val)
        end
    end

    if (t.V₁ᵒ  != 0)
        t.n = t.V₁ᵒ / t.V₂ᵒ

        R = t.P₁ˢ / (t.I₁ˢ)^2
        L = sqrt((t.V₁ˢ*t.I₁ˢ)^2 - t.P₁ˢ^2) / t.ω / (t.I₁ˢ)^2
        t.Rₚ = R / 2
        t.Rₛ = R / 2 / t.n^2
        t.Lₚ = L / 2
        t.Lₛ = L / 2 / t.n^2

        t.Rₘ = (t.V₁ᵒ)^2 / t.P₁ᵒ
        t.Lₘ = (t.V₁ᵒ)^2 / sqrt((t.V₁ᵒ*t.I₁ᵒ)^2 - t.P₁ᵒ^2) / t.ω
    end

    # determine ABCD parameters
    s = symbols(:s)
    Zᵖ_winding = convert(Array{Basic},Diagonal([1 for i in 1:2]))
    Zˢ_winding = convert(Array{Basic},Diagonal([1 for i in 1:2]))
    Y_turn = convert(Array{Basic},Diagonal([1 for i in 1:2]))
    Y_iron = convert(Array{Basic},Diagonal([1 for i in 1:2]))
    Z_stray = convert(Array{Basic},Diagonal([1 for i in 1:2]))

    Zᵖ_winding[1,2] += s*t.Lₚ + t.Rₚ
    Zˢ_winding[1,2] += s*t.Lₛ + t.Rₛ
    Y_turn[2,1] += s*t.Cₜ
    N_tr = [t.n 0; 0 1/t.n]

    # Y_iron
    (t.Lₘ != 0) ? Y_iron[2,1] += 1/s/t.Lₘ : nothing
    (t.Rₘ != 0) ? Y_iron[2,1] += 1/t.Rₘ : nothing
    (t.Cₛ != 0) ? Z_stray = [1 1/s/t.Cₛ; 0 1] : nothing
    Z = Zᵖ_winding * Y_iron * N_tr * Zˢ_winding

    Z = connect_parallel!(Z, Z_stray)

    if (t.pins == 1)
        t.ABCD = Z
    else
        if (t.organization == :YY)
            Z₃ₚ = zeros(Basic, 6, 6)
            for i in 1:3
                Z₃ₚ[i,i] = Z[1,1]
                Z₃ₚ[i,3+i] = Z[1,2]
                Z₃ₚ[3+i,i] = Z[2,1]
                Z₃ₚ[3+i,3+i] = Z[2,2]
            end

            t.ABCD = Z₃ₚ
        else
            Zˢ_winding = convert(Array{Basic},Diagonal([1 for i in 1:6]))
            Y_turn = convert(Array{Basic},Diagonal([1 for i in 1:6]))
            Z_stray = convert(Array{Basic},Diagonal([1 for i in 1:6]))

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

            if (t.Cₛ != 0)
                Z = connect_parallel!(Z_inner, Z_stray)
            else
                Z = Z_inner
            end

            t.ABCD = Z
        end
    end

    elem = Element(input_pins = t.pins, output_pins = t.pins, element_value = t)
end

function eval_abcd(t :: Transformer, s :: Complex)
    return N.(subs.(t.ABCD, symbols(:s), s))
end
