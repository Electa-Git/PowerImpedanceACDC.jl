export crossbonded_cable

@with_kw mutable struct Crossbonded_cable <: Transmission_line
    nₛ :: Int = 1            # number of the minor sections
    mₛ :: Int = 1            # number of the major sections
    Zᶜᵇ :: Union{Int, Float64, Basic} = 0 # cross-bonding impedance
    sections :: Dict{Symbol, Cable} = Dict{Symbol, Cable}()
end

"""
    function crossbonded_cable(;args...)
This function creates cross bonded three phase cable component. The component consists
of dictionary cable definitions. It creates one component by concatenating cable sections
with applied layers transposition matrix. Cable definitions are the same as for the
cable impementation.

Arguments:
- nₛ - number of minor sections
- mₛ - number of major sections
- Zᶜᵇ - cross-bonding impedance
- sections - consists of cable segments

After providing number of section `nₛ`, the user can provide either `nₛ` cable definitions
or one definition. In case that the one cable definition is provided, it is considered
that all the sections are the same.

Cables should be called using designator `C` with indices 1, 2, 3, ..., nₛ.

Transposition matrix transposes a-b-c to c-a-b.
"""
function crossbonded_cable(;args...)
    cable = Crossbonded_cable()
    pins = 3
    transformation = false
    for (key, val) in kwargs_pairs(args)
        if in(key, propertynames(cable))
            setfield!(cable, key, val)
        elseif isa(val.element_value, Cable)
            val.element_value.eliminate = false
            cable.sections[key] = val.element_value
        elseif (key == :transformation)
            transformation = val
        else
            throw(ArgumentError("Unknown property $(key) of the cross-bonded cable."))
        end
    end

    n = length(filter(p -> occursin("C", string(p)), keys(cable.sections)))
    if (n < cable.nₛ) && (n > 1)
        throw(ArgumentError("Definitions should be given for only one or for
            $(cable.nₛ) minor cable sections."))
    end

    elem = Element(input_pins = pins, output_pins = pins, element_value = cable,
            transformation = transformation)
end

function eval_abcd(m :: Crossbonded_cable, s :: Complex)
    n = length(m.sections[:C1].positions)       # number of cables, it should be 3
    nₗ = length(m.sections[:C1].conductors)      # number of conducting layers, normally 2

    order = [i for i in 1:2n*nₗ]
    for i in 1:nₗ
        order[1+(i-1)n:i*n] = [i + (j-1)nₗ for j in 1:n]
    end
    order[1+n*nₗ:2n*nₗ] = order[1:n*nₗ] + n*nₗ*ones(Int, n*nₗ, 1)

    Mᶜᵇ = convert(Array{Complex}, Diagonal([1 for i in 1:(2n*nₗ)]))
    isa(m.Zᶜᵇ, Basic) ? Z = N(subs(m.Zᶜᵇ, symbols(:s), s)) : Z = m.Zᶜᵇ
    Mᶜᵇ[n+1:2n, n*nₗ+n+1:n*nₗ+2n] = Diagonal([2Z for i in 1:n])

    t_order = [i for i in 1:2n*nₗ]
    for i in 2:2:4
        t_order[(i-1)n+1:i*n] = (t_order[(i-1)n+1:i*n])[[3;1;2]]
    end

    abcd = 0
    for (key, val) in m.sections
        abcd_new = reshape(eval_abcd(m.sections[key], s)[order, :], 2n*nₗ, 2n*nₗ)
        abcd_new = reshape(abcd_new[:, order], 2n*nₗ, 2n*nₗ)
        if (abcd == 0)
            abcd = abcd_new
        else
            abcd_new = abcd_new[t_order, :]
            abcd_new = abcd_new[:, t_order]
            abcd = abcd * Mᶜᵇ * abcd_new
        end
    end
    if length(m.sections) == 1
        abcd_new = abcd[t_order, :]
        abcd_new = abcd[:, t_order]
        abcd_new = Mᶜᵇ * abcd
        for i in 2:m.nₛ
            abcd *= abcd_new
        end
    end
    abcd = convert(Array{Complex}, abcd)

    no_eliminate = [i for i in 1:n]
    abcd = kron_abcd(abcd, no_eliminate)

    total_abcd = convert(Array{Complex}, Diagonal([1 for dummy in 1:2n]))
    for i in 1:m.mₛ
        total_abcd *= abcd
    end

    return total_abcd
end

function eval_y(m :: Crossbonded_cable, s :: Complex)
    abcd_to_y(eval_abcd(m,s))
end
