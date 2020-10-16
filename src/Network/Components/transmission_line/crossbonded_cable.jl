export crossbonded_cable #export the crossbonded_cable function (first letter lowercase-> function)

@with_kw mutable struct Crossbonded_cable <: Transmission_line #generate the mutable structure Crossbonded_cable ()
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
cable implementation.

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
function crossbonded_cable(;args...) #variable arguments function -> "varargs" functions, they take an arbitrary number of arguments (Dictionary cable definition)
    cable = Crossbonded_cable() #zero-arguments function TBC better
    pins = 3 #definition of the number of pins for a crossbonded_cable
    transformation = false #transformation (Park) false as default
    for (key, val) in kwargs_pairs(args)  #this cycle goes through all the tuple of keywords and values (key, val) in the "variable" args
        if in(key, propertynames(cable)) #if key is in propertynames(cable) -> enter. propertynames -> report public properties of variable cable
            setfield!(cable, key, val) #Assign val to key of the composite type cable
        elseif isa(val.element_value, Cable) #determin if val.element_value is of the type Cable. If yes-> enter
            val.element_value.eliminate = false #from cable mutable struct in cable.jl
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

    order = [i for i in 1:2n*nₗ] #in case of n=3 and nₗ=2 -> order= vector of 12 elements from 1 to 12
    for i in 1:nₗ
        order[1+(i-1)n:i*n] = [i + (j-1)nₗ for j in 1:n] #if n=3 and nₗ=2 -> order = [1 3 5 2 4 6 7 8 9 10 11 12] (first all core voltages then all sheats voltages, currents still not sorted)
    end
    order[1+n*nₗ:2n*nₗ] = order[1:n*nₗ] + n*nₗ*ones(Int, n*nₗ, 1)# if n=3 and nₗ=12 -> order = [1 3 5 2 4 6 7 9 11 8 10 12] (fist all core voltages then all sheats voltages, core currents then sheats currents)

    Mᶜᵇ = convert(Array{Complex}, Diagonal([1 for i in 1:(2n*nₗ)])) #for n=3 nₗ=2 Mᶜᵇ= [12x12] matrix
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
    abcd = kron_abcd(abcd, no_eliminate)#Reduction of the ABCD matrix explained in section 2.6??

    total_abcd = convert(Array{Complex}, Diagonal([1 for dummy in 1:2n]))
    for i in 1:m.mₛ #i that goes from 1 to mₛ. Where mₛ is the number of major sections
        total_abcd *= abcd #Computation of the equivalent cable ABCD parameters -> bottom part pag 25 simulator_tutorial
    end

    return total_abcd
end

function eval_y(m :: Crossbonded_cable, s :: Complex)
    abcd_to_y(eval_abcd(m,s))
end
