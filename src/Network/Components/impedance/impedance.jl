export impedance

@with_kw mutable struct Impedance
    value :: Array{Basic} = []      # impedance value
    ABCD :: Array{Basic} = []
end


"""
    impedance(;z :: Union{Int, Float64, Basic, Array{Basic}} = 0, pins :: Int = 0)
Creates impedance with specified number of input/output pins `pins`. The impedance expression
 `exp` has to be given in Ω and can have both numerical and symbolic value (example: `z = s-2`).

Pins are named: `1.1`, `1.2`, ..., `1.pins` and `2.1`, `2.2`, ..., `2.pins`

In the case of 1×1 impedance, parameter `z` has only one value.
Example: `impedance(z = 1000, pins = 1)`

If the impedance is multiport, then its value
is given as an array `[val]` with one, `pins` or `pins × pins` number of elements. In the case of one element,
then the impedance has only diagonal nonzero values equal to `val`. In the case of `pins` number of elements,
the impedance has only diagonal nonzero values equal to the values in an array `val`. If the array has the size
`pin × pin`, impedance matrix has the size `pin × pin` and all its values defined.
Examples:
```julia
impedance(z = [s], pins = 3)    # 3×3 impedance with diagonal values equal s
impedance(z = [2,s,s/2], pins = 3) # 3×3 impedance with diagonal values equal 2, s, 0.5s, respectively
impedance(z = [1,s,3,4], pins = 2) # 2×2 impedance with all values defined
```

"""
function impedance(;z :: Union{Int, Float64, Basic, Array{Basic}} = 0, pins :: Int = 0,
        transformation = false)
    if !isempty(z)
        if (pins != 0)
            if (length(z) === pins)
                (pins == 1) ? z = convert(Array{Basic}, Diagonal([z])) : z = convert(Array{Basic}, Diagonal(z))
            elseif (length(z) === 1)
                z = convert(Array{Basic}, Diagonal([z[1], z[1], z[1]]))
            elseif length(z) === pins*pins
                z = reshape(z,pins,pins)
            else
                error("invalid element specification
                ⥤ number of impedance parameters must be 1, $pins or $(pins*pins)")
            end
        else
            pins = Int(sqrt(length(z)))
            z = reshape(z,pins,pins)
        end
        imp = Impedance(value = z)

        # determine ABCD
        m1 = zeros(Basic, 2pins, 2pins)
        m2 = zeros(Basic, 2pins, 2pins)

        # from the matrices according to MNA
        for i in 1:pins
          for j in 1:pins
            z = imp.value[i,j]
            if (z != 0)
              m1[i,i] += 1/z
              m1[pins + j, i] -= 1/z
              m2[i,j] += 1/z
              m2[pins + j, j] -= 1/z
            end
            m1[i, pins + i] = -1
            m2[pins + i, pins + i] = -1
          end
        end
        imp.ABCD = inv(m1)*m2

        element = Element(element_value = imp, input_pins = pins, output_pins = pins,
            transformation = transformation)
    end
    element
end


function eval_abcd(imp :: Impedance, s :: Complex)
    abcd = N.(subs.(imp.ABCD, symbols(:s), s))
    abcd = convert(Array{Float64}, real(abcd)) + 1im*convert(Array{Float64}, imag(abcd))
    abcd = convert(Array{Complex}, abcd)
    return abcd
end

function eval_y(imp :: Impedance, s :: Complex)
    return abcd_to_y(eval_abcd(imp, s))
end

# POWER FLOW
function make_power_flow_ac!(imp :: Impedance , dict :: Dict{String, Any},
            global_dict :: Dict{String, Any})
    key = string(length(dict["branch"]))
    ((dict["branch"])[string(key)])["transformer"] = false
    ((dict["branch"])[string(key)])["tap"] = 1
    ((dict["branch"])[string(key)])["shift"] = 0
    ((dict["branch"])[string(key)])["c_rating_a"] = 1

    abcd = eval_abcd(imp, global_dict["omega"] * 1im)
    n = 3
    Z = (abcd[1:n,n+1:end])[1,1] / global_dict["Z"]

    ((dict["branch"])[string(key)])["br_r"] = real(Z)
    ((dict["branch"])[string(key)])["br_x"] = imag(Z)
    ((dict["branch"])[string(key)])["g_fr"] = 0
    ((dict["branch"])[string(key)])["b_fr"] = 0
    ((dict["branch"])[string(key)])["g_to"] = 0
    ((dict["branch"])[string(key)])["b_to"] = 0

end

function make_power_flow_dc!(imp :: Impedance, dict :: Dict{String, Any},
                        global_dict :: Dict{String, Any})
    key = length(dict["branchdc"])
    ((dict["branchdc"])[string(key)])["l"] = 0
    ((dict["branchdc"])[string(key)])["c"] = 0

    abcd = eval_abcd(imp, 1e-6*1im)
    Z = abcd[1,2] / global_dict["Z"]
    ((dict["branchdc"])[string(key)])["r"] = real(Z)
end
