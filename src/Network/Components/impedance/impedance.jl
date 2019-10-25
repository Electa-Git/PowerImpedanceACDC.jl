export impedance

@with_kw mutable struct Impedance
    value :: Array{Basic} = []      # impedance value
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
```
impedance(z = [s], pins = 3)    # 3×3 impedance with diagonal values equal s
impedance(z = [2,s,s/2], pins = 3) # 3×3 impedance with diagonal values equal 2, s, 0.5s, respectively
impedance(z = [1,s,3,4], pins = 2) # 2×2 impedance with all values defined
```

"""
function impedance(;z :: Union{Int, Float64, Basic, Array{Basic}} = 0, pins :: Int = 0)
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
        z = Impedance(value = z)
        element = Element(element_value = z, input_pins = pins, output_pins = pins)
    end
    element
end


function create_abcd!(elem :: Element, imp :: Impedance)
    m1 = zeros(Basic, np(elem), np(elem))
    m2 = zeros(Basic, np(elem), np(elem))

    # from the matrices according to MNA
    for i in 1:nip(elem)
      for j in 1:nop(elem)
        z = imp.value[i,j]
        if (z != 0)
          m1[i,i] += 1/z
          m1[nip(elem) + j, i] -= 1/z
          m2[i,j] += 1/z
          m2[nip(elem) + j, j] -= 1/z
        end
        m1[i, nop(elem) + i] = -1
        m2[nip(elem) + i, nop(elem) + i] = -1
      end
    end
    return m1\m2
end

function eval_abcd(element :: Element, imp :: Impedance, s :: Complex)
    return N.(subs.(element.ABCD, symbols(:s), s))
end
