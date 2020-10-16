export impedance

@with_kw mutable struct Impedance
    value :: Array{Basic} = []      # impedance value #default value for array
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
        transformation = false)  #z define the value of the impedance pins define the n of pin -> integer
    if !isempty(z) #if z has a value -> enter
        if (pins != 0) #if pins has a value different from 0 -> enter
            if (length(z) === pins)  #if the n of element in z is = to n of pin -> enter
                (pins == 1) ? z = convert(Array{Basic}, Diagonal([z])) : z = convert(Array{Basic}, Diagonal(z)) #ternary operator
                #ternary operator: if pins is =1 do the 1st otherwise to the option after :
            elseif (length(z) === 1) #if z has just 1 value
                z = convert(Array{Basic}, Diagonal([z[1], z[1], z[1]])) #Diagonal-> construct a matrix with diagonal values the one inserted
                #convert(T,x) -> Convert x to a value of type T
            elseif length(z) === pins*pins #if the length of z is equal to the number pins*pins->enter
                z = reshape(z,pins,pins) #construct a square matrix z with the value in z with dimension pins*pins
            else
                error("invalid element specification
                ⥤ number of impedance parameters must be 1, $pins or $(pins*pins)")
                #prints an error message with the statement and the number corresponding to
                #pins or pins*pins
            end
        else
            pins = Int(sqrt(length(z))) #
            #What if not integer sqrt(length(z))??.
            z = reshape(z,pins,pins) #construct a square matrix with the value in z with dimension pins*pins
        end
        imp = Impedance(value = z) #put in the function impedance in value the Z obtained in these passages

        # determine ABCD
        m1 = zeros(Basic, 2pins, 2pins) #create square matrix m1 with dimension 2pins*2pins data type Basic
        m2 = zeros(Basic, 2pins, 2pins) #same as above

        # from the matrices according to MNA ->FROM PAG;17 SIMULATION TUTORIAL
        for i in 1:pins  #for cycle that goes from 1 to the number that pins has
          for j in 1:pins #nested for loop, e.g. if pins=3 it does : i=1,j=1 enter. i=1; j=2 enter. i=1 j=3 enter. i=2,j=1 enter etc.
            z = imp.value[i,j] #access to the value in impedance at the position (i,j)
            if (z != 0) #if the impedance matrix in position (i,j is not void it enters)
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
