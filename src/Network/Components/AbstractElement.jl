__precompile__()

export closing_impedance, save_data, plot_data

include("plot.jl")

const ABCD_parameters = Union{Array{Basic}, Array{Array{Basic}}}

"""
Struct guarantees representation of the component like a multiport
network using ABCD parameters. It consists of:
- element unique symbol inside constructed network - `symbol`
- dictionary that maps pins inside the network - `pins`
- number of input pins - `input_pins`
- number of output pins - `output_pins`
- component definition - `element_value`
- ABCD multiport parameters - `ABCD`

It generates ABCD matrix that maps:
`output_pins(voltage, current) = [ABCD matrix] × input_pins(voltage, current)``
"""
mutable struct Element
  symbol::Symbol
  pins :: Dict{Symbol, Symbol}
  input_pins :: Int
  output_pins :: Int
  element_value :: Any  # component defined type

  ABCD :: Union{Array{Basic}, Array{Array{Basic}}}

  function Element(;args...)
    elem = new()
    for (key, val) in kwargs_pairs(args)
      if (key in propertynames(elem))
        setfield!(elem, key, val)
      else
        if key === :ports       # add not connected pins
          pins = Dict{Symbol,Symbol}(branch => Symbol() for branch in val)
          val = pins
          key = :pins
          setfield!(elem, key, val)
        end
      end
    end

    # definition of pins
    if !isdefined(elem, :pins)
      elem.pins = merge(Dict{Symbol, Symbol}(Symbol(string("1.",i)) => Symbol() for i in 1:nip(elem)),
                        Dict{Symbol, Symbol}(Symbol(string("2.",i)) => Symbol() for i in 1:nop(elem)))
    end

    if !(isdefined(elem, :ABCD))
        tf = create_abcd!(elem, elem.element_value)
        add!(elem, :ABCD, tf)
    end
    elem
  end
end

for (n,m) in Dict(:nip => :input_pins, :nop => :output_pins)
  @eval ($n)(e::Element) = e.$m
end
np(e::Element) = nip(e) + nop(e) # number pins

function add!(elem::Element, sym::Symbol, value::Any)
  if (sym in propertynames(elem))
    setfield!(elem, sym, value)
  end
end

function get_abcd(element::Element)
    return element.ABCD
end

function get_abcd(element::Element, s::Complex)
    return eval_abcd(element, element.element_value, s)
end

nip_abcd(e::Element) = Int(size(get_abcd(e),1)/2)
nop_abcd(e::Element) = Int(size(get_abcd(e),2)/2)
np_abcd(e::Element) = nip_abcd(e) + nop_abcd(e) # number pins

function eval_abcd(element::Element, dict::Dict{Any,Any})
  s = symbols(:s)
  abcd = eval_abcd(element.element_value, dict[s])
  for i in 1:length(abcd)
      dict[Basic(string(element.symbol,i))] = abcd[i]
  end
end

function connect_serial!(a::ABCD_parameters, b::ABCD_parameters)
    return a*b
end

connect_serial!(a::Element, b::Element) = connect_serial!(a.ABCD, b.ABCD)
connect_serial!(a::Element, b::ABCD_parameters) = connect_serial!(a.ABCD, b)
connect_serial!(a::ABCD_parameters, b::Element) = connect_serial!(a, b.ABCD)


function connect_parallel!(ABCD₁::ABCD_parameters, ABCD₂::ABCD_parameters)
    n = Int(size(ABCD₁, 1)/2)

    (a₁, b₁, c₁, d₁) = (ABCD₁[1:n,1:n], ABCD₁[1:n,n+1:end], ABCD₁[n+1:end,1:n], ABCD₁[n+1:end, n+1:end])
    (a₂, b₂, c₂, d₂) = (ABCD₂[1:n,1:n], ABCD₂[1:n,n+1:end], ABCD₂[n+1:end,1:n], ABCD₂[n+1:end, n+1:end])

    if (n == 1)
        a = (b₂[1] * a₁[1] + b₁[1] * a₂[1]) / (b₁[1] + b₂[1])
        b = b₁[1] * b₂[1] / (b₁[1] + b₂[1])
        c = c₁[1] + c₂[1] + (d₂[1] - d₁[1]) * (a₁[1] - a₂[1]) / (b₁[1] + b₂[1])
        d = d₁[1] + (d₂[1] - d₁[1]) * b₁[1] / (b₁[1] + b₂[1])
    else
        I = convert(Array{Basic}, Diagonal([1 for dummy in 1:n]))
        if all(b₁[i] == 0 for i in 1:n)
            a = a₂
            b = zeros(Basic, n, n)
            c = c₁ + c₂ + (d₂ - d₁) * (b₂ \ (a₁ - a₂))
            d = d₁
        elseif all(b₂[i] == 0 for i in 1:n)
            a = a₁
            b = zeros(Basic, n, n)
            c = c₁ + c₂ + (d₂ - d₁) * (b₁ \ (a₁ - a₂))
            d = d₁
        else
            a = inv(inv(b₁) + inv(b₂))  * (inv(b₁) * a₁ + inv(b₂) * a₂)
            b = inv(inv(b₁) + inv(b₂))
            c = c₁ + c₂ + (d₂ - d₁) * inv(b₁ + b₂) * (a₁ - a₂)
            d = d₁ + (d₂ - d₁) * inv(b₁ + b₂) * b₁
        end
    end

    ABCD = vcat(hcat(a,b), hcat(c,d))
    return ABCD
end
connect_parallel!(a::Element, b::Element) = connect_parallel!(a.ABCD, b.ABCD)
connect_parallel!(a::Element, b::ABCD_parameters) = connect_parallel!(a.ABCD, b)
connect_parallel!(a::ABCD_parameters, b::Element) = connect_parallel!(a, b.ABCD)


function closing_impedance(ABCD :: ABCD_parameters, Zₜ :: Union{Array{Basic}, Int, Float64}, direction = :output)
    n = Int(size(ABCD, 1)/2)
    (a, b, c, d) = (ABCD[1:n,1:n], ABCD[1:n,n+1:end], ABCD[n+1:end,1:n], ABCD[n+1:end, n+1:end])

    Zₑ = 0
    if (n == 1)
        if (direction == :output)
            Zₑ = (a .* Zₜ + b) ./ (c .* Zₜ + d)
        else
            Zₑ = (d .* Zₜ - b) ./ (c .* Zₜ - a)
        end
    else
        I = convert(Array{Basic}, Diagonal([1 for dummy in 1:n]))
        if (direction == :output)
            Zₑ = (a * Zₜ + b) * inv(c * Zₜ + d)
        else
            Zₑ = inv(Zₜ * c - a) * (Zₜ * d - b)
        end
    end
    return Zₑ
end

function closing_impedance(ABCD :: Array{Complex}, Zₜ :: Union{Array{Complex}, Int, Float64, Complex}, direction = :output)
    n = Int(size(ABCD, 1)/2)
    (a, b, c, d) = (ABCD[1:n,1:n], ABCD[1:n,n+1:end], ABCD[n+1:end,1:n], ABCD[n+1:end, n+1:end])

    Zₑ = 0
    if (n == 1)
        if (direction == :output)
            Zₑ = (a .* Zₜ + b) ./ (c .* Zₜ + d)
        else
            Zₑ = (d .* Zₜ - b) ./ (c .* Zₜ - a)
        end
    else
        I = convert(Array{Complex}, Diagonal([1 for dummy in 1:n]))
        if (direction == :output)
            Zₑ = (a * Zₜ + b) * pinv(c * Zₜ + d)
        else
            Zₑ = pinv(Zₜ * c - a) * (Zₜ * d - b)
        end
    end
    return Zₑ
end


closing_impedance(elem :: Element, Zₜ :: Union{Array{Basic}, Array{Complex}, Int, Float64, Complex}, direction = :output) =
    closing_impedance(elem.ABCD, Zₜ, direction)

function get_nodes(element::Element)
    return values(element.pins)
end

# gets
function get_nodes(element::Element, pin::Symbol)
    array = Symbol[]
    for (key, val) in element.pins
        !occursin(string(pin)[1:2], string(key)) && push!(array, val)
    end
    return array
end

"""
    function save_data(element :: Element, file_name :: String; omega_range = (-3, 5, 1000),
        scale = :log)
Used for saving data of the specific component defined as `Element`. The function calls component
specific function. The data is saved in csv textual file.

Additional parameter `omega_range` provides possibility to manually add the
frequency scale for saving data. Scale can be given as logarithmic (`:log`) and
linear (`:lin`).
"""
function save_data(element :: Element, file_name :: String; omega_range = (-3, 5, 1000),
    scale = :log)
    (min_ω, max_ω, n_ω) = omega_range
    n = (max_ω - min_ω) / n_ω
    if (scale == :log)
        omegas = [exp10(min_ω)*10^(i*n) for i in 1:Int(n_ω)]
    else
        omegas = [min_ω+i*n for i in 1:Int(n_ω)]
    end

    file_name = string("./files/",  file_name)
    save_data(element.element_value, file_name, omega_range = omega_range, scale = scale)
end

"""
    function plot_data(element :: Element)
Component defined function for plotting data. It differs from component to component.
"""
function plot_data(element :: Element; omega_range = (-3, 5, 1000),
    scale = :log)
    (min_ω, max_ω, n_ω) = omega_range
    n = (max_ω - min_ω) / n_ω
    if (scale == :log)
        omegas = [exp10(min_ω)*10^(i*n) for i in 1:Int(n_ω)]
    else
        omegas = [min_ω+i*n for i in 1:Int(n_ω)]
    end

    plot_data(element.element_value, omegas)
end

include("../Network.jl")
include("element_types.jl")
