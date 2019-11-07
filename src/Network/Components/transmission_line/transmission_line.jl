export transmission_line
export Transmission_line, Conductors, Groundwires


@with_kw mutable struct Conductors
    nᵇ :: Int = 1                       # number of bundles (phases)
    nˢᵇ :: Int = 1                      # number of subconductors per bundle
    yᵇᶜ :: Union{Int, Float64}  = 0     # height above the ground of the lowest bundle  [m]
    Δyᵇᶜ :: Union{Int, Float64} = 0     # vertical offset between the bundles   [m]
    Δxᵇᶜ :: Union{Int, Float64} = 0     # horizontal offset between the lowest bundles  [m]
    Δ̃xᵇᶜ :: Union{Int, Float64} = 0     # horizontal offset in group of bundles    [m]
    dˢᵃᵍ :: Union{Int, Float64} = 0     # sag offset    [m]
    dˢᵇ :: Union{Int, Float64} = 0      # subconductor spacing (symmetric)  [m]
    rᶜ :: Union{Int, Float64} = 0       # conductor radius  [m]
    Rᵈᶜ :: Union{Int, Float64} = 0      # DC resistance for the entire conductor [Ω/m]
    gᶜ :: Union{Int, Float64} = 1e-11   # shunt conductance
    μᵣᶜ :: Union{Int, Float64} = 1      # relative conductor permeability
    positions :: Tuple{Vector{Union{Int, Float64}}, Vector{Union{Int, Float64}}} = ([],[])   # add absolute positions manually
    organization :: Symbol = Symbol()
end

@with_kw mutable struct Groundwires
    nᵍ :: Int = 0                        # number of groundwires (typically 0 or 2)
    Δxᵍ :: Union{Int, Float64} = 0       # horizontal ofsset between groundwires [m]
    Δyᵍ :: Union{Int, Float64} = 0       # vertical offset between the lowest conductor and groundwires  [m]
    rᵍ :: Union{Int, Float64} = 0        # ground wire radius  [m]
    dᵍˢᵃᵍ ::  Union{Int, Float64} = 0    # sag offset [m]
    Rᵍᵈᶜ :: Union{Int, Float64} = 0      # groundwire DC resistance [Ω/m]
    μᵣᵍ :: Union{Int, Float64} = 1       # relative groundwire permeability
    positions :: Tuple{Vector{Union{Int, Float64}}, Vector{Union{Int, Float64}}} = ([],[])    # add absolute positions manually
    eliminate :: Bool = true
end

@with_kw mutable struct Transmission_line
    length :: Union{Int, Float64} = 0       # line length [km]
    conductors :: Conductors = Conductors()
    groundwires :: Groundwires = Groundwires()
    earth_parameters :: NTuple{N, Union{Int,Float64}} where N = (1,1,1) # (μᵣ_earth, ϵᵣ_earth, ρ_earth) in units ([], [], [Ωm])

    # formed parameters
    x_array :: Vector{Union{Int, Float64}} = []
    y_array :: Vector{Union{Int, Float64}} = []
    r_array :: Vector{Union{Int, Float64}} = []
    ρ_array :: Vector{Union{Int, Float64}} = []
    μ_array :: Vector{Union{Int, Float64}} = []
end

"""
    transmission_line(;args...)
Generates the element `elem` with the  `element_value` of the type `Transmission_line`. Arguments should be given in the
form of struct `Transmission_line` fields:
...
- length - line length [m]
- conductors - defined in the
```julia
struct Conductors
    nᵇ :: Int = 1                       # number of bundles (phases)
    nˢᵇ :: Int = 1                      # number of subconductors per bundle
    yᵇᶜ :: Union{Int, Float64}  = 0     # height above the ground of the lowest bundle  [m]
    Δyᵇᶜ :: Union{Int, Float64} = 0     # vertical offset between the bundles   [m]
    Δxᵇᶜ :: Union{Int, Float64} = 0     # horizontal offset between the lowest bundles  [m]
    Δ̃xᵇᶜ :: Union{Int, Float64} = 0     # horizontal offset in group of bundles    [m]
    dˢᵃᵍ :: Union{Int, Float64} = 0     # sag offset    [m]
    dˢᵇ :: Union{Int, Float64} = 0      # subconductor spacing (symmetric)  [m]
    rᶜ :: Union{Int, Float64} = 0       # conductor radius  [m]
    Rᵈᶜ :: Union{Int, Float64} = 0      # DC resistance for the entire conductor [Ω/m]
    gᶜ :: Union{Int, Float64} = 1e-11   # shunt conductance
    μᵣᶜ :: Union{Int, Float64} = 1      # relative conductor permeability
    positions :: Tuple{Vector{Union{Int, Float64}}, Vector{Union{Int, Float64}}} = ([],[])   # add absolute positions manually
    organization :: Symbol = Symbol()
end
```
- groundwires - defined in the
```julia
struct Groundwires
    nᵍ :: Int = 0                        # number of groundwires (typically 0 or 2)
    Δxᵍ :: Union{Int, Float64} = 0       # horizontal ofsset between groundwires [m]
    Δyᵍ :: Union{Int, Float64} = 0       # vertical offset between the lowest conductor and groundwires  [m]
    rᵍ :: Union{Int, Float64} = 0        # ground wire radius  [m]
    dᵍˢᵃᵍ ::  Union{Int, Float64} = 0    # sag offset [m]
    Rᵍᵈᶜ :: Union{Int, Float64} = 0      # groundwire DC resistance [Ω/m]
    μᵣᵍ :: Union{Int, Float64} = 1       # relative groundwire permeability
    positions :: Tuple{Vector{Union{Int, Float64}}, Vector{Union{Int, Float64}}} = ([],[])    # add absolute positions manually
    eliminate :: Bool = true
end
```
- earth\\_parameters - with default value `(1,1,1)` and meaning (μᵣ\\_earth, ϵᵣ\\_earth, ρ\\_earth) in units ([], [], [Ωm])

Example:
```
transmission_line(length = 227e3, conductors = Conductors(nᵇ = 2, nˢᵇ = 2, organization = :flat,
Rᵈᶜ = 0.06266, rᶜ = 0.01436, yᵇᶜ = 27.5, Δxᵇᶜ = 11.8, dˢᵇ = 0.4572, dˢᵃᵍ = 10),
earth_parameters = (1,1,100),
groundwires = Groundwires(nᵍ = 2, Δxᵍ = 6.5, Δyᵍ = 7.5, Rᵍᵈᶜ = 0.9196, rᵍ = 0.0062, dᵍˢᵃᵍ = 10))
```
"""
function transmission_line(;args...)
    tl = Transmission_line()
    for (key, val) in kwargs_pairs(args)
        if in(key, propertynames(tl))
            setfield!(tl, key, val)
        end
    end

    # check definitions and calculate parameters
    function estimate_flat(c :: Conductors)
        if (c.nᵇ == 2)
            return ([-c.Δxᵇᶜ/2 c.Δxᵇᶜ/2], [c.yᵇᶜ c.yᵇᶜ])
        elseif (c.nᵇ == 3)
            return ([-c.Δxᵇᶜ 0 c.Δxᵇᶜ], [c.yᵇᶜ c.yᵇᶜ c.yᵇᶜ])
        elseif (c.nᵇ == 6)
            return ([-c.Δxᵇᶜ 0 c.Δxᵇᶜ -c.Δxᵇᶜ 0 c.Δxᵇᶜ], [c.yᵇᶜ c.yᵇᶜ c.yᵇᶜ c.yᵇᶜ+c.Δyᵇᶜ c.yᵇᶜ+c.Δyᵇᶜ c.yᵇᶜ+c.Δyᵇᶜ])
        else
            throw(ArgumentError("Invalid definition of flat conductor organization."))
        end
    end

    function estimate_vertical(c :: Conductors)
        return ([0 for i in 1:nᵇ], [c.yᵇᶜ+i*c.Δyᵇᶜ for i in 1:nᵇ])
    end

    function estimate_delta(c :: Conductors)
        if (c.nᵇ % 3 == 0)
            if (c.nᵇ == 3)
                return ([-c.Δxᵇᶜ/2 c.Δxᵇᶜ/2 0], [c.yᵇᶜ c.yᵇᶜ+c.Δyᵇᶜ c.yᵇᶜ])
            elseif (c.nᵇ == 6)
                return ([-c.Δxᵇᶜ/2-c.Δ̃xᵇᶜ -c.Δxᵇᶜ/2-c.Δ̃xᵇᶜ/2 -c.Δxᵇᶜ/2 c.Δxᵇᶜ/2 c.Δxᵇᶜ/2+c.Δ̃xᵇᶜ/2 c.Δxᵇᶜ/2+c.Δ̃xᵇᶜ],
                        [c.yᵇᶜ c.yᵇᶜ+c.Δyᵇᶜ c.yᵇᶜ c.yᵇᶜ c.yᵇᶜ+c.Δyᵇᶜ c.yᵇᶜ])
            end
        else
            throw(ArgumentError("Delta cannot be constructed from $(c.nᵇ) conductors."))
        end
    end

    function estimate_concentric(c :: Conductors)
        if (c.nᵇ % 3 == 0)
            if (c.nᵇ == 3)
                return ([-c.Δ̃xᵇᶜ 0 0], [c.yᵇᶜ+c.Δyᵇᶜ c.yᵇᶜ c.yᵇᶜ+2c.Δyᵇᶜ])
            elseif (c.nᵇ == 6)
                return ([-c.Δxᵇᶜ/2-c.Δ̃xᵇᶜ -c.Δxᵇᶜ/2 -c.Δxᵇᶜ/2 c.Δxᵇᶜ/2 c.Δxᵇᶜ/2 c.Δxᵇᶜ/2+c.Δ̃xᵇᶜ],
                        [c.yᵇᶜ+c.Δyᵇᶜ c.yᵇᶜ c.yᵇᶜ+2c.Δyᵇᶜ c.yᵇᶜ c.yᵇᶜ+2c.Δyᵇᶜ c.yᵇᶜ+c.Δyᵇᶜ])
            end
        else
            throw(ArgumentError("Delta cannot be constructed from $(c.nᵇ) conductors."))
        end
    end

    function estimate_offset(c :: Conductors)
        if (c.nᵇ % 3 == 0)
            if (c.nᵇ == 3)
                return ([-c.Δ̃xᵇᶜ 0 0], [c.yᵇᶜ+c.Δyᵇᶜ c.yᵇᶜ c.yᵇᶜ+2c.Δyᵇᶜ])
            elseif (c.nᵇ == 6)
                return ([-c.Δxᵇᶜ/2-c.Δ̃xᵇᶜ -c.Δxᵇᶜ/2 -c.Δxᵇᶜ/2 c.Δxᵇᶜ/2 c.Δxᵇᶜ/2 c.Δxᵇᶜ/2+c.Δ̃xᵇᶜ],
                        [c.yᵇᶜ+c.Δyᵇᶜ c.yᵇᶜ c.yᵇᶜ+2c.Δyᵇᶜ c.yᵇᶜ c.yᵇᶜ+2c.Δyᵇᶜ c.yᵇᶜ+c.Δyᵇᶜ])
            end
        else
            throw(ArgumentError("Delta cannot be constructed from $(c.nᵇ) conductors."))
        end
    end

    dict_organization = Dict(:flat => (c :: Conductors) -> estimate_flat(c),
                             :vertical => (c :: Conductors) -> estimate_vertical(c),
                             :delta => (c :: Conductors) -> estimate_delta(c),
                             :concentric => (c :: Conductors) -> estimate_concentric(c),
                             :offset => (c :: Conductors) -> estimate_offset(c),
                             Symbol() => (c :: Conductors) -> (0, c.yᵇᶜ))

    function bundle_position(nᵇ::Int, dᵇᶜ::Union{Int,Float64})
         if (nᵇ == 1)
             return (0,0)
         else
             ϕ = 2π/nᵇ
             r = dᵇᶜ/2/sin(ϕ/2)
             ϕₛ = π/2
             if iseven(nᵇ)
                 ϕₛ += ϕ/2
             end
             xˢᵇ = []
             yˢᵇ = []
             for i in 1:nᵇ
                 push!(xˢᵇ, r*cos(ϕₛ))
                 push!(yˢᵇ, r*sin(ϕₛ))
                 ϕₛ += ϕ
             end
             return (xˢᵇ,yˢᵇ)
        end
    end

    μ₀ = 4*π*1e-7
    if !(in(tl.conductors.organization, keys(dict_organization)) || isempty(tl.conductors.positions))
        throw(ArgumentError("Conductor positions are not defined."))
    else
        (x,y) = (0,0)
        if in(tl.conductors.organization, keys(dict_organization))
            (x,y) = (dict_organization[tl.conductors.organization])(tl.conductors)
        else
            (x,y) = tl.conductors.positions
        end
        (xˢᵇ, yˢᵇ) = bundle_position(tl.conductors.nˢᵇ, tl.conductors.dˢᵇ)
        for i in 1:tl.conductors.nᵇ
            for j in 1:tl.conductors.nˢᵇ
                push!(tl.x_array, x[i] + xˢᵇ[j])
                push!(tl.y_array, y[i] + yˢᵇ[j] -2/3*tl.conductors.dˢᵃᵍ)
                push!(tl.r_array, tl.conductors.rᶜ)
                push!(tl.ρ_array, tl.conductors.Rᵈᶜ*1e-3)
                push!(tl.μ_array, tl.conductors.μᵣᶜ*μ₀)
            end
        end
    end

    # Ground wires
    n = tl.groundwires.nᵍ
    (x,y) = ([tl.groundwires.Δxᵍ*i for i in -(n-1)/2:1:(n-1)/2], [tl.conductors.yᵇᶜ+tl.groundwires.Δyᵍ for i in 1:n])
    for i in 1:n
        push!(tl.x_array, x[i])
        push!(tl.y_array, y[i]-2/3*tl.groundwires.dᵍˢᵃᵍ)
        push!(tl.r_array, tl.groundwires.rᵍ)
        push!(tl.ρ_array, tl.groundwires.Rᵍᵈᶜ*1e-3)
        push!(tl.μ_array, tl.groundwires.μᵣᵍ*μ₀)
    end

    elem = Element(input_pins = tl.conductors.nᵇ, output_pins = tl.conductors.nᵇ, element_value = tl)
end

function create_abcd!(element :: Element, tl :: Union{Transmission_line, Cable})
    n = Int(np(element))
    symbol = "tl"
    abcd_tf =  reshape([Basic(string(symbol,i)) for i in 1:n^2], n, n)
    return abcd_tf
end


function eval_parameters(tl :: Transmission_line, s :: Complex)
    (μᵣ_earth, ϵᵣ_earth, ρ_earth) = tl.earth_parameters
    μ₀ = 4*π*1e-7
    ϵ = 8.85e-12
    μ_earth = μᵣ_earth*μ₀
    ϵ_earth = ϵᵣ_earth*ϵ
    σ_earth = 1/ρ_earth

    Num = tl.conductors.nᵇ * tl.conductors.nˢᵇ + tl.groundwires.nᵍ
    # matrix initialization
    P = zeros(Complex, Num,Num)
    Z = zeros(Complex, Num,Num)

    dₑ = sqrt(1/(s * μ_earth * (σ_earth + s*ϵ_earth)))       # depth of penetration
    for iPhase in 1:Num
        (xᵢ, yᵢ, rᵢ) = (tl.x_array[iPhase], tl.y_array[iPhase], tl.r_array[iPhase])
        (ρ, μ) = (tl.ρ_array[iPhase]*(π*rᵢ^2), tl.μ_array[iPhase])
        m = sqrt(s*μ/ρ)

        for jPhase in 1:Num
            (xⱼ, yⱼ) = (tl.x_array[jPhase], tl.y_array[jPhase])
            (Dᵢⱼ, D̂ᵢⱼ, dᵢⱼ) = (sqrt((xᵢ-xⱼ)^2+(yᵢ+yⱼ)^2), sqrt((xᵢ-xⱼ)^2+(yᵢ+yⱼ+2*dₑ)^2), sqrt((xᵢ-xⱼ)^2+(yᵢ-yⱼ)^2))
            # jPhase==1 ? println(Dᵢⱼ, D̂ᵢⱼ, dᵢⱼ) : nothing

            if (iPhase != jPhase)
                Z[iPhase,jPhase] += s*log(D̂ᵢⱼ/dᵢⱼ)*μ₀/(2*π)
                P[iPhase,jPhase] = 1/(2*π*ϵ)*log(Dᵢⱼ/dᵢⱼ)
            else
                Z[iPhase,jPhase] += s*μ₀/(2*π)*log(D̂ᵢⱼ/rᵢ) + m*ρ/(2*π*rᵢ)*coth(0.733*m*rᵢ) + .3179*ρ/(π*rᵢ^2)
                P[iPhase,jPhase] = 1/(2*π*ϵ)*log(Dᵢⱼ/rᵢ)
            end
        end
    end

    #apply Kron elimination to bundled subconductors
    if (tl.groundwires.nᵍ + tl.conductors.nˢᵇ != 0)
        if (tl.conductors.nˢᵇ != 0)
            cond_noElim = [(i-1)*tl.conductors.nˢᵇ + 1 for i in 1:tl.conductors.nᵇ]
            cond_Elim = setdiff(1:Num, cond_noElim)
            n_noElim = length(cond_noElim)
            for iPhase in 1:tl.conductors.nᵇ
                cond_noElim_curr = tl.conductors.nˢᵇ*(iPhase - 1) + 1
                for iCond in tl.conductors.nˢᵇ*(iPhase-1)+2 : tl.conductors.nˢᵇ*iPhase
                    # println(iCond)
                    Z[:,iCond] -= Z[:,cond_noElim_curr]
                    Z[iCond,:] -= Z[cond_noElim_curr,:]

                    P[:,iCond] -= P[:,cond_noElim_curr]
                    P[iCond,:] -= P[cond_noElim_curr,:]
                end
            end
            Z = Z[[cond_noElim; cond_Elim],:]
            Z = Z[:,[cond_noElim; cond_Elim]]
            P = P[[cond_noElim; cond_Elim],:]
            P = P[:,[cond_noElim; cond_Elim]]
        end
        # matrix to be reduced
        Z = Z[1:n_noElim,1:n_noElim] - Z[1:n_noElim,1+n_noElim:end]*
        inv(Z[1+n_noElim:end,1+n_noElim:end])*Z[1+n_noElim:end,1:n_noElim]
        P = P[1:n_noElim,1:n_noElim] - P[1:n_noElim,1+n_noElim:end]*
        inv(P[1+n_noElim:end,1+n_noElim:end])*P[1+n_noElim:end,1:n_noElim]
    end

    #invert to get the Y matrix
    Y = s*pinv(P) + Diagonal([tl.conductors.gᶜ for dummy in 1:size(P,1)])

    return (Z, Y)
end


"""
    eval_abcd(tl::Transmission_line, s :: Complex)
Form ABCD representation from known values for Y and Z and write values in the dictionary
    Γ = √ZY
    Yᶜ = Z⁻¹γ
    ABCD = [cosh(Γl)    Yᶜ⁻¹sinh(Γl)
            Yᶜsinh(Γl)  cosh(Γl)]
"""
function eval_abcd(tl :: Transmission_line, s :: Complex)
    (Z, Y) = eval_parameters(tl, s)
    γ = sqrt(Z*Y)
    Yc = inv(Z) * γ

    n = tl.conductors.nᵇ
    abcd = zeros(Complex, 2n, 2n)

    abcd[1:n, 1:n] = cosh(γ*tl.length)
    abcd[1:n,n+1:end] = inv(Yc) * sinh(γ*tl.length)
    abcd[n+1:end,1:n] = Yc * sinh(γ*tl.length)
    abcd[n+1:end, n+1:end] = cosh(γ*tl.length)

    return abcd
end


function save_data(tl :: Union{Transmission_line, Cable}, file_name :: String, omegas)
    open(string(file_name, "_z.txt"), "w") do f
        for omega in omegas
            (Z,Y) = eval_parameters(tl, 1im*omega)
            writedlm(f, [omega reshape(Z, 1, length(Z))], ",")
        end
    end

    open(string(file_name, "_y.txt"), "w") do f
        for omega in omegas
            (Z,Y) = eval_parameters(tl, 1im*omega)
            writedlm(f, [omega reshape(Y, 1, length(Y))], ",")
        end
    end

    open(string(file_name, "_abcd.txt"), "w") do f
        for omega in omegas
            abcd = eval_abcd(tl, 1im*omega)
            writedlm(f, [omega reshape(abcd, 1, length(abcd))], ",")
        end
    end
end

function plot_data(tl :: Union{Transmission_line, Cable}, omegas)
    Z_series = []
    Y_shunt = []
    for omega in omegas
        (Z,Y) = eval_parameters(tl, 1im*omega)
        push!(Z_series, Z)
        push!(Y_shunt, Y)
    end
    n = tl.conductors.nᵇ

    bode(Z_series, omega = omegas, titles = reshape([string("Z_{", i, "," , j, "}")
                    for j in 1:n for i in 1:n],n,n))
    save_plot("files/z_series")
    bode(Y_shunt, omega = omegas, titles = reshape([string("Y_{", i, "," , j, "}")
                    for j in 1:n for i in 1:n],n,n))
    save_plot("files/dc_side")
end
