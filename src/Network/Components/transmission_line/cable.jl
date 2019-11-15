export cable
export Cable, Insulator, Conductor

include("transmission_line.jl")

@with_kw mutable struct Conductor
    rᵢ :: Union{Int, Float64} = 0              # inner radius
    rₒ :: Union{Int, Float64} = 0              # outer radius
    ρ  :: Union{Int, Float64} = 0              # conductor resistivity [Ωm]
    μᵣ  :: Union{Int, Float64} = 1             # relative permeability

    A :: Union{Int, Float64} = 0               # nominal area
end

@with_kw mutable struct Insulator
    rᵢ :: Union{Int, Float64} = 0               # inner radius
    rₒ :: Union{Int, Float64} = 0               # outer radius
    ϵᵣ :: Union{Int, Float64} = 1               # relative permittivity
    μᵣ :: Union{Int, Float64} = 1               # relative permeability

    a :: Union{Int, Float64} = 0                # inner semiconductor outer radius
    b :: Union{Int, Float64} = 0                # outer semiconductor inner radius
end

@with_kw mutable struct Cable <: Transmission_line
    length :: Union{Int, Float64} = 0                    # line length [km]
    conductors :: OrderedDict{Symbol, Conductor} = OrderedDict{Symbol, Conductor}()
    insulators :: OrderedDict{Symbol, Insulator} = OrderedDict{Symbol, Insulator}()
    positions :: Vector{Tuple{Real,Real}} = []
    earth_parameters :: NTuple{N, Union{Int,Float64}} where N = (1,1,1) # (μᵣ, ϵᵣ, ρ) in units ([], [], [Ωm])
    configuration :: Symbol = :coaxial
    type :: Symbol = :underground   # or aerial
end

"""
    cable(;args...)
Generates the element `elem` with the  `element_value` of the type `Cable`. Arguments should be given in the
form of struct `Cable` fields:
- length - length of the cable in [m]
- earth\\_parameters - (μᵣᵍ, ϵᵣᵍ, ρᵍ) in units ([], [], [Ωm])
meaning ground (earth) relative premeability, relative permittivity and
ground resistivity
- conductors - dictionary with the key symbol being: C1, C2, C3 or C4, and the value
given with the struct `Conductor`. If the sheath consists of metalic screen and sheath,
then add screen with a key symbol SC and sheath with C2.
```julia
struct Conductor
    rᵢ :: Union{Int, Float64} = 0              # inner radius
    rₒ :: Union{Int, Float64} = 0              # outer radius
    ρ  :: Union{Int, Float64} = 0              # conductor resistivity [Ωm]
    μᵣ  :: Union{Int, Float64} = 1             # relative premeability

    A :: Union{Int, Float64} = 0               # nominal area

    screen_r :: Union{Int, Float64} = 0        # metalic screen outer radius
    screen_ρ :: Union{Int, Float64} = 0        # metalic screen resistivity [Ωm]
end
```
- insulators - dictionary with the key being symbol: I1, I2, I3 and I4, and the value
given with the struct `Insulator`. For the insulator 2 the semiconducting layers
can be added by specifying outer radius of the inner semiconducting layer and
inner radius of the outer semiconducting layer.
```julia
struct Insulator
    rᵢ :: Union{Int, Float64} = 0               # inner radius
    rₒ :: Union{Int, Float64} = 0               # outer radius
    ϵᵣ :: Union{Int, Float64} = 1               # relative permittivity
    μᵣ :: Union{Int, Float64} = 1               # relative permeability

    a :: Union{Int, Float64} = 0                # inner semiconductor outer radius
    b :: Union{Int, Float64} = 0                # outer semiconductor inner radius
end
```
- positions - given as an array in (x,y) format
- configuration - symbol with two possible values: coaxial (default) and pipe-type
- type - symbol representing underground or aerial cable
"""
function cable(;args...)
    c = Cable()

    for (key, val) in kwargs_pairs(args)
        if key == :positions
            for v in val
                push!(c.positions, v)
            end
        elseif in(key, propertynames(c))
            setfield!(c, key, val)
        elseif isa(val, Conductor)
            c.conductors[key] = val
        elseif isa(val, Insulator)
            c.insulators[key] = val
        end
    end

    # conversion procedure
    # core outer radius
    (c.conductors[:C1].A != 0) ? c.conductors[:C1].ρ = c.conductors[:C1].ρ * π *
                        c.conductors[:C1].rₒ^2 / c.conductors[:C1].A : nothing
    # add metalic screen conversions, equivalent sheat layer
    if in(:SC, keys(c.conductors))
        !in(:C2, keys(c.conductors)) && throw(ArgumentError("There must be present sheath together with screen layer."))
        if (c.conductors[:SC].A != 0)
            c.conductors[:C2].rᵢ = sqrt(c.conductors[:SC].rₒ^2 - c.conductors[:SC].A / π)
        else
            c.conductors[:C2].rᵢ = c.conductors[:SC].rᵢ
        end
        c.conductors[:C2].rₒ = sqrt((c.conductors[:C2].rₒ^2 - c.conductors[:SC].rₒ^2) *
                    c.conductors[:SC].ρ / c.conductors[:C2].ρ + c.conductors[:SC].rₒ^2)
        delete!(c.conductors, :SC)

        # change Insulator 1
        c.insulators[:I1].rₒ = c.conductors[:C2].rᵢ
        # change Insulator 2
        if in(:I2, keys(c.insulators))
            x = log(c.insulators[:I2].rₒ / c.conductors[:C2].rₒ) / log(c.insulators[:I2].rₒ / c.insulators[:I2].rᵢ)
            c.insulators[:I2].rᵢ = c.conductors[:C2].rₒ
            c.insulators[:I2].ϵᵣ *=  x
            c.insulators[:I2].μᵣ /= x
        end
    end
    # semiconductor configuration
    if in(:I1, keys(c.insulators)) && (c.insulators[:I1].a != 0)
        x = log(c.insulators[:I1].rₒ / c.insulators[:I1].rᵢ) / log(c.insulators[:I1].b / c.insulators[:I1].a)
        c.insulators[:I1].ϵᵣ *=  x
        N = 1.4
        c.insulators[:I1].μᵣ = c.insulators[:I1].μᵣ * (1 + 2π^2 * N^2 * (c.insulators[:I1].rₒ^2 - c.insulators[:I1].rᵢ^2) / log(c.insulators[:I1].rₒ / c.insulators[:I1].rᵢ))
    end

    n = length(c.positions)
    elem = Element(input_pins = n, output_pins = n, element_value = c)
end

function eval_parameters(c :: Cable, s :: Complex)
    (μᵣᵍ, ϵᵣᵍ, ρᵍ) = c.earth_parameters

    μ₀ = 4π*1e-7
    ϵ₀ = 8.85e-12
    μᵍ = μᵣᵍ*μ₀
    ϵᵍ = ϵᵣᵍ*ϵ₀
    σᵍ = 1/ρᵍ
    γ = 0.5772156649
    g = 1e-11

    nₗ = length(c.conductors)       # number of cable layers
    n = length(c.positions)        # number of cables
    Z = zeros(Complex,n*nₗ,n*nₗ)    # series impedance matrix
    P = zeros(Complex,n*nₗ,n*nₗ)

    # make series impedance
    i = 1
    for key in keys(c.conductors)
        (rᵢ, rₒ, μ, ρ) = (c.conductors[key].rᵢ, c.conductors[key].rₒ, c.conductors[key].μᵣ*μ₀, c.conductors[key].ρ)
        m = sqrt(s*μ/ρ)
        Δr = rₒ - rᵢ
        if (rᵢ != 0)
            Zᵃᵃ = ρ*m/(2π*rᵢ)*coth(m*Δr) - ρ/(2π*rᵢ*(rᵢ+rₒ))
            Zᵇᵇ = ρ*m/(2π*rₒ)*coth(m*Δr) + ρ/(2π*rₒ*(rᵢ+rₒ))
        else
            Zᵇᵇ = ρ*m/(2π*rₒ)*coth(0.733m*rₒ) + 0.3179ρ/(π*rₒ^2)
        end
        Zᵃᵇ = ρ*m/(π*(rₒ+rᵢ))*csch(m*Δr)

        Z[i,i] += Zᵇᵇ
        if (i > 1)
            Z[i,i-1] += -Zᵃᵇ
            Z[i-1,i] += -Zᵃᵇ
            Z[i-1,i-1] += Zᵃᵃ
        end
        if (i == nₗ)
            m = sqrt(s*μᵍ/ρᵍ)
            H = 2c.positions[1][2]
            dᵢⱼ = max(maximum([r.rₒ for r in values(c.conductors)]), maximum([r.rₒ for r in values(c.insulators)]))
            x = dᵢⱼ
            Zᵍ = s*μᵍ/(2π) * (-log(γ*m*dᵢⱼ/2) + 0.5 - 2*m*H/3)
            Z[i,i] += Zᵍ
        end
        i += 1
    end

    # make shunt admittance
    i = 1
    for key in keys(c.insulators)
        (rᵢ, rₒ, μ, ϵ) = (c.insulators[key].rᵢ, c.insulators[key].rₒ, c.insulators[key].μᵣ*μ₀, c.insulators[key].ϵᵣ*ϵ₀)
        Zⁱ = s*μ/(2π) * log(rₒ/rᵢ)
        Pⁱ = log(rₒ/rᵢ) / (2π*ϵ)

        Z[i,i] += Zⁱ
        P[1:i,1:i] += ones(i,i) * Pⁱ
        i += 1
    end


    for i in 1:n
        Z[(i-1)*nₗ+1:i*nₗ, (i-1)*nₗ+1:i*nₗ] = Z[1:nₗ,1:nₗ]
        P[(i-1)*nₗ+1:i*nₗ, (i-1)*nₗ+1:i*nₗ] = P[1:nₗ,1:nₗ]
        for j in i+1:n
            # adding earth return impedance
            m = sqrt(s*μᵍ/ρᵍ)
            H = c.positions[i][2] + c.positions[j][2]
            dᵢⱼ = sqrt((c.positions[i][1] - c.positions[j][1])^2 + (c.positions[i][2] - c.positions[j][2])^2)
            x = abs(c.positions[i][1] - c.positions[j][1])
            Zᵍ = s*μᵍ/(2π) * (-log(γ*m*dᵢⱼ/2) + 0.5 - 2*m*H/3)
            Z[i*nₗ, j*nₗ] += Zᵍ
            Z[j*nₗ, i*nₗ] += Zᵍ
        end
    end

    # reduction for represention core, sheath and armor
    for k in 1:n
        for l in 1:n
            for i in nₗ-1:-1:1
                for j in 1:i
                    Z[(l-1)nₗ+1:l*nₗ, (k-1)*nₗ + j] += Z[(l-1)nₗ+1:l*nₗ, (k-1)*nₗ + i+1]
                end
            end
            for i in nₗ-1:-1:1
                for j in 1:i
                    Z[(k-1)*nₗ + j, (l-1)nₗ+1:l*nₗ] += Z[(k-1)*nₗ + i+1, (l-1)nₗ+1:l*nₗ]
                end
            end
        end
    end

    if (c.type == :underground)
        for i in 1:n
            for j in 1:n
                H = c.positions[i][2] + c.positions[j][2]
                x = abs(c.positions[i][1] - c.positions[j][1])
                y = abs(c.positions[i][2] - c.positions[j][2])
                (i == j) ? D₁ = max(maximum([r.rₒ for r in values(c.conductors)]), maximum([r.rₒ for r in values(c.insulators)])) : D₁ = sqrt(x^2 + y^2)
                (i == j) ? D₂ = H : D₂ = sqrt(x^2 + H^2)
                Pᵢⱼ = log(D₂ / D₁) / (2π*ϵ₀)

                P[(i-1)nₗ+1:i*nₗ , (j-1)nₗ+1:j*nₗ] += ones(nₗ, nₗ) * Pᵢⱼ
            end
        end
    end

    #apply Kron elimination to additional layers
    cond_noElim = [(i-1)*nₗ + 1 for i in 1:n]
    cond_Elim = setdiff(1:n*nₗ, cond_noElim)
    n_noElim = length(cond_noElim)

    Z = Z[[cond_noElim; cond_Elim],:]
    Z = Z[:,[cond_noElim; cond_Elim]]
    P = P[[cond_noElim; cond_Elim],:]
    P = P[:,[cond_noElim; cond_Elim]]

    Z = Z[1:n_noElim,1:n_noElim] - Z[1:n_noElim,1+n_noElim:end]*inv(Z[1+n_noElim:end,1+n_noElim:end])*Z[1+n_noElim:end,1:n_noElim]
    P = P[1:n_noElim,1:n_noElim] - P[1:n_noElim,1+n_noElim:end]*inv(P[1+n_noElim:end,1+n_noElim:end])*P[1+n_noElim:end,1:n_noElim]

    Y = s*inv(P)

    return (Z,Y)
end

function eval_abcd(c :: Cable, s :: Complex)
    (Z, Y) = eval_parameters(c, s)
    γ = sqrt(Z*Y)
    Yc = inv(Z) * γ

    n = length(c.positions)
    abcd = zeros(Complex, 2n, 2n)

    abcd[1:n, 1:n] = cosh(γ*c.length)
    abcd[1:n,n+1:end] = inv(Yc) * sinh(γ*c.length)
    abcd[n+1:end,1:n] = Yc * sinh(γ*c.length)
    abcd[n+1:end, n+1:end] = cosh(γ*c.length)

    return abcd
end
