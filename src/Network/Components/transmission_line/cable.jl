export cable
export Cable, Insulator, Conductor

@with_kw mutable struct Conductor
    rᵢ :: Union{Int, Float64} = 0              # inner radius
    rₒ :: Union{Int, Float64} = 0              # outer radius
    ρ  :: Union{Int, Float64} = 0              # conductor resistivity [Ωm]
    μᵣ  :: Union{Int, Float64} = 1             # relative permeability
end

@with_kw mutable struct Insulator
    rᵢ :: Union{Int, Float64} = 0               # inner radius
    rₒ :: Union{Int, Float64} = 0               # outer radius
    ϵᵣ :: Union{Int, Float64} = 1               # relative permittivity
    μᵣ :: Union{Int, Float64} = 1               # relative permeability
end

@with_kw mutable struct Cable
    length :: Union{Int, Float64} = 0                    # line length [km]
    conductors :: OrderedDict{Symbol, Conductor} = OrderedDict{Symbol, Conductor}()
    insulators :: OrderedDict{Symbol, Insulator} = OrderedDict{Symbol, Insulator}()
    positions :: Vector{Tuple{Real,Real}} = []
    earth_parameters :: NTuple{N, Union{Int,Float64}} where N = (1,1,1) # (μᵣ_earth, ϵᵣ_earth, ρ_earth) in units ([], [], [Ωm])
    configuration :: Symbol = :coaxial

    # formed parameters
    r_array :: Vector{Union{Int, Float64}} = []
    ρ_array :: Vector{Union{Int, Float64}} = []
    μ_array :: Vector{Union{Int, Float64}} = []
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
given with the struct `Conductor`
```julia
struct Conductor
    rᵢ :: Union{Int, Float64} = 0              # inner radius
    rₒ :: Union{Int, Float64} = 0              # outer radius
    ρ  :: Union{Int, Float64} = 0              # conductor resistivity [Ωm]
    μᵣ  :: Union{Int, Float64} = 1             # relative premeability
end
```
- insulators - dictionary with the key being symbol: I1, I2, I3 and I4, and the value
given with the struct `Insulator`
```julia
struct Insulator
    rᵢ :: Union{Int, Float64} = 0               # inner radius
    rₒ :: Union{Int, Float64} = 0               # outer radius
    ϵᵣ :: Union{Int, Float64} = 1               # relative permittivity
    μᵣ :: Union{Int, Float64} = 1               # relative permeability
end
```
- positions - given as an array in (x,y) format
- configuration - symbol with two possible values: coaxial (default) and pipe-type

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
    n = length(c.positions)
    elem = Element(input_pins = n, output_pins = n, element_value = c)
end

#     cond_dict = Dict("Cu" => (1.72e-8, 1),
#                      "Pb" => (22e-8, 1),
#                      "Al" => (2.83e-8, 1),
#                      "steel" => (18e-8, 10))
#     insulator_dict = Dict("XLPE" => (2.3, 1),
#                           "mass-impreganted" => (4.2, 1),
#                           "fluid-filled" => (3.5, 1))


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

    # add conductor series impedance
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
            # Zᵍ = s*μᵍ/(2π) * (-bessely0(m*dᵢⱼ) + 2/(4 + (m*x)^2) * exp(-m*H))
            Z[i,i] += Zᵍ
        end
        i += 1
    end

    # add insulators
    i = 1
    for key in keys(c.insulators)
        (rᵢ, rₒ, μ, ϵ) = (c.insulators[key].rᵢ, c.insulators[key].rₒ, c.insulators[key].μᵣ*μ₀, c.insulators[key].ϵᵣ*ϵ₀)
        Zⁱ = s*μ/(2π) * log(rₒ/rᵢ)
        Pⁱ = log(rₒ/rᵢ) / (2π*ϵ)

        Z[i,i] += Zⁱ
        for j in 1:i
            for k in 1:i
                P[j,k] += Pⁱ
            end
        end
        i += 1
    end

    # adding earth return impedance
    for i in 1:n
        Z[(i-1)*nₗ+1:i*nₗ, (i-1)*nₗ+1:i*nₗ] = Z[1:nₗ,1:nₗ]
        P[(i-1)*nₗ+1:i*nₗ, (i-1)*nₗ+1:i*nₗ] = P[1:nₗ,1:nₗ]
        for j in i+1:n
            m = sqrt(s*μᵍ/ρᵍ)
            H = c.positions[i][2] + c.positions[j][2]
            dᵢⱼ = sqrt((c.positions[i][1] - c.positions[j][1])^2 + (c.positions[i][2] - c.positions[j][2])^2)
            x = abs(c.positions[i][1] - c.positions[j][1])
            Zᵍ = s*μᵍ/(2π) * (-log(γ*m*dᵢⱼ/2) + 0.5 - 2*m*H/3)
            # Zᵍ = s*μᵍ/(2π) * (-bessely0(m*dᵢⱼ) + 2/(4 + (m*x)^2) * exp(-m*H))
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

    for i in 1:n
        for j in 1:n
            H = c.positions[i][2] + c.positions[j][2]
            x = abs(c.positions[i][1] - c.positions[j][1])
            y = abs(c.positions[i][2] - c.positions[j][2])
            (i == j) ? D₁ = max(maximum([r.rₒ for r in values(c.conductors)]), maximum([r.rₒ for r in values(c.insulators)])) : D₁ = sqrt(x^2 + y^2)
            (i == j) ? D₂ = H : D₂ = sqrt(x^2 + H^2)
            Pᵢⱼ = log(D₂ / D₁) / (2π*ϵ₀)

            for k in 1:nₗ
                for l in 1:nₗ
                    P[(i-1)nₗ+k, (j-1)nₗ+l] += Pᵢⱼ
                end
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
