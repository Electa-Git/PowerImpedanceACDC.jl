abstract type Transmission_line end

function eval_y(tl :: Transmission_line, s :: Complex)
    return abcd_to_y(eval_abcd(tl, s))
end

function make_power_flow_ac!(tl :: Transmission_line , dict :: Dict{String, Any},
            global_dict :: Dict{String, Any})
    key = string(length(dict["branch"]))
    ((dict["branch"])[string(key)])["transformer"] = false
    ((dict["branch"])[string(key)])["tap"] = 1
    ((dict["branch"])[string(key)])["shift"] = 0
    ((dict["branch"])[string(key)])["c_rating_a"] = 1

    abcd = eval_abcd(tl, global_dict["omega"] * 1im)
    n = Int(size(abcd, 1)/2)
    Z = (abcd[1:n,n+1:end])[1,1] / global_dict["Z"]
    Y = (abcd[n+1:end,1:n] * inv(abcd[n+1:end,n+1:end]))[1,1] * global_dict["Z"]

    ((dict["branch"])[string(key)])["br_r"] = real(Z)
    ((dict["branch"])[string(key)])["br_x"] = imag(Z)
    ((dict["branch"])[string(key)])["g_fr"] = real(Y)
    ((dict["branch"])[string(key)])["b_fr"] = imag(Y)
    ((dict["branch"])[string(key)])["g_to"] = real(Y)
    ((dict["branch"])[string(key)])["b_to"] = imag(Y)

end

function make_power_flow_dc!(tl :: Transmission_line, dict :: Dict{String, Any},
            global_dict :: Dict{String, Any})
    key = length(dict["branchdc"])
    ((dict["branchdc"])[string(key)])["l"] = 0
    ((dict["branchdc"])[string(key)])["c"] = 0

    abcd = eval_abcd(tl, 1e-6*1im)
    n = Int(size(abcd, 1)/2)
    Z = (abcd[1:n,n+1:end])[1,1] / global_dict["Z"]
    Y = (abcd[n+1:end,1:n] * inv(abcd[n+1:end,n+1:end]))[1,1] * global_dict["Z"]
    ((dict["branchdc"])[string(key)])["r"] = real(Z)
end

function save_data(tl :: Transmission_line, file_name :: String, omegas)
    println(typeof(tl))
    open(string(file_name, "_z.txt"), "w") do f
        for omega in omegas
            (Z,Y) = eval_parameters(tl, 1im*omega)
            writedlm(f, [omega reshape(Z, 1, length(Z))], ",")
        end
    end

    open(string(file_name, "_y.txt"), "w") do f
        for omega in omegas #omegas= 2*pi* 10 .^range(min_ω, max_ω, length= n_ω)
            (Z,Y) = eval_parameters(tl, 1im*omega) #receives in input the transmission line mutable struct (could be also cable or OHL) and S and computes the Y and Z pul matrices
            writedlm(f, [omega reshape(Y, 1, length(Y))], ",") #reshape changes the order of the elements e.g. assuming Y [3x3]= [1 2 3; 4 5 6; 7 8 9]-> reshape(Y,1,length(Y)) -> Y= [1 4 7; 2 5 8; 3 6 9] 
        end
    end

    open(string(file_name, "_abcd.txt"), "w") do f
        for omega in omegas
            abcd = eval_abcd(tl, 1im*omega)
            writedlm(f, [omega reshape(abcd, 1, length(abcd))], ",")
        end
    end
end

function plot_data(tl :: Transmission_line, omegas)
    Z_series = []
    Y_shunt = []
    for omega in omegas
        (Z,Y) = eval_parameters(tl, 1im*omega)
        push!(Z_series, Z)
        push!(Y_shunt, Y)
    end
    n = Int(sqrt(length(Z_series[1, :])))

    p = bode(Z_series, omega = omegas, titles = reshape([string("Z_{", i, "," , j, "}")
                    for j in 1:n for i in 1:n],n,n))
    save_plot(p, "files/z_series_bode")
    p = bode(Y_shunt, omega = omegas, titles = reshape([string("Y_{", i, "," , j, "}")
                    for j in 1:n for i in 1:n],n,n))
    save_plot(p, "files/y_shunt_bode")
end
