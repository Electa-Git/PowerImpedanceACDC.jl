export nyquistplot

function nyquistplot(L, omega; zoom :: String = "", SM :: String = "", title :: String = "Nyquist plot")

    # Determine eigenvalues of the loop gain matrix at each frequency point 
    Λ = eigvals.(L)

    # Number of frequency points
    omegaₙ = size(Λ, 1)
    # Number of eigenvalues for each frequency point
    λₙ = size(Λ[1], 1)
    # L has dimensions λₙxλₙxomegaₙ 

    # Sorting algorithm of eigenvalues to create continuous eigenloci in Nyquist plot 
    # Create empty matrix for sorted eigenvalues 
    Λₛ = zeros(Complex{Float64}, omegaₙ, λₙ)
    # First row of eigenvalues remains the same 
    Λₛ[1,:] = Λ[1]

    for i in 2:omegaₙ
        
        # Row of eigenvalues at previous frequency point 
        Λ₁ = Λ[i-1]
        # Row of eigenvalues at current frequency point
        Λ₂ = Λ[i]

        # Algorithm: 
        # see Munkres.jl https://juliapackages.com/p/munkres
        # see eigenshuffle.m https://nl.mathworks.com/matlabcentral/fileexchange/22885-eigenshuffle
        # Sorting eigenvalues of row i to obtain minimum distance between row i and row i-1

        Λ₁_Re = real(Λ₁)
        Λ₂_Re = real(Λ₂)

        Λ₁_Re_Grid = [i for i in Λ₁_Re, j in 1:length(Λ₂_Re)]
        Λ₂_Re_Grid = [j for i in 1:length(Λ₁_Re), j in Λ₂_Re]

        Λ_Re = abs.(Λ₁_Re_Grid - Λ₂_Re_Grid)

        Λ₁_Im = imag(Λ₁)
        Λ₂_Im = imag(Λ₂)

        Λ₁_Im_Grid = [i for i in Λ₁_Im, j in 1:length(Λ₂_Im)]
        Λ₂_Im_Grid = [j for i in 1:length(Λ₁_Im), j in Λ₂_Im]

        Λ_Im = abs.(Λ₁_Im_Grid - Λ₂_Im_Grid)

        # Possible distances between eigenvalues
        D = sqrt.(Λ_Re.^2 + Λ_Im.^2)
        # Sequence of sorted eigenvalues
        S = munkres(D)

        Λ[i] = Λ[i][S]
        # Insert sorted eigenvalues in matrix Λₛ
        Λₛ[i, :] = Λ[i]

    end

    Λ = Λₛ

    # Plotting of eigenloci

    x = real(Λ)
    y = imag(Λ)
    
    # plotly() # To activate interactive plot
    # Automatic setting of plot limits 
    x_min = minimum(x)
    x_max = maximum(x)
    y_abs = abs.(y)
    y_max = maximum(y_abs)
    # Zoom in on (-1, 0j) point 
    if zoom == "yes" 
        xlims = (-2.6, 0.6)
        ylims = (-1.6, 1.6)
    else
        xlims = (x_min - 0.1*abs(x_min), x_max + 0.1*abs(x_max))
        ylims = (-1.1*y_max, 1.1*y_max)
    end
    # Indicating unity circle in red 
    θ = 0:0.01:2*π
    np = plot(cos.(θ), sin.(θ), c = :red, linewidth = 3, line = :dash, label = :none)
    plot!(xlabel = "Real axis", ylabel = "Imaginary axis", title = title, xlims = xlims, ylims = ylims, framestyle = :box, legend=:topleft) 
    vline!([0], c = :black, line = :dash, label = :none)
    hline!([0], c = :black, line = :dash, label = :none)
    # Indicating (-1, 0j) point with red cross 
    scatter!([-1], [0]; markercolor = :red, markershape = :xcross, markersize = 10, label = :none)
    # Loci from fmin to fmax are indicated with solid line, Loci from -fmin to -fmax are indicated with dashed line
    for i in 1:size(Λ,2) 
            λ = Λ[:,i]
            x = real(λ)
            y = imag(λ)
            y_abs = abs.(y)
            y_max = maximum(y_abs)
            color = palette(:default) #palette(:tab10)
            plot!(x, y, linewidth = 3, c = color[i], label = "Lambda " * string(i))
            plot!(x, -y, linewidth = 3, c = color[i], linestyle = :dash, label = :none)
            index = findall(x -> x == y_max, y_abs)
            index = index[1]
            # Arrows added to show direction of each eigenloci
                    if index != length(λ)
                            if x[index+1] < x[index]
                                    scatter!([x[index]], [y[index]]; markercolor = color[i], markershape = :ltriangle, markersize = 10, label = :none)
                                    scatter!([x[index]], [-y[index]]; markercolor = color[i], markershape = :rtriangle, markersize = 10, label = :none)
                            else
                                    scatter!([x[index]], [y[index]]; markercolor = color[i], markershape = :rtriangle, markersize = 10, label = :none)
                                    scatter!([x[index]], [-y[index]]; markercolor = color[i], markershape = :ltriangle, markersize = 10, label = :none)
                            end
                    end 
    end
    plot!(size=(1000,1000))
    display(np)

    x = real(Λ)
    y = imag(Λ)

    # Counting of clockwise and counterclockwise encirclements 
    cw = []
    ccw = []

    for i in 1:λₙ
        cwi = 0      #Number of clockwise encirclements
        ccwi = 0     #Number of counterclockwise encirclements  
        for j in 2:omegaₙ
            if (y[j-1, i] < 0 && y[j, i] > 0 && x[j, i] < -1)
                cwi += 1
            elseif (y[j-1, i] > 0 && y[j, i] < 0  && x[j-1, i] < -1)
                ccwi += 1
            end
        end
        push!(cw, cwi)
        push!(ccw, ccwi)
    end

    # Nyquist stability criterion
    # N: Net number of clockwise encirclements
    # P: Number of RHP poles of loop gain (matrix)
    # Z: number of RHP poles of closed-loop system
    # N = Z - P 
    # If Z = N + P > 0, system is unstable 
    # P is equal to or greater than zero => System is unstable if N > 0

    N = sum(cw) - sum(ccw) 

    println("")
    if N > 0
        println("Result stability assessement: Unstable system \n")
    elseif N < 0
        println("Result stability assessement: Unstable subsystem \n")
    else
        println("Result stability assessement: Stable system if subsystems are stable \n")
    end

    ##### ----- STABILITY MARGINS ----- #####

    # Use stabilitymargin.jl to calculate stability margins of eigenloci
    # Possibilities for SM are "PM", "GM" and "VM"
    if SM != "no"
        for i in 1:λₙ
            println("Stability margins for λ", i, ":")
            stabilitymargin(Λ[:,i], omega, SM = SM)
        end
    end

end


