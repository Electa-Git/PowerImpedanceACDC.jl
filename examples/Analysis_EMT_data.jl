# Scrip that applies the Generalized Nyquist Criteria and evaluates passivity to data obtained from a frequency scan in an EMT software
using DelimitedFiles
using LinearAlgebra
using Plots

# Import the file
file = "./files/test_EMT.txt"
L, vars = readdlm(file,'\t',ComplexF64,'\n',header=true) # Read the file contents

freq = real(L[:,1]) # Frequency at the first column
L = L[:,2:end] # The matrix is flat at each frequency as [11 12 13 ... 21 22 23 ...]
N = round(Int64,sqrt(size(L,2))) # Dimension of the matrix at every frequency

Lreshape = [] # Reshape as a vector of many square matrices as frequencies
for i in 1:size(L,1)
    push!(Lreshape, transpose(reshape(L[i,:],(N,N)))) # Reshape transposes the values, so we need to undo it
end

# Apply Nyquist 
nyquist_external_data = nyquistplot(Lreshape, 2*pi.*freq, zoom = "yes", SM = "PM")
# savefig("nyquist_external_data.png")

# Compute passivity
passivity(Lreshape, 2*pi.*freq, "Loop-gain")

# Compute matrix condition number
cond_num_L = []
for i in 1:size(L,1)
    push!(cond_num_L, cond(Lreshape[i])) 
end
conditioning = plot( freq, cond_num_L, label = "Loop gain", c = :black, linewidth = 3, xaxis = :log10,yaxis = :log10)
plot!(xlabel = "Frequency [Hz]", ylabel = "Condition number", title = "Condition number", xlims = (freq[1],freq[end]),framestyle = :box, legend=:topright) 
display(conditioning)
