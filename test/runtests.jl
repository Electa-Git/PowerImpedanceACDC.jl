using JLD2
using PowerImpedanceACDC
using Test


@testset "PowerImpedanceACDC" begin
        
        include("power_flow_test.jl")
        include("imp_test.jl")
        include("adm_MMC_test.jl")
        include("solvers_test.jl")
        include("adm_OHL_test.jl")

end

