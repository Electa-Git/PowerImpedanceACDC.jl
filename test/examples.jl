
@testset "impedance" begin
    include("tests/impedance_test.jl")
    @test isapprox(round(abs(imp[end][])), round(abs(1.5*omega[end]*1im + 2)))
end

@testset "transformer" begin
    include("tests/transformer_test.jl")
    @test isapprox(net.elements[:t].element_value.Rₚ, 0.7399359605650204)
end
