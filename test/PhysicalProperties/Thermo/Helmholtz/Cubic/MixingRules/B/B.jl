

@testset "B from crit" begin
    @test ProcessSim.BFromCrit(0.8664, 120.0, 12.0e5) / 8.3144626 ≈ 0.8664 * 120.0 / 12.0e5
end

include("Linear.jl")
