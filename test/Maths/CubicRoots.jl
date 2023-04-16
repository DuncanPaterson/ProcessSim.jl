@testset "CubicRoots" begin
    @test all(ProcessSim.CubicRoots(-3.0, 3.0, -1.0) .≈ [1.0, 1.0, 1.0])
    @test all(ProcessSim.CubicRoots(-9000.001, -9.999991e6, 10000.0) .≈ [-1.0e3, 0.001, 1.0e4])
end