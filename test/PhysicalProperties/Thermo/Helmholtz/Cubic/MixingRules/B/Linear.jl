

@testset "Linear B" begin
    b_c = [0.5, 0.4, 0.3]
    params = ProcessSim.Linear_B(b_c)
    n = [0.3, 0.3, 0.4]
    @test ProcessSim.CubicB(n, params) ≈ n' * b_c
end