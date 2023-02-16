
struct MockD end

struct MockB end

function ProcessSim.CubicD(n, T, Model::MockD)
    return 1.0
end

function ProcessSim.CubicB(n, Model::MockB)
    return 1.0
end


@testset "General Cubic" begin
    params = GeneralCubic(1.0, 0.0, MockD(), MockB())

    T = 100.0
    V = 2.0
    @test ProcessSim.Helmholtz(T, V, n, params) ≈ n' * ((sqrt.(a_c * a_c')) * n)
end