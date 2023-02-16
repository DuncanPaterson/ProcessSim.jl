
struct MockAlpha

end

function ProcessSim.Alpha(T::Float64, Model::MockAlpha)
    return [1.0, 1.0, 1.0]
end


@testset "vdw1f D" begin
    a_c = [0.5, 0.4, 0.3]
    params = ProcessSim.vdw1f_D(a_c, [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0], MockAlpha())
    n = [0.3, 0.3, 0.4]
    T = 120.0
    f =  ProcessSim.Alpha(T, MockAlpha())
    @test ProcessSim.CubicD(T, n, params) ≈ n' * ((sqrt.(a_c * a_c')) * n)
end