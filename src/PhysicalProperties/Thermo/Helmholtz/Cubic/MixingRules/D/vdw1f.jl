
struct vdw1f_D{TVector, TMatrix, Alpha} <: Cubic_D
    a_c::TVector
    kij::TMatrix
    alpha::Alpha
end

function CubicD(T, n, Model::vdw1f_D)
    sqrt_a = sqrt.(Model.a_c .* Alpha(T, Model.alpha))

    return n' * (((sqrt_a * sqrt_a') .* (1.0 .- Model.kij)) * n)
end