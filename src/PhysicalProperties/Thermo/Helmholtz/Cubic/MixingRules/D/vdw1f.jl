
struct vdw1f_D{TVector, TMatrix, Alpha} <: Cubic_D
    a_c::TVector
    kij::TMatrix
    alpha::Alpha
end

function CubicD(n, T, Model::Linear_B)

    sqrt_a = sqrt.(Model.a_c .* CubicAlpha(T, Model::Alpha))

    return n' * (((sqrt_a * sqrt_a') .* (1.0 - kij)) * n)
end