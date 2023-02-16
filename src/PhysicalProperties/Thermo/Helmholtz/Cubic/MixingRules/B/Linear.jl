
struct Linear_B{TVector} <: Cubic_B
    b_c::TVector
end

function CubicB(n, Model::Linear_B)
    return sum(n .* Model.b_c)
end