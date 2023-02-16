

struct SoaveAlphaPolynomial{TMatrix,TVector} <: AlphaModel
    c::TMatrix
    T_c::TVector
    a_c::TVector
end

function Alpha(T, Model::SoaveAlphaPolynomial)
    alpha = similar(Model.T_c)
    one_min_tr = similar(alpha)

    alpha .= 1.0
    one_min_tr .= 1.0 .- sqrt.(T ./ Model.T_c)

    for i in range(size(Model.c, 2))
        ci = Model.c[:, i]
        alpha .= alpha + ci * one_min_tr
    end

    alpha .= alpha .* alpha

    return alpha
end


function SoaveAlpham(ω)
    return 0.48 + 1.574 * ω - 0.175 * ω^2
end

function PRAlpham(ω)
    return 0.17464 + 1.54226 * ω - 0.26992 * ω^2
end