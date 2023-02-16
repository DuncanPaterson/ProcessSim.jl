

struct SoaveAlphaPolynomial{TMatrix,TVector} <: AlphaModel
    c::TMatrix
    T_c::TVector
end

function Alpha(T, Model::SoaveAlphaPolynomial)

    tr  = sqrt.(T ./ Model.T_c)
    one_min_tr = 1.0 .- tr
    alpha = similar(tr)
    alpha .= 1.0

    for i in range(1, size(Model.c, 2))
        ci = Model.c[:, i]
        alpha .= alpha .+ ci .* one_min_tr.^i
    end

    alpha .= alpha.^2
    return alpha
end


function SoaveAlpham(ω)
    return 0.48 + 1.574 * ω - 0.175 * ω^2
end

function PRAlpham(ω)
    return 0.17464 + 1.54226 * ω - 0.26992 * ω^2
end