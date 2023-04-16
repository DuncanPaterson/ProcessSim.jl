function QuadraticRoots(b, c)
    if (c[1] == 0)
        neg_c_inv_a = -c/a
        if (neg_c_inv_a < 0)
            return SVector{2, typeof(b)}(NaN, NaN)
        end
        res = sqrt(neg_c_inv_a)
        return SVector{2, typeof(b)}(-res, res)
    end
    discriminant = b^2 - 4 * c
    x1 = (-b - sqrt(discriminant)) / 2
    x2 = (-b + sqrt(discriminant)) / 2
    if (x1 < x2)
        return SVector{2, typeof(b)}(x1, x2)
    else
        return SVector{2, typeof(b)}(x2, x1)
    end
end