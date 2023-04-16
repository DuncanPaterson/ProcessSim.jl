function CubicRoots(a, b, c, d)
    roots = MVector{3,typeof(a)}(NaN, NaN, NaN)
    if (a == 0)
        roots[1], roots[2] = quadratic_roots(c/b, d/b)
        return roots
    else
        return CubicRoots(b/a, c/a, d/a)
    end
end

function CubicRoots(b, c, d)
    roots = MVector{3,typeof(b)}(NaN, NaN, NaN)
    if (d == 0) 
        roots[1], roots[2] = quadratic_roots([c[2], c[3]])
        roots[3] = 0
        sort!(roots)
        return roots
    end
    Q = (b * b - 3 * c) / 9
    R = (2 * b * b * b - 9 * b * c + 27 * d) / 54
    arg = R * R - Q^3
    if (arg < 0)
        rtQ = sqrt(Q)
        theta = acos(R / (Q * rtQ)) / 3
        st = sin(theta)
        ct = cos(theta)
        roots[1] = -2 * rtQ * ct - b / 3
        roots[2] = -rtQ * (-ct + sqrt(Real(3)) * st) - b / 3
        roots[3] = rtQ * (ct + sqrt(Real(3)) * st) - b / 3
    else
        A = (R >= 0 ? -1.0 : 1.0) * cbrt(abs(R) + sqrt(arg))
        B = 0.0
        if (A != 0) 
            B = Q / A
        end
        roots[1] = A + B - b / 3
        if (A == B || arg == 0.0) 
            roots[2] = -A - b / 3
            roots[3] = -A - b / 3
        end
    end
    # Polish roots
    a = typeof(b)(1)
    for r in roots
        f = fma(a, r, b)
        f = fma(f, r, c)
        f = fma(f, r, d)
        df = fma(3 * a, r, 2 * b)
        df = fma(df, r, c)
        if (df != 0) 
            d2f = fma(6 * a, r, 2 * b)
            denom = 2 * df * df - f * d2f
            if (denom != 0) 
                r -= 2 * f * df / denom
            else
                r -= f / df
            end
        end
    end
    sort!(roots)
    return roots
end