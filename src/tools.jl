myround(x; kwargs...) = round(x; kwargs...)
myround(s::Symbol; kwargs...) = s
myround(i::Int; kwargs...) = i
myround(b::Bool; kwargs...) = b

function showstruct(io::IO, this)
    for name in fieldnames(typeof(this))
        println(io, "$(lpad(name,20)) = $(myround.(getfield(this,name),sigdigits=5))")
    end
end


"""
    rexp(x;trunc=500.0)

Regularized exponential. Linear continuation for `x>trunc`,  
returns 1/rexp(-x) for `x<-trunc`.
"""
function rexp(x; trunc = 500.0)
    return if x < -trunc
        1.0 / rexp(-x; trunc)
    elseif x <= trunc
        exp(x)
    else
        exp(trunc) * (x - trunc + 1)
    end
end

"""
    rlog(u; eps=1.0e-40)

Regularized logarithm. Smooth linear continuation for `x<eps`.
This means we can calculate a "logarithm"  of a small negative number.
"""
function rlog(x; eps = 1.0e-40)
    if x < eps
        return log(eps) + (x - eps) / eps
    else
        return log(x)
    end
end
