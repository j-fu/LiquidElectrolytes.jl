abstract type RFunction<:Function end

myround(x; kwargs...) = round(x; kwargs...)
myround(s::Symbol; kwargs...) = s
myround(i::Int; kwargs...) = i
myround(b::Bool; kwargs...) = b
myround(f::Function; kwargs...) = string(f)


function showstruct(io::IO, this)
    for name in fieldnames(typeof(this))
        println(io, "$(lpad(name, 20)) = $(myround.(getfield(this, name), sigdigits = 5))")
    end
    return
end

"""
    RExp(trunc)

Functor struct for regularized exponential. Linear continuation for `x>trunc`,  
returns 1/rexp(-x) for `x<-trunc`. Objects of this type are meant to replace
the exponential function.
"""
Base.@kwdef struct RExp{T<:AbstractFloat} <: RFunction
    trunc::T=500.0
end

function (rexp::RExp)(x)
    (;trunc) = rexp
    if x < -trunc
        return 1.0 / rexp(-x; trunc)
    elseif x <= trunc
        return exp(x)
    else
        return exp(trunc) * (x - trunc + 1)
    end
end



"""
    RLog(eps)

Functor struct for regularized  logarithm. Smooth linear continuation for `x<eps`.
This means we can calculate a "logarithm"  of a small negative number.
Objects of this type are meant to replace the logarithm function.
"""
Base.@kwdef struct RLog{T<:AbstractFloat} <: RFunction
    eps::T=1.0e-40
end

function (rlog::RLog)(x)
    (;eps)=rlog
    if x < eps
        return log(eps) + (x - eps) / eps
    else
        return log(x)
    end
end
