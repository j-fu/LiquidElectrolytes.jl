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
struct RExp{T<:AbstractFloat} <: RFunction
    trunc::T
end

"""
    (::RExp)(x)

Evaluation of functor at `x`
"""
function (rexp::RExp)(x)
    (;trunc) = rexp
    if x < -trunc
        return 1.0 / rexp(-x)
    elseif x <= trunc
        return exp(x)
    else
        return exp(trunc) * (x - trunc + 1)
    end
end

"""
    RExp(T::type)

Return `RExp(log(sqrt(floatmax(T))))`
"""
RExp(::Type{T}) where T=RExp{T}(log(sqrt(floatmax(T))))

"""
    RExp()
Return RExp(Float64)
"""
RExp() = RExp(Float64)

"""
    rexp(x::Any)

Exponential function used in LiquidElectrolytes.jl.
Calls `Base.exp` by default.

This can be overwritten  by defining e.g, `LiquidElectrolytes.rexp(x::Number)=myexp(x)`
where `myexp=RExp()`.
"""
rexp(x::Any)=Base.exp(x)


"""
    RLog(eps)

Functor struct for regularized  logarithm. Smooth linear continuation for `x<eps`.
This means we can calculate a "logarithm"  of a small negative number.
Objects of this type are meant to replace the logarithm function.
"""
struct RLog{T<:AbstractFloat} <: RFunction
    eps::T
end

function (rlog::RLog)(x)
    (;eps)=rlog
    if x < eps
        return log(eps) + (x - eps) / eps
    else
        return log(x)
    end
end



"""
    RLog(T::type)

Return `RLog(eps(T)^4)`
"""
RLog(::Type{T}) where T=RLog{T}(eps(T)^4)

"""
    RLog()
Return RLog(Float64)
"""
RLog() = RLog(Float64)

"""
    rlog(x::Any)

Logarithm function used in LiquidElectrolytes.jl.
Calls `Base.log` by default.

This can be overwritten  by defining e.g, `LiquidElectrolytes.rlog(x::Number)=mylog(x)`
where `mylog=RLog()`.
"""
rlog(x::Any)=Base.log(x)
