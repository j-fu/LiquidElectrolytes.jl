"""
    myround(x; kwargs...)

Rounding for use in [`showstruct`](@ref).
"""
function myround end

myround(x; kwargs...) = round(x; kwargs...)
myround(x::Vector; kwargs...) = round.(x; kwargs...)
myround(s::Symbol; kwargs...) = ":$(s)"
myround(i::Int; kwargs...) = i
myround(b::Bool; kwargs...) = b
myround(::Nothing; kwargs...) = "nothing"
myround(f::Function; kwargs...) = string(f)
myround(c::DiffCache; kwargs...) = typeof(c)


"""
    showstruct(io::IO, this)

Print struct with field names and field values.
"""
function showstruct(io::IO, this)
    println(io, "$(typeof(this))(")
    for name in fieldnames(typeof(this))
        println(io, "$name = $(myround(getfield(this, name), sigdigits = 5)), ")
    end
    println(io, ")")
    return nothing
end


"""
    RExp(trunc)

Callable struct for regularized exponential. Linear continuation for `x>trunc`,  
returns 1/rexp(-x) for `x<-trunc`. Objects `myexp::RExp`  of this type are meant to replace
the exponential function and allow to invoke `myexp(x)`.
"""
struct RExp{T <: AbstractFloat}
    trunc::T
end

function (myexp::RExp)(x)
    (; trunc) = myexp
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
RExp(::Type{T}) where {T} = RExp{T}(log(sqrt(floatmax(T))))

"""
    RExp()
Return RExp(Float64)
"""
RExp() = RExp(Float64)

"""
    RLog(eps)

Callable struct for regularized  logarithm. Smooth linear continuation for `x<eps`.
This means we can calculate a "logarithm"  of a small negative number.
Objects `mylog::RLog`  of this type are meant to replace
the logarithm function and allow to invoke `mylog(x)`.
"""
struct RLog{T <: AbstractFloat}
    eps::T
end

function (mylog::RLog)(x)
    (; eps) = mylog
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
RLog(::Type{T}) where {T} = RLog{T}(eps(T)^4)

"""
    RLog()
Return RLog(Float64)
"""
RLog() = RLog(Float64)


"""
    _splitz(range::AbstractRange)

If range contains zero, split it into two parts, one with values <=0 and one with values >=0.
Otherwise, return the range or its reverse, such that first value always is the one with the smallest absolute value.
"""
function _splitz(range::AbstractRange)
    if range[1] >= 0
        return [range]
    elseif range[end] <= 0
        return [reverse(range)]
    else
        [0:(-step(range)):range[1], 0:step(range):range[end]]
    end
end

"""
    _splitz(range::Vector)

Version of [`_splitz(range::AbstractRange)`](@ref) for vectors.
"""
function _splitz(range::Vector)
    if range[1] >= 0
        return [vcat([0.0], range)]
    elseif range[end] <= 0
        return [vcat([0.0], reverse(range))]
    else
        for i in 1:length(range)
            if range[i] ≈ 0.0
                return [reverse(range[1:i]), range[i:end]]
            elseif i > 1 && range[i - 1] < 0.0 && range[i] > 0.0
                return [vcat([0.0], reverse(range[1:(i - 1)])), vcat([0.0], range[i:end])]
            end
        end
    end
end
