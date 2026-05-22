"""
    myround(x; kwargs...)

Rounding for use in [`showstruct`](@ref).
"""
function myround end

myround(x; kwargs...) = round(x; kwargs...)
myround(x::Vector; kwargs...) = myround.(x; kwargs...)
myround(s::Symbol; kwargs...) = ":$(s)"
myround(i::Int; kwargs...) = i
myround(b::Bool; kwargs...) = b
myround(::Nothing; kwargs...) = "nothing"
myround(f::Function; kwargs...) = string(f)
myround(c::DiffCache; kwargs...) = string(nameof(typeof(c)))


"""
    showstruct(io::IO, this)

Print struct with field names and field values.
"""
function showstruct(io::IO, this)
    println(io, "$(string(nameof(typeof(this))))(")
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

function (mylog::RLog)(x::Any)
    (; eps) = mylog
    return ifelse(x < eps, log(eps) + (x - eps) / eps, log(x))
end

function (mylog::RLog)(x::Union{AbstractFloat, ForwardDiff.Dual})
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
    splitc(range::AbstractRange; center=0)

If range contains `center`, split it into two parts, one with values <=center and one with values >=center.
Otherwise, return the range or its reverse, such that first value always is the one with the smallest absolute value.
"""
function splitc(range::AbstractRange; center = 0)
    if range[1] >= center
        return [range]
    elseif range[end] <= center
        return [reverse(range)]
    else
        [center:(-step(range)):range[1], center:step(range):range[end]]
    end
end

"""
    splitc(range::Vector; center=0)

Version of [`splitc(range::AbstractRange)`](@ref) for vectors.
"""
function splitc(range::Vector; center = 0)
    if range[1] >= center
        return [vcat([center], range)]
    elseif range[end] <= center
        return [vcat([center], reverse(range))]
    else
        for i in 1:length(range)
            if range[i] ≈ center
                return [reverse(range[1:i]), range[i:end]]
            elseif i > 1 && range[i - 1] < center && range[i] > center
                return [vcat([center], reverse(range[1:(i - 1)])), vcat([center], range[i:end])]
            end
        end
    end
end

# Large value of gap capacitance enforcing
# Dirichlet BC
const C_large = 1.0e15

myvalue(::Any) = 0.0
myvalue(x::ForwardDiff.Dual) = value(x)
myvalue(x::AbstractFloat) = x
