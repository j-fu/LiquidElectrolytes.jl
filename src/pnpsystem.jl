"""
    pnpstorage(f, u, node, electrolyte)            

Finite volume storage term
"""
function pnpstorage(f, u, node, electrolyte)
    f[electrolyte.iϕ] = zero(eltype(u))
    f[electrolyte.ip] = zero(eltype(u))
    for ic in 1:(electrolyte.nc)
        f[ic] = u[ic]
    end
    return
end

"""
    pnpbstorage(f, u, node, electrolyte)

Finite volume boundary storage term
"""
function pnpbstorage(f, u, node, electrolyte)
    for ia in (electrolyte.nc + 1):(electrolyte.nc + electrolyte.na)
        f[ia] = u[ia]
    end
    return
end

"""
    pnpreaction(f, u, node, electrolyte)            

Finite volume reaction term
"""
function pnpreaction(f, u, node, electrolyte)
    ## Charge density
    f[electrolyte.iϕ] = -chargedensity(u, electrolyte)
    if solvepressure(electrolyte)
        f[electrolyte.ip] = 0
    else
        f[electrolyte.ip] = u[electrolyte.ip]
    end

    for ic in 1:(electrolyte.nc)
        f[ic] = 0
    end
    return
end


"""
     dμex(γk, γl, electrolyte)

Calculate differences of excess chemical potentials from activity coefficients
"""
@inline function dμex(γk, γl, electrolyte)
    return (rlog(γk) - rlog(γl)) * (electrolyte.RT)
end

"""
    sflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte; evelo=0.0)

 Sedan flux,  see Gaudeul/Fuhrmann 2022

Appearantly first described by Yu, Zhiping  and Dutton, Robert, SEDAN III, www-tcad.stanford.edu/tcad/programs/sedan3.html

 see also the 198? Fortran code available via
 https://web.archive.org/web/20210518233152/http://www-tcad.stanford.edu/tcad/programs/oldftpable.html

Verification calculation is in the paper.
"""
function sflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte; evelo = 0.0)
    (; D, z, F, RT) = electrolyte
    bp, bm = fbernoulli_pm(z[ic] * dϕ * F / RT + dμex(γk, γl, electrolyte) / RT - evelo / D[ic])
    return D[ic] * (bm * ck - bp * cl)
end

#=

B(b-a)/B(a-b)= (b-a)*(exp(a-b)-1)/ (a-b)*(exp(b-a)-1)
             = (1-exp(a-b))/(exp(b-a)-1)
             = exp(-b)*(exp(b)-exp(a)) / exp(-a)*(exp(b)-exp(a))
             = exp(a)/exp(b)

c/barc= C
ck/cl = bp/bm = exp(z ϕk*F/RT + μex_k/RT)/exp(z ϕl*F/RT + μex_l/RT)
=#

"""
    aflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte; evelo=0)

Flux expression based on  activities, see Fuhrmann, CPC 2015
??? Do we need to divide the velocity by the inverse activity coefficient ?

"""
function aflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte; evelo = 0.0)
    (; D, z, F, RT) = electrolyte
    Dx = D[ic] * (1 / γk + 1 / γl) / 2
    bp, bm = fbernoulli_pm(z[ic] * dϕ * F / RT - evelo / D[ic])
    return Dx * (bm * ck * γk - bp * cl * γl)
end

#=
ck/cl= bp/betaK  / bm/betal
     =  exp(z\phik *F/RT + \muexK/RT)/ exp(z\phil *F/RT + \muexL/RT)
=#

"""
    cflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte; evelo = 0)

Flux expression based on central differences, see Gaudeul/Fuhrmann 2022
"""
function cflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte; evelo = 0.0)
    (; D, z, F, RT) = electrolyte
    μk = rlog(ck) * RT
    μl = rlog(cl) * RT
    return D[ic] * 0.5 * (ck + cl) * ((μk - μl + dμex(γk, γl, electrolyte) + z[ic] * F * dϕ) / RT - evelo / D[ic])
end
#=

 ck/cl = exp(\muex_k + zF\phiK)/exp(\muex_l + zF\phil)
=#
"""
    pnpflux(f, u, edge, electrolyte)

Finite volume flux. It calls either [`sflux`](@ref), [`cflux`](@ref) or [`aflux`](@ref).
"""
function pnpflux(f, u, edge, electrolyte)
    iϕ = electrolyte.iϕ # index of potential
    ip = electrolyte.ip
    (; ip, iϕ, v0, v, M0, M, κ, ε_0, ε, RT, nc, eneutral, pscale, p_bulk, scheme) = electrolyte
    evelo = edgevelocity(electrolyte, edge.index)

    pk, pl = u[ip, 1] * pscale - p_bulk, u[ip, 2] * pscale - p_bulk
    ϕk, ϕl = u[iϕ, 1], u[iϕ, 2]

    @views qk, ql = chargedensity(u[:, 1], electrolyte), chargedensity(u[:, 2], electrolyte)
    @views c0k, bar_ck = c0_barc(u[:, 1], electrolyte)
    @views c0l, bar_cl = c0_barc(u[:, 2], electrolyte)

    dϕ = ϕk - ϕl

    f[iϕ] = ε * ε_0 * dϕ * !eneutral

    if solvepressure(electrolyte)
        f[ip] = u[ip, 1] - u[ip, 2] + (qk + ql) * dϕ / (2 * pscale)
    end

    for ic in 1:nc
        ck, cl = u[ic, 1], u[ic, 2]
        Mrel = M[ic] / M0 + κ[ic]
        barv = v[ic] + κ[ic] * v0
        tildev = barv - Mrel * v0
        γk = rexp(tildev * pk / RT) * (bar_ck / c0k)^Mrel * (1 / (v0 * bar_ck))
        γl = rexp(tildev * pl / RT) * (bar_cl / c0l)^Mrel * (1 / (v0 * bar_cl))

        if scheme == :μex
            f[ic] = sflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte; evelo)
        elseif electrolyte.scheme == :act
            f[ic] = aflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte; evelo)
        elseif electrolyte.scheme == :cent
            f[ic] = cflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte; evelo)
        else
            error("no such scheme: $(scheme)")
        end
    end
    return
end

"""
    PNPSystem(grid;
             celldata=ElectrolyteData(),
             bcondition=(f, u, n, e)-> nothing,
             reaction=(f, u, n, e)-> nothing,
             kwargs...)

Create VoronoiFVM system for generalized Poisson-Nernst-Planck. Input:
- `grid`: discretization grid
- `celldata`: composite struct containing electrolyte data
- `bcondition`: boundary condition
- `reaction` : reactions of the bulk species
- `kwargs`: Keyword arguments of VoronoiFVM.System
"""
function PNPSystem(
        grid;
        celldata = ElectrolyteData(),
        bcondition = (f, u, n, e) -> nothing,
        reaction = (f, u, n, e) -> nothing,
        kwargs...
    )

    function _pnpreaction(f, u, node, electrolyte)
        pnpreaction(f, u, node, electrolyte)
        reaction(f, u, node, electrolyte)
        return nothing
    end

    sys = VoronoiFVM.System(
        grid;
        data = celldata,
        flux = pnpflux,
        reaction = _pnpreaction,
        storage = pnpstorage,
        bcondition,
        species = [1:(celldata.nc)..., celldata.iϕ, celldata.ip],
        kwargs...
    )
    for ia in (celldata.nc + 1):(celldata.nc + celldata.na)
        enable_boundary_species!(sys, ia, [celldata.Γ_we])
    end
    return sys
end

"""
    electrolytedata(sys)
Extract electrolyte data from system.
"""
electrolytedata(sys) = sys.physics.data

"""
    pnpunknowns(sys)

Return vector of unknowns initialized with bulk data.
"""
function pnpunknowns(sys)
    (; iϕ, ip, nc, na, c_bulk, Γ_we) = electrolytedata(sys)
    u = unknowns(sys)
    @views u[iϕ, :] .= 0
    @views u[ip, :] .= 0
    for ic in 1:nc
        @views u[ic, :] .= c_bulk[ic]
    end
    for ia in (nc + 1):(nc + na)
        @views u[ia, :] .= 0
    end
    return u
end
