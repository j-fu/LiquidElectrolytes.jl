"""
    pnpstorage(f, u, node, electrolyte)            

Finite volume storage term
"""
function pnpstorage(f, u, node, electrolyte)
    f[electrolyte.iϕ] = zero(eltype(u))
    f[electrolyte.ip] = zero(eltype(u))
    for ic = 1:(electrolyte.nc)
        f[ic] = u[ic]
    end
end

"""
    pnpbstorage(f, u, node, electrolyte)

Finite volume boundary storage term
"""
function pnpbstorage(f, u, node, electrolyte)
    for ia = (electrolyte.nc + 1):(electrolyte.nc + electrolyte.na)
        f[ia] = u[ia]
    end
end

"""
    pnpreaction(f, u, node, electrolyte)            

Finite volume reaction term
"""
function pnpreaction(f, u, node, electrolyte)
    ## Charge density
    f[electrolyte.iϕ] = -charge(u, electrolyte)
    f[electrolyte.ip] = 0
    for ic = 1:(electrolyte.nc)
        f[ic] = 0
    end
end

"""
    default_bcondition(f,u,bnode,electrolyte)

Default boundary condition amounts to `nothing`
"""
default_bcondition(f, u, bnode, electrolyte) = nothing

"""
     dμex(βk, βl, electrolyte)

Calculate differences of excess chemical potentials from reciprocal activity coefficients
"""
@inline function dμex(βk, βl, electrolyte)
    if βk > βl
        rlog(βk / βl, electrolyte) * (electrolyte.RT)
    else
        -rlog(βl / βk, electrolyte) * (electrolyte.RT)
    end
end

"""
    sflux(ic,dϕ,ck,cl,βk,βl,bar_ck,bar_cl,electrolyte)

 Sedan flux,  see Gaudeul/Fuhrmann 2022

Appearantly first described by Yu, Zhiping  and Dutton, Robert, SEDAN III, www-tcad.stanford.edu/tcad/programs/sedan3.html

 see also the 198? Fortran code available via
 https://web.archive.org/web/20210518233152/http://www-tcad.stanford.edu/tcad/programs/oldftpable.html

Verification calculation is in the paper.
"""
function sflux(ic, dϕ, ck, cl, βk, βl, bar_ck, bar_cl, electrolyte)
    (; D, z, F, RT) = electrolyte
    bp, bm = fbernoulli_pm(z[ic] * dϕ * F / RT + dμex(βk, βl, electrolyte) / RT)
    D[ic] * (bm * ck - bp * cl)
end

"""
    aflux(ic,dϕ,ck,cl,βk,βl,bar_ck,bar_cl,electrolyte)

Flux expression based on reciprocal activity coefficents, see Fuhrmann, CPC 2015
"""
function aflux(ic, dϕ, ck, cl, βk, βl, bar_ck, bar_cl, electrolyte)
    (; D, z, F, RT) = electrolyte
    bp, bm = fbernoulli_pm(z[ic] * dϕ * F / RT)
    D[ic] * (bm * ck * βk - bp * cl * βl) * (1 / βk + 1 / βl) / 2
end

"""
    aflux(ic,dϕ,ck,cl,βk,βl,bar_ck,bar_cl,electrolyte)

Flux expression based on central differences, see Gaudeul/Fuhrmann 2022
"""
function cflux(ic, dϕ, ck, cl, βk, βl, bar_ck, bar_cl, electrolyte)
    (; D, z, F, RT) = electrolyte
    μk = rlog(ck / bar_ck, electrolyte) * RT
    μl = rlog(cl / bar_cl, electrolyte) * RT
    D[ic] * 0.5 * (ck + cl) * (μk - μl + dμex(βk, βl, electrolyte) + z[ic] * F * dϕ) / RT
end

"""
    pnpflux(f, u, edge, electrolyte)

Finite volume flux. It calls either [`sflux`](@ref), [`cflux`](@ref) or [`aflux`](@ref).
"""
function pnpflux(f, u, edge, electrolyte)
    iϕ = electrolyte.iϕ # index of potential
    ip = electrolyte.ip
    (; ip, iϕ, v0, v, M0, M, κ, ε_0, ε, RT, nc, eneutral, pscale, scheme) = electrolyte

    pk, pl = u[ip, 1] * pscale, u[ip, 2] * pscale
    ϕk, ϕl = u[iϕ, 1], u[iϕ, 2]

    @views qk, ql = charge(u[:, 1], electrolyte), charge(u[:, 2], electrolyte)
    @views c0k, bar_ck = c0_barc(u[:, 1], electrolyte)
    @views c0l, bar_cl = c0_barc(u[:, 2], electrolyte)

    dϕ = ϕk - ϕl
    dp = pk - pl

    f[iϕ] = ε * ε_0 * dϕ * !eneutral
    f[ip] = dp + (qk + ql) * dϕ / 2

    βk, βl = 1.0, 1.0
    bikerman = !iszero(v)

    for ic = 1:nc
        f[ic] = 0.0
        ## Regularize ck,cl so they don't become zero
        ck, cl = u[ic, 1], u[ic, 2]
        barv = 0.0

        ## Calculate the reciprocal activity coefficients first,
        ## as these expressions are less degenerating.
        if bikerman
            Mrel = M[ic] / M0
            barv = v[ic] + (κ[ic] - Mrel) * v0
            βk = exp(barv * pk / (RT)) * (bar_ck / c0k)^Mrel
            βl = exp(barv * pl / (RT)) * (bar_cl / c0l)^Mrel
        end

        if scheme == :μex
            f[ic] = sflux(ic, dϕ, ck, cl, βk, βl, bar_ck, bar_cl, electrolyte)
        elseif electrolyte.scheme == :act
            f[ic] = aflux(ic, dϕ, ck, cl, βk, βl, bar_ck, bar_cl, electrolyte)
        elseif electrolyte.scheme == :cent
            f[ic] = cflux(ic, dϕ, ck, cl, βk, βl, bar_ck, bar_cl, electrolyte)
        else
            error("no such scheme: $(scheme)")
        end
    end
end

"""
    PNPSystem(grid;
             celldata=ElectrolyteData(),
             bcondition=default_bcondition,
             kwargs...)

Create VoronoiFVM system. Input:
- `grid`: discretization grid
- `celldata`: composite struct containing electrolyte data
- `bcondition`: boundary condition
- `kwargs`: Keyword arguments of VoronoiFVM.System
"""
function PNPSystem(grid; celldata = ElectrolyteData(), bcondition = default_bcondition, kwargs...)
    sys = VoronoiFVM.System(grid;
                            data = deepcopy(celldata),
                            flux = pnpflux,
                            reaction = pnpreaction,
                            storage = pnpstorage,
                            bcondition,
                            species = [1:(celldata.nc)..., celldata.iϕ, celldata.ip],
                            kwargs...)
    for ia = (celldata.nc + 1):(celldata.nc + celldata.na)
        enable_boundary_species!(sys, ia, [celldata.Γ_we])
    end
    sys
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
    for ic = 1:nc
        @views u[ic, :] .= c_bulk[ic]
    end
    for ia = (nc + 1):(nc + na)
        @views u[ia, :] .= 0
    end
    u
end
