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
    default_reaction(f, u, bnode, electrolyte)

Default reaction amounts to `nothing`
"""
default_reaction(f, u, node, electrolyte) = nothing


"""
     dμex(γk, γl, electrolyte)

Calculate differences of excess chemical potentials from activity coefficients
"""
@inline function dμex(γk, γl, electrolyte)
    if γk > γl
        rlog(γk / γl, electrolyte) * (electrolyte.RT)
    else
        -rlog(γl / γk, electrolyte) * (electrolyte.RT)
    end
end

"""
    sflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte)

 Sedan flux,  see Gaudeul/Fuhrmann 2022

Appearantly first described by Yu, Zhiping  and Dutton, Robert, SEDAN III, www-tcad.stanford.edu/tcad/programs/sedan3.html

 see also the 198? Fortran code available via
 https://web.archive.org/web/20210518233152/http://www-tcad.stanford.edu/tcad/programs/oldftpable.html

Verification calculation is in the paper.
"""
function sflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte)
    (; D, z, F, RT) = electrolyte
    bp, bm = fbernoulli_pm(z[ic] * dϕ * F / RT + dμex(γk, γl, electrolyte) /RT)
    D[ic]*(bm * ck - bp * cl)
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
    aflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte)

Flux expression based on  activities, see Fuhrmann, CPC 2015
"""
function aflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte)
    (; D, z, F, RT) = electrolyte
    bp, bm = fbernoulli_pm(z[ic] * dϕ * F / RT)
    D[ic] * (bm * ck * γk - bp * cl * γl) * (1 / γk + 1 / γl) / 2
end

#=
ck/cl= bp/betaK  / bm/betal
     =  exp(z\phik *F/RT + \muexK/RT)/ exp(z\phil *F/RT + \muexL/RT)
=#

"""
    cflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte)

Flux expression based on central differences, see Gaudeul/Fuhrmann 2022
"""
function cflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte)
    (; D, z, F, RT) = electrolyte
    μk = rlog(ck, electrolyte) * RT
    μl = rlog(cl, electrolyte) * RT
    D[ic] * 0.5 * (ck+ cl) * (μk - μl +  dμex(γk, γl, electrolyte) + z[ic] * F * dϕ) / RT
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

    pk, pl = u[ip, 1] * pscale-p_bulk, u[ip, 2] * pscale-p_bulk
    ϕk, ϕl = u[iϕ, 1], u[iϕ, 2]

    @views qk, ql = charge(u[:, 1], electrolyte), charge(u[:, 2], electrolyte)
    @views c0k, bar_ck = c0_barc(u[:, 1], electrolyte)
    @views c0l, bar_cl = c0_barc(u[:, 2], electrolyte)

    dϕ = ϕk - ϕl
    dp = pk - pl

    f[iϕ] = ε * ε_0 * dϕ * !eneutral
    f[ip] = dp + (qk + ql) * dϕ / 2

    γk, γl = 1.0, 1.0
    bikerman = !iszero(v)

    for ic = 1:nc
        f[ic] = 0.0
        ## Regularize ck,cl so they don't become zero
        ck, cl = u[ic, 1], u[ic, 2]
        barv = 0.0

        ## Calculate the  activity coefficients first,
        ## as these expressions are less degenerating.
        if bikerman
            Mrel = M[ic] / M0
            barv=v[ic] + κ[ic]*v0
            tildev=barv - Mrel*v0
            γk = exp(tildev * pk / (RT)) * (bar_ck / c0k)^Mrel*(1/bar_ck)
            γl = exp(tildev * pl / (RT)) * (bar_cl / c0l)^Mrel*(1/bar_cl)
        end

        if scheme == :μex
            f[ic] = sflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte)
        elseif electrolyte.scheme == :act
            f[ic] = aflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte)
        elseif electrolyte.scheme == :cent
            f[ic] = cflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte)
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

Create VoronoiFVM system for generalized Poisson-Nernst-Planck. Input:
- `grid`: discretization grid
- `celldata`: composite struct containing electrolyte data
- `bcondition`: boundary condition
- `reaction` : reactions of the bulk species
- `kwargs`: Keyword arguments of VoronoiFVM.System
"""
function PNPSystem(grid; celldata = ElectrolyteData(), bcondition = default_bcondition, reaction = default_reaction, kwargs...)

    function _pnpreaction(f, u, node, electrolyte)
        pnpreaction(f, u, node, electrolyte)
        reaction(f, u, node, electrolyte)
        nothing
    end

    sys = VoronoiFVM.System(grid;
                            data = celldata,
                            flux = pnpflux,
                            reaction = _pnpreaction,
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
