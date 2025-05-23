"""
    pnpstorage!(f, u, node, electrolyte)            

Finite volume storage term
"""
function pnpstorage!(f, u, node, electrolyte)
    (; nc) = electrolyte
    for ic in 1:nc
        f[ic] = u[ic]
    end
    return
end

"""
    pnpbstorage!(f, u, node, electrolyte)

Finite volume boundary storage term
"""
function pnpbstorage!(f, u, node, electrolyte)
    (; nc, na) = electrolyte
    for ia in (nc + 1):(nc + na)
        f[ia] = u[ia]
    end
    return
end

"""
    pnpreaction!(f, u, node, electrolyte)            

Finite volume reaction term
"""
function pnpreaction!(f, u, node, electrolyte)
    (; iϕ, ip) = electrolyte
    f[iϕ] = -chargedensity(u, electrolyte)
    if solvepressure(electrolyte)
        f[ip] = 0
    else
        f[ip] = u[ip]
    end
    return
end


"""
     dμex(γk, γl, electrolyte)

Calculate differences of excess chemical potentials from activity coefficients
"""
@inline function dμex(γk, γl, electrolyte)
    (; RT, rlog) = electrolyte
    return (rlog(γk) - rlog(γl)) * RT
end

"""
    μex_flux!(f,dϕ,ck,cl,γk,γl,electrolyte; evelo=0.0)

Excess chemical potential (Sedan) upweind flux, based on modified
Scharfetter-Gummel scheme, see Gaudeul/Fuhrmann 2022.

Appearantly first described by Yu, Zhiping  and Dutton, Robert, SEDAN III, www-tcad.stanford.edu/tcad/programs/sedan3.html

see also the 198? Fortran code available via
https://web.archive.org/web/20210518233152/http://www-tcad.stanford.edu/tcad/programs/oldftpable.html

Verification calculation is in the paper.
"""
function μex_flux!(f, dϕ, ck, cl, γk, γl, electrolyte; evelo = 0.0)
    (; D, z, F, RT, nc) = electrolyte
    for ic in 1:nc
        bp, bm = fbernoulli_pm(z[ic] * dϕ * F / RT + dμex(γk[ic], γl[ic], electrolyte) / RT - evelo / D[ic])
        f[ic] = D[ic] * (bm * ck[ic] - bp * cl[ic])
    end
    return nothing
end


"""
    act_flux!(ic,dϕ,ck,cl,γk,γl,electrolyte; evelo=0)

Scharfetter-Gummel upwind flux expression based on  activities, see Fuhrmann, CPC 2015.
As shown in Chainais-Hilliaret, Cances, Fuhrmann, Gaudeul, 2019, this is
consistent to thermodynamic equilibrium but shows inferior convergence behavior.
"""
function act_flux!(f, dϕ, ck, cl, γk, γl, electrolyte; evelo = 0.0)
    (; D, z, F, RT, nc) = electrolyte
    for ic in 1:nc
        Dx = D[ic] * (1 / γk[ic] + 1 / γl[ic]) / 2
        bp, bm = fbernoulli_pm(z[ic] * dϕ * F / RT - evelo / D[ic])
        f[ic] = Dx * (bm * ck[ic] * γk[ic] - bp * cl[ic] * γl[ic])
    end
    return nothing
end

"""
    cent_flux!(ic,dϕ,ck,cl,γk,γl,electrolyte; evelo = 0)

Flux expression based on central differences, see Gaudeul/Fuhrmann 2022.
As shown in the paper  this is
covergent,  consistent to thermodynamic equilibrium but may show 
inferior convergence behavior.
"""
function cent_flux!(f, dϕ, ck, cl, γk, γl, electrolyte; evelo = 0.0)
    (; D, z, F, RT, nc, rlog) = electrolyte
    for ic in 1:nc
        μk = rlog(ck[ic]) * RT
        μl = rlog(cl[ic]) * RT
        f[ic] = D[ic] * 0.5 * (ck[ic] + cl[ic]) * ((μk - μl + dμex(γk[ic], γl[ic], electrolyte) + z[ic] * F * dϕ) / RT - evelo / D[ic])
    end
    return nothing
end

"""
    pnpflux(f, u, edge, electrolyte)

Finite volume flux. It calls either [`μex_flux!`](@ref), [`cent_flux!`](@ref) or [`act_flux!`](@ref).
"""
function pnpflux!(f, u, edge, electrolyte)
    (;
        nc, ip, iϕ,
        ε_0, ε,
        eneutral, pscale, p_bulk,
        upwindflux!, actcoeff!,
        γk_cache, γl_cache,
    ) = electrolyte

    evelo = edgevelocity(electrolyte, edge.index)

    pk, pl = u[ip, 1] * pscale - p_bulk, u[ip, 2] * pscale - p_bulk
    ϕk, ϕl = u[iϕ, 1], u[iϕ, 2]
    ck, cl = view(u, :, 1), view(u, :, 2)
    qk, ql = chargedensity(ck, electrolyte), chargedensity(cl, electrolyte)
    γk, γl = get_tmp(γk_cache, u), get_tmp(γl_cache, u)
    actcoeff!(γk, ck, pk, electrolyte)
    actcoeff!(γl, cl, pl, electrolyte)

    dϕ = ϕk - ϕl
    f[iϕ] = ε * ε_0 * dϕ * !eneutral
    if solvepressure(electrolyte)
        f[ip] = u[ip, 1] - u[ip, 2] + (qk + ql) * dϕ / (2 * pscale)
    end

    upwindflux!(f, dϕ, ck, cl, γk, γl, electrolyte; evelo)

    return
end

struct PNPSystem <: AbstractElectrochemicalSystem
    vfvmsys::VoronoiFVM.System
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
        grid::ExtendableGrid;
        celldata = ElectrolyteData(),
        bcondition = (f, u, n, e) -> nothing,
        reaction = (f, u, n, e) -> nothing,
        kwargs...
    )
    update_derived!(celldata)

    function _pnpreaction!(f, u, node, electrolyte)
        pnpreaction!(f, u, node, electrolyte)
        reaction(f, u, node, electrolyte)
        return nothing
    end

    sys = VoronoiFVM.System(
        grid;
        data = celldata,
        flux = pnpflux!,
        reaction = _pnpreaction!,
        storage = pnpstorage!,
        bcondition,
        species = [1:(celldata.nc)..., celldata.iϕ, celldata.ip],
        kwargs...
    )
    for ia in (celldata.nc + 1):(celldata.nc + celldata.na)
        enable_boundary_species!(sys, ia, [celldata.Γ_we])
    end
    return PNPSystem(sys)
end

"""
    electrolytedata(sys)
Extract electrolyte data from system.
"""
electrolytedata(sys::AbstractElectrochemicalSystem) = sys.vfvmsys.physics.data

"""
    unknowns(sys)

Return vector of unknowns initialized with bulk data.
"""
function VoronoiFVM.unknowns(esys::AbstractElectrochemicalSystem)
    sys = esys.vfvmsys
    (; iϕ, ip, nc, na, c_bulk, Γ_we) = electrolytedata(esys)
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
