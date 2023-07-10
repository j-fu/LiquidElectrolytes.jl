const NCMAX = 10
using StaticArrays

"""
      sum y_i=1
      sum v_i c_i=1
      cbar=1/sum v_i y_i
"""

molfrac(μ, data) = exp(μ / R * data.T)

function conc!(c, u, data)
    y0 = 1
    sumvy = 0
    for i = 1:(data.nc)
        c[i] = molfrac(u[i], data)
        sumvy += c[i] * data.v[i]
        y0 -= c[i]
    end
    sumvy += y0 * data.v0
    cbar = 1 / sumvy
    for i = 1:(data.nc)
        c[i] *= cbar
    end
end

function cpstorage(f, u, node, data)
    f[data.iϕ] = zero(eltype(u))
    conc!(f, u, data)
end

function charge(c, data)
    q = zero(eltype(c))
    for ic = 1:(data.nc)
        q += c[ic] * data.z[ic]
    end
    q * F
end

charge(u, i, data) = @views charge(u[:, i], data)

function nppreaction(f, u, node, data)
    ## Charge density
    conc!(f, u, data)
    f[data.iϕ] = -charge(f, data)
    for ic = 1:(data.nc)
        f[ic] = 0
    end
end

default_bcondition(f, u, bnode, data) = nothing

"""
after "centered flux"
"""
function cpcflux(f, u, edge, data)
    iϕ = data.iϕ # index of potential
    ip = data.ip
    ## Poisson flux
    dϕ = u[iϕ, 1] - u[iϕ, 2]
    dp = u[ip, 1] - u[ip, 2]

    muavg = @MVector zeros(NCMAX)
    cavg = @MVector zeros(NCMAX)
    for ic = 1:nc
        muavg[i] = 0.5 * (u[ic, 1] + u[ic, 2])
    end
    conc!(cavg, muavg, data)
    q = charge(cavg, data)
    f[iϕ] = data.ε * ε_0 * dϕ
    f[ip] = dp + q * dϕ

    mu0k = 1.0
    mu0l = 1.0
    for ic = 1:nc
        mu0k -= 1.0 - molfrac(u[ic, 1], data)
        mu0l -= 1.0 - molfrac(u[ic, 2], data)
    end
    mu0k = log(mu0k)
    mu0l = log(mu0l)
    for ic = 1:(data.nc)
        if !iszero(data.v)
            ## Caclculate excess chemical potentials
            muk = -log_bar_ck - (log_c0k - log_bar_ck) * data.v[ic] / data.v0
            mul = -log_bar_cl - (log_c0l - log_bar_cl) * data.v[ic] / data.v0
        else
            muk = 0.0
            mul = 0.0
        end

        ## Combine potential gradient and excess chemical gradient
        arg = data.Z[ic] * (u[iϕ, 1] - u[iϕ, 2]) + (muk - mul) / (R * data.T)

        ## Call Bernoulli function
        bp, bm = fbernoulli_pm(arg)

        ## Calculate drift-diffusion flux
        f[ic] = data.D[ic] * (bm * u[ic, 1] - bp * u[ic, 2])
    end
end

# for ic = 1:data.nc
#     if !iszero(data.v)
#         ## Caclculate excess chemical potentials
#         muk = -log_bar_ck - (log_c0k - log_bar_ck) * data.v[ic] / data.v0
#         mul = -log_bar_cl - (log_c0l - log_bar_cl) * data.v[ic] / data.v0
#     else
#         muk=0.0
#         mul=0.0
#     end

#     ## Combine potential gradient and excess chemical gradient
#     arg = data.Z[ic] * (u[iϕ,1] - u[iϕ,2]) + (muk - mul)/(R*data.T)

#     ## Call Bernoulli function
#     bp, bm = fbernoulli_pm(arg)

#     ## Calculate drift-diffusion flux
#     f[ic] = data.D[ic] * (bm * u[ic,1] - bp * u[ic,2])
# end

function NPPSystem(grid; electrolyte = nothing, bcondition = default_bcondition, kwargs...)
    sys = VoronoiFVM.System(grid;
                            data = electrolyte,
                            flux = nppflux,
                            reaction = nppreaction,
                            storage = nppstorage,
                            bcondition,
                            species = [1:(electrolyte.nc)..., electrolyte.iϕ, electrolyte.ip],
                            kwargs...)
end

electrolytedata(sys)::ElectrolyteData = sys.physics.data

boundarydata(sys) = sys.physics.data.bdata

function nppunknowns(sys)
    (; iϕ, ip, nc, c_bulk) = electrolytedata(sys)
    u = unknowns(sys)
    @views u[iϕ, :] .= 0
    @views u[ip, :] .= 0
    for ic = 1:nc
        @views u[ic, :] .= c_bulk[ic]
    end
    u
end
