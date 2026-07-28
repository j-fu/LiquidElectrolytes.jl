### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 93a15e56-7ba2-4c77-951f-6c08ecef0d5d
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "docs"))
    using Revise
    using CairoMakie
    using GridVisualize
    default_plotter!(CairoMakie)
    using Colors
    set_theme!(theme_dark())
    using Test
    using LiquidElectrolytes
    using ExtendableGrids
    using LessUnitful
    using PlutoUI
    using VoronoiFVM
    using LinearAlgebra
    VoronoiFVM.log_output!()
end

# ╔═╡ de0e3435-2736-4f75-8706-aade9016bc1f
begin
    scanrate = 10
    Ls = [0.1, 0.5, 1, 1.5] * ufac"mm"
    nperiods = 2
    C_elec = 0.01ufac"mol/dm^3"
    C_O = 1.0e-3 * ufac"mol / dm^3"
    C_R = 1.0e-10 * ufac"mol / dm^3"
    ircompfactor = 0.95
    z = [3, 2, 1, -1]
    f_κ = 2
    vmin = -0.5ufac"V"
    vmax = 0.1ufac"V"
    const k0::Float64 = 13.5 * 100 + ufac"cm / s"
    const α = 0.45
    const E0::Float64 = -0.173ufac"V"
    const k0::Float64 = 13.5 * ufac"cm / s"
    const iO = 1
    const iR = 2
    damp_initial = 0.1
    tol_round = 1.0e-9
    max_round = 3
end

# ╔═╡ bf07c780-2c90-4bce-a18a-3f9a7514dfc2
begin
    ε1 = 6.0 * ph"VacuumElectricPermittivity"
    ε2 = 30.0 * ph"VacuumElectricPermittivity"
    x2 = 0.29ufac"nm"
    x3 = 0.59ufac"nm"
    d1 = x2
    d2 = x3 - x2
    C_gap = 1.0 / (d1 / ε1 + d2 / ε2)
    C_gap / ufac"μF/cm^2"
end

# ╔═╡ 08470000-576c-478c-9528-878552cec9a5
md"""
Sawtooth voltage control:
"""

# ╔═╡ d956f42f-f5c4-401e-873b-d076fb27a887
sawtooth = SawTooth(; scanrate, vmin, vmax)

# ╔═╡ e28357fd-29b5-44e3-abb3-182f9a4ff936
md"""
Dielectric decrement: this is in the moment an experiemental feature due to lack of an implementation of a better model.
"""

# ╔═╡ e5ab8c3c-3614-4add-aab9-689e8720d1c3
function ε_dec(x)
    if x[1] < 1.0e-9
        return 0.01
    else
        return 1.0
    end
end

# ╔═╡ ffb2f6e4-4a9c-4b9e-95a8-4bfb7150e502
md"""
Redox rection function using the rate expressions according to Landstorfer et al.:
"""

# ╔═╡ a5b7c828-c97d-49eb-bd7f-16462c9ebd50
function redoxreaction_dirichlet(f, u, bnode, data)
    (; ip, iϕ, v, v0, F, RT, κ) = data
    c0, barc = c0_barc(u, data)
    μR = chemical_potential(u[iR], barc, u[ip], v[iR] + κ[iR] * v0, data)
    μO = chemical_potential(u[iO], barc, u[ip], v[iO] + κ[iO] * v0, data)
    A = (μR - μO - F * E0) / RT
    j_O = rrate(k0, α, A)
    f[iO] -= j_O
    f[iR] += j_O
    return
end

# ╔═╡ abc2f8e7-8357-47c4-bb3e-bbbbb5b816e9
function redoxreaction_robin(y, u, bnode, data)
    (; iϕ, F, RT) = data
    # E is the applied electrode potential (vs. the reference electrode).
    ϕ_PET = u[iϕ]
    ϕ_L = 0
    E = u[iϕ]
    iO = 1
    iR = 2
    η = (E - E0 - (ϕ_PET - ϕ_L)) * F / RT
    k_f = k0 * exp(-α * η)
    k_b = k0 * exp((1 - α) * η)
    j_O = k_b * u[iR] - k_f * u[iO]
    y[iO] -= j_O
    y[iR] += j_O

    return
end

# ╔═╡ 30698409-7cbd-4328-937a-e76f26d15151
redoxreaction = redoxreaction_robin

# ╔═╡ c32610f4-88e2-499d-abbc-03d6513806d3
edata_unc = ElectrolyteData(;
    nc = 4,
    z,
    c_bulk = [C_O, C_R, C_elec, C_elec + z[1] * C_O + z[2] * C_R],
    κ = [10, 10, 10, 10.0] * f_κ,
    D = [
        7.5e-6,
        10.4e-6,
        2.0e-5,
        1.9e-5,
    ] * ufac"cm^2/s",
    ε = 80.0,
    T = 293.15 * ufac"K",
    # Position of "pseudo reference electrode" for pseudopotentiostat approach
    x_ref = [10.0 * ufac"nm"],
    # Dielectric decrement
    ε_dec,
    # Redox reaction function to be evaluated during ohmic drop compensation
    redoxreaction,
    # IR compensation factor used during ohmic drop compensation
    ircompfactor,
    C_gap,
    # alternative activity coefficients
    # actcoeff! = pbi_gamma!,
)


# ╔═╡ 3489aed7-dd28-4bfa-878c-e6128a6c49f8
@test iselectroneutral(edata_unc.c_bulk, edata_unc)

# ╔═╡ 18c45fc4-f73e-4c06-a997-1c171094ba0b
R_u = Ls ./ LiquidElectrolytes.conductivity(edata_unc, edata_unc.c_bulk)

# ╔═╡ f892d2bb-4864-4ebb-9934-7a7d2c63f857
md"""
Half cell boundary condition. In order to accomodate the different types of IR compensation, we need to be careful when to call evaluate the redox reaction, and when to set the working electrode voltage.
"""

# ╔═╡ 03ee1bfe-264e-4ee9-adb2-9966f90790c1
function halfcellbc(f, u, bnode, data)
    (; Γ_we, Γ_bulk, ϕ_we, iϕ, ircompensation) = data
    bulkbcondition(f, u, bnode, data; region = Γ_bulk)
    if bnode.region == Γ_we
        if ircompensation == :none
            # With IR compensation, working electrode voltage is set by a generic
            # operator which adds the necessary compensation value to ϕ_we which is
            # internally defined by the sawtooth function
            potentialbcondition!(f, u, bnode, data, ϕ_we)
        end
        if ircompensation != :ohmicdrop
            # Ohmic drop compensation needs to evaluate the faradaic current, so the
            # reaction expression is invoked in the generic operator
            redoxreaction(f, u, bnode, data)
        end
    end
    return nothing
end

# ╔═╡ 8fd659d3-f5cc-4d3a-aa61-e6f932d5823b
function run(; ircompensation = :pseudopotentiostat, kwargs...)
    tsols = []
    cvresults = []
    grids = []
    celldatas = []
    for L in Ls

        celldata = deepcopy(edata_unc)
        celldata.ircompensation = ircompensation
        celldata.Ru = L / LiquidElectrolytes.conductivity(celldata, celldata.c_bulk)
        X = geomspace(0, L, 0.1 * ufac"nm", 0.1 * L)
        grd = simplexgrid(X)
        pnpcell = PNPSystem(grd; bcondition = halfcellbc, celldata)
        cvresult = LiquidElectrolytes.cvsweep(
            pnpcell;
            voltages = sawtooth,
            nperiods,
            store_solutions = true,
            handle_exceptions = true,
            damp_initial,
            tol_round,
            max_round,
            kwargs...
        )
        push!(cvresults, cvresult)
        push!(grids, grd)
        push!(celldatas, celldata)
        push!(tsols, cvresult.tsol)
    end
    return cvresults, tsols, grids, celldatas
end

# ╔═╡ 7a13cc37-6389-4b04-aa5b-7b2fcaa5624c
cvr_irc, tsols_irc, grids_irc, celldatas_irc = run(; ircompensation = :pseudopotentiostat);

# ╔═╡ a356f994-e64d-4430-8027-f69fe0748a86
cvr_unc, tsols_unc, grids_unc, celldatas_unc = run(; ircompensation = :none);

# ╔═╡ 170f0933-19a2-492d-83b2-3a38811b16e9
cvr_odr, tsols_odr, grids_odr, celldatas_odr = run(; ircompensation = :ohmicdrop);

# ╔═╡ b8236e74-612d-4d14-9da3-51d191e307cc
function cv2plot(voltages = :voltages)
    fig = Figure(size = (1200, 400))
    ax = Axis(fig[1, 1], title = "pseudopotentiostat, $voltages")
    xlims!(ax, -0.8, 0.8)
    for (cvr, L) in zip(cvr_irc, Ls)
        v = getfield(cvr, voltages)
        i_F = [j[1] for j in cvr.j_reaction] * edata_unc.F
        i_C = cvr.j_cap
        lines!(ax, v, i_F + i_C / ufac"mA/cm^2", label = "L=$(L / ufac"mm")mm")
    end
    axislegend(ax)

    ax = Axis(fig[1, 2], title = "no ir compensation, $voltages")
    xlims!(ax, -0.8, 0.8)
    for (cvr, L) in zip(cvr_unc, Ls)
        v = getfield(cvr, voltages)
        i_F = [j[1] for j in cvr.j_reaction] * edata_unc.F
        i_C = cvr.j_cap
        lines!(ax, v, i_F + i_C / ufac"mA/cm^2", label = "L=$(L / ufac"mm")mm")
    end
    axislegend(ax)

    ax = Axis(fig[1, 3], title = "ohmic drop, $voltages")
    xlims!(ax, -0.8, 0.8)
    for (cvr, L) in zip(cvr_odr, Ls)
        v = getfield(cvr, voltages)
        i_F = [-j[1] for j in cvr.j_we] * edata_unc.F
        i_C = cvr.j_cap
        lines!(ax, v, i_F + i_C / ufac"mA/cm^2", label = "L=$(L / ufac"mm")mm")
    end
    axislegend(ax)

    return fig
end

# ╔═╡ bb5efb53-9ba1-4e43-8519-d3e09047d8ea
cv2plot(:voltages)

# ╔═╡ 7c386967-b39f-4192-923b-c2aa0f905774
cv2plot(:dlvoltages)

# ╔═╡ 36e35875-a035-4c60-b9a7-024565129e8a
cv2plot(:sawtooth)

# ╔═╡ 4b24bf27-e2af-401a-a834-260c18924920
cmax = 0.005

# ╔═╡ 462ca6e2-7221-4cb5-8ee8-2bbdddb80d6d
figscale = 0.6

# ╔═╡ 6888ef1c-3871-4a8a-b149-29fb874ffeb3
md"""
Boundary layer thickness for ir compensated case:
"""

# ╔═╡ 619d66b8-22ce-41b3-96c3-7a04bdc6318f
md"""
Boundary layer thickness for uncompensated case:
"""

# ╔═╡ e27281a3-8a8d-415b-a877-c08a275ca5fb
"""
    blthickness(grd, celldata, tsol; species=1) -> Float64

Estimate the diffusion layer thickness from a transient solution.

Scans every time snapshot in `tsol` for the given `species` index and finds the
outermost grid node where the concentration deviates from the bulk value
`cbulk[species]` by more than 0.1 mol/m³. Returns the maximum such position over
all time steps (in SI units, i.e. meters).
"""
function blthickness(grd, celldata, tsol; species = 1)
    X = grd[Coordinates][1, :]
    xbl = 0
    for it in 1:length(tsol.t)
        u = tsol[species, :, it]
        i = findlast(c -> abs(c - celldata.c_bulk[species]) > 1.0e-1, u)
        if isnothing(i)
            i = 1
        end
        xbl = max(xbl, X[i])
    end
    return xbl
end

# ╔═╡ c521b0bf-5546-426d-aabb-5f2627ab8891
[blthickness(grid, cd, tsol) for (grid, cd, tsol) in zip(grids_irc, celldatas_irc, tsols_irc)] / ufac"nm"

# ╔═╡ 36b73fe7-48c2-40ea-a080-ded83f9751f7
[blthickness(grid, cd, tsol) for (grid, cd, tsol) in zip(grids_unc, celldatas_unc, tsols_unc)] / ufac"nm"

# ╔═╡ e666e908-fbc3-404b-8c8c-f75316fb1296
[blthickness(grid, cd, tsol) for (grid, cd, tsol) in zip(grids_odr, celldatas_odr, tsols_odr)] / ufac"nm"

# ╔═╡ 1a773c61-9af0-4f23-8eff-20ac4316182e
function plottsol(
        grd, celldata, tsol; xsplit = 5ufac"nm",
        species = celldata.iϕ,
        save = nothing, limits = nothing,
        figscale = 1
    )
    label = "C/(mol/dm^3)"

    X = grd[XCoordinates]
    wl = 150 * figscale
    wr = log10(X[end] / ufac"nm") * 60 * figscale
    figsize = (wl + wr + 200 * figscale, 300 * figscale)

    fig = Figure(size = figsize)
    axl = Axis(
        fig[1, 1],
        xlabel = L"x/nm",
        ylabel = L"t/s"
    )
    xlims!(axl, 0, 10)
    axr = Axis(
        fig[1, 2],
        yticksvisible = false,
        yticklabelsvisible = false,
        xscale = log10,
        xlabel = L"x/nm"
    )

    colgap!(fig.layout, 2)
    colsize!(fig.layout, 1, Fixed(wl * figscale))
    colsize!(fig.layout, 2, Fixed(wr * figscale))

    xlims!(axr, 10, X[end] / ufac"nm")
    #	hidedecorations!(axl)
    #	hidedecorations!(axr)
    colorscale = identity
    colormap = :terrain
    scale = 1 / ufac"mol/dm^3"
    if species == celldata.iϕ
        colorscale = identity
        colormap = :seismic
        scale = 1
    end
    T = tsol.t
    nx = length(X)
    nt = length(T)
    u = [tsol[species, ix, it] for ix in 1:nx, it in 1:nt ] * scale
    x = [X[ix] / ufac"nm" for ix in 1:nx, it in 1:nt ]
    y = [T[it] for ix in 1:nx, it in 1:nt ]
    if isnothing(limits)
        limits = extrema(u)
    end
    levels = range(limits..., length = 25)

    @time contourf!(axl, x, y, u; levels, colormap, colorscale)
    @time cr = contourf!(axr, x, y, u; levels, colormap, colorscale)
    Colorbar(fig[1, 3], cr; label)
    colsize!(fig.layout, 3, Fixed(40 * figscale))

    if !isnothing(save)
        CairoMakie.save(datadir(save), fig)
        @info "saved to $(datadir(save))"
    end
    return fig
end


# ╔═╡ 8c9af918-0d4e-411f-b520-880f591a5ed1
plots_irc = [plottsol(grd, celldata, tsol; figscale, species = 1, limits = (0, cmax)) for (grd, celldata, tsol) in zip(grids_irc, celldatas_irc, tsols_irc)];

# ╔═╡ dbb604e0-a3ae-4558-ba73-a460aaf99c02
plots_unc = [plottsol(grd, celldata, tsol; figscale, species = 1, limits = (0, cmax)) for (grd, celldata, tsol) in zip(grids_unc, celldatas_unc, tsols_unc)];

# ╔═╡ 91a78db3-e464-48b8-a266-4020b21d6f5d
plots_odr = [plottsol(grd, celldata, tsol; figscale, species = 1, limits = (0, cmax)) for (grd, celldata, tsol) in zip(grids_odr, celldatas_odr, tsols_odr)];

# ╔═╡ 4cd52496-43be-46f0-b773-cc401fe2f820
plots = hcat(plots_irc, plots_unc, plots_odr)

# ╔═╡ c78a70f3-ea86-48b6-ada8-15dc8b039933
PlutoUI.ExperimentalLayout.grid(plots)

# ╔═╡ 784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
html"""<style>.dont-panic{ display: none }</style>"""


# ╔═╡ 1a44bd29-9966-44f6-a51a-2773dc8851b8
html"""
<style>
main {
    max-width: 1500px;
    margin: 0 auto;
}
</style>
"""

# ╔═╡ Cell order:
# ╠═93a15e56-7ba2-4c77-951f-6c08ecef0d5d
# ╠═de0e3435-2736-4f75-8706-aade9016bc1f
# ╠═30698409-7cbd-4328-937a-e76f26d15151
# ╠═bf07c780-2c90-4bce-a18a-3f9a7514dfc2
# ╟─08470000-576c-478c-9528-878552cec9a5
# ╠═d956f42f-f5c4-401e-873b-d076fb27a887
# ╟─e28357fd-29b5-44e3-abb3-182f9a4ff936
# ╠═e5ab8c3c-3614-4add-aab9-689e8720d1c3
# ╟─ffb2f6e4-4a9c-4b9e-95a8-4bfb7150e502
# ╠═a5b7c828-c97d-49eb-bd7f-16462c9ebd50
# ╠═abc2f8e7-8357-47c4-bb3e-bbbbb5b816e9
# ╠═c32610f4-88e2-499d-abbc-03d6513806d3
# ╠═3489aed7-dd28-4bfa-878c-e6128a6c49f8
# ╠═18c45fc4-f73e-4c06-a997-1c171094ba0b
# ╟─f892d2bb-4864-4ebb-9934-7a7d2c63f857
# ╠═03ee1bfe-264e-4ee9-adb2-9966f90790c1
# ╠═8fd659d3-f5cc-4d3a-aa61-e6f932d5823b
# ╠═7a13cc37-6389-4b04-aa5b-7b2fcaa5624c
# ╠═a356f994-e64d-4430-8027-f69fe0748a86
# ╠═170f0933-19a2-492d-83b2-3a38811b16e9
# ╠═b8236e74-612d-4d14-9da3-51d191e307cc
# ╠═bb5efb53-9ba1-4e43-8519-d3e09047d8ea
# ╠═7c386967-b39f-4192-923b-c2aa0f905774
# ╠═36e35875-a035-4c60-b9a7-024565129e8a
# ╠═4b24bf27-e2af-401a-a834-260c18924920
# ╠═462ca6e2-7221-4cb5-8ee8-2bbdddb80d6d
# ╠═8c9af918-0d4e-411f-b520-880f591a5ed1
# ╠═dbb604e0-a3ae-4558-ba73-a460aaf99c02
# ╠═91a78db3-e464-48b8-a266-4020b21d6f5d
# ╠═4cd52496-43be-46f0-b773-cc401fe2f820
# ╠═c78a70f3-ea86-48b6-ada8-15dc8b039933
# ╟─6888ef1c-3871-4a8a-b149-29fb874ffeb3
# ╠═c521b0bf-5546-426d-aabb-5f2627ab8891
# ╟─619d66b8-22ce-41b3-96c3-7a04bdc6318f
# ╠═36b73fe7-48c2-40ea-a080-ded83f9751f7
# ╠═e666e908-fbc3-404b-8c8c-f75316fb1296
# ╠═e27281a3-8a8d-415b-a877-c08a275ca5fb
# ╠═1a773c61-9af0-4f23-8eff-20ac4316182e
# ╠═784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╠═1a44bd29-9966-44f6-a51a-2773dc8851b8
