### A Pluto.jl notebook ###
# v0.20.25

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 670585dc-bee1-4d2f-88b7-282a791badc8
# ╠═╡ skip_as_script = true
#=╠═╡
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "docs"))
    using Revise
    using CairoMakie
    using GridVisualize
    default_plotter!(CairoMakie)
    using Colors
    set_theme!(theme_dark())

end
  ╠═╡ =#

# ╔═╡ 93a15e56-7ba2-4c77-951f-6c08ecef0d5d
begin
    using Test
    using LiquidElectrolytes
    using ExtendableGrids
    using LessUnitful
    using PlutoUI
    using VoronoiFVM
    using LinearAlgebra
end

# ╔═╡ fb05cabb-026f-4f2e-a180-9a7b29e65451
begin
    scanrate = 100ufac"V/s"
    L = 1ufac"mm"
    nperiods = 2
    C_elec = 0.01ufac"mol/dm^3"
    C_O = 1.0e-3 * ufac"mol / dm^3"
    C_R = 1.0e-10 * ufac"mol / dm^3"
    ircompfactor = 0.95
    z = [3, 2, 1, -1]
    f_κ = 2

end;

# ╔═╡ 502e3108-72cd-40c0-b0f1-6fa922446ef9
begin
    const k0::Float64 = 13.5 * ufac"cm / s"
    const α = 0.45
    const E0::Float64 = -0.173ufac"V"
    const k0::Float64 = 13.5 * ufac"cm / s"

    const iO = 1
    const iR = 2
end;


# ╔═╡ 147576d8-e9c4-4afd-a579-0fd285b09bd7
function ε_dec(x)
    if x[1] < 1.0e-9
        return 0.01
    else
        return 1.0
    end
end

# ╔═╡ 6ef4e7d7-f268-4547-ad46-ae01bed4e800
function pbi_gamma!(γ, c, p, electrolyte)
    sumc = zero(eltype(c))
    for i in 1:electrolyte.nc
        sumc += c[i]
    end
    g = 1.0 / (1.0 - sumc * electrolyte.v0)
    for i in 1:electrolyte.nc
        γ[i] = g
    end
    return nothing
end


# ╔═╡ 96599aed-df1f-4776-8c91-e8c65110d840
function pbo_gamma!(γ, c, p, electrolyte)
    (; cspecies) = electrolyte
    for ic in cspecies
        γ[ic] = 1
    end
    return nothing
end

# ╔═╡ 1043a71d-1cc8-40ad-a335-0e949ea283eb
begin
    hmin = 0.01ufac"nm"
    hmax = L / 50
    X = geomspace(0, L, hmin, hmax)
    grid = simplexgrid(X)
end

# ╔═╡ 16490d65-5f1b-4428-ae90-42b21d6a48bd
sawtooth = SawTooth(
    scanrate = scanrate,
    vmin = -0.5, vmax = 0.1
)


# ╔═╡ 76cd84dd-505e-4ee4-9ef7-eafaf122eaa9
function redoxreaction(f, u, bnode, data)
    (; ip, iϕ, v, v0, F, RT, κ) = data
    c0, barc = c0_barc(u, data)
    μR = chemical_potential(u[iR], barc, u[ip], v[iR] + κ[iR] * v0, data)
    μO = chemical_potential(u[iO], barc, u[ip], v[iO] + κ[iO] * v0, data)
    A = (μR - μO - F * E0) / RT
    j_O = rrate(k0, α, A)
    f[iO] -= j_O
    return f[iR] += j_O
end

# ╔═╡ 55ce733f-82ef-4c66-94e6-b948e9865385
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

    x_ref = [10.0 * ufac"nm"],
    ε_dec,
    redoxreaction,
    #actcoeff! = pbi_gamma!,
    #		rlog=RLog(),
    #	rexp=RExp()

)

# ╔═╡ 1896cb6f-ccf9-4914-9827-2c0ba4bf5498
iselectroneutral(edata_unc.c_bulk, edata_unc)

# ╔═╡ bb3d0553-79c3-423f-9ae0-63efa5a0c25c
R_u = L / LiquidElectrolytes.conductivity(edata_unc, edata_unc.c_bulk)

# ╔═╡ 1089d336-2dce-4cb0-a6de-feb7008f30c9
function halfcellbc(f, u, bnode, data)
    (; Γ_we, Γ_bulk, ϕ_we, iϕ, ircompensation) = data
    bulkbcondition(f, u, bnode, data; region = Γ_bulk)
    if bnode.region == Γ_we
        if !isa(bnode, NamedTuple) && ircompensation == :none
            boundary_dirichlet!(f, u, bnode; region = Γ_we, species = iϕ, value = ϕ_we)
        end
        if ircompensation != :ohmicdrop
            redoxreaction(f, u, bnode, data)
        end
    end
    return nothing
end

# ╔═╡ 9977fef4-3ba8-4456-a01a-3d174c434c21
function simulate(
        grid, celldata;
        damp_initial = 0.5,
        tol_round = 1.0e-9,
        max_round = 3,
        kwargs...

    )

    pnpcell = PNPSystem(grid; bcondition = halfcellbc, celldata)
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
    return pnpcell, cvresult
end

# ╔═╡ 9d7921c9-1b19-4399-aab0-e43dbc974117
cell_unc, cvresult_unc = simulate(grid, edata_unc)

# ╔═╡ b43f5473-1624-4fb6-bbe6-c23a1322841c
begin
    edata_pts = deepcopy(edata_unc)
    edata_pts.ircompensation = :pseudopotentiostat
    edata_odr = deepcopy(edata_unc)
    edata_odr.ircompensation = :ohmicdrop
    edata_odr.Ru = R_u
end

# ╔═╡ 7c64dcde-9530-4494-bf07-dce3e9e6d247
cell_pts, cvresult_pts = simulate(grid, edata_pts)

# ╔═╡ 2efc764c-f7a0-4afe-9bca-4bcb32232510
cell_odr, cvresult_odr = simulate(grid, edata_odr)

# ╔═╡ 023b176b-6174-4b87-b716-6435062ebbb0
md"""
t/s: $@bind t PlutoUI.Slider(range(cvresult_unc.tsol.t[begin], cvresult_unc.tsol.t[end], length = 1001), show_value = true)
"""

# ╔═╡ 28c3bed8-f225-4405-bbdd-90951ef6211f
#=╠═╡
function xplotsol(grid, sol, celldata; xcut = 10ufac"nm")
    !isdefined(Main, :PlutoRunner) && return

    X = grid[XCoordinates]
    (; iϕ) = celldata

    function pp(axl, axr)
        ylims!(axl, -0.5, 0.5)
        ylims!(axr, 1.0e-5, 10)
        lines!(
            axr, X / ufac"nm", sol[iO, :] / (ufac"mol / dm^3"),
            color = :red, label = "O"
        )
        lines!(
            axr, X / ufac"nm", sol[iR, :] / (ufac"mol / dm^3"),
            color = :blue, label = "R"
        )
        lines!(axl, X / ufac"nm", sol[iϕ, :], color = :green, label = "ϕ")
        return nothing
    end

    fig = Figure(size = (700, 300))
    Label(
        fig[0, 1:2],
        "ircomp=$(celldata.ircompensation) scanrate=$(sawtooth.scanrate)V/s, E=$(round(sol[iϕ, 1], sigdigits = 3))V, "
    )

    colgap!(fig.layout, 2)

    axl = Axis(fig[1, 1]; xscale = identity, xlabel = "x/nm", ylabel = "U/V")
    axr = Axis(
        fig[1, 1], yaxisposition = :right; xscale = identity,
        yticksvisible = false,
        yticklabelsvisible = false,
        yscale = log10
    )
    xlims!(axl, 0, xcut / ufac"nm")
    xlims!(axr, 0, xcut / ufac"nm")
    linkxaxes!(axl, axr)
    pp(axl, axr)
    axislegend(axl, position = :lt, backgroundcolor = RGBA(1, 1, 1, 0.5))

    axl = Axis(
        fig[1, 2]; xscale = log10, xlabel = "x/nm",
        yticksvisible = false,
        yticklabelsvisible = false,
    )
    axr = Axis(
        fig[1, 2], yaxisposition = :right; xscale = log10, ylabel = "c/(mol/dm^3)",
        yscale = log10
    )
    xlims!(axl, xcut / ufac"nm", X[end] / ufac"nm")
    xlims!(axr, xcut / ufac"nm", X[end] / ufac"nm")
    linkxaxes!(axl, axr)
    pp(axl, axr)
    axislegend(axr, position = :rt, backgroundcolor = RGBA(1, 1, 1, 0.5))
    return fig
end

  ╠═╡ =#

# ╔═╡ f68da385-91a5-4878-8e71-11d380fa7d9b
#=╠═╡
xplotsol(grid, cvresult_unc.tsol(t), edata_unc; xcut = 1ufac"nm")
  ╠═╡ =#

# ╔═╡ 9095f5dd-4e2f-4c2d-89c3-e7aaf2c872d3
#=╠═╡
xplotsol(grid, cvresult_pts.tsol(t), edata_pts; xcut = 1ufac"nm")
  ╠═╡ =#

# ╔═╡ ef27ab55-3502-4747-960e-be25240c5df8
#=╠═╡
xplotsol(grid, cvresult_odr.tsol(t), edata_odr; xcut = 1ufac"nm")
  ╠═╡ =#

# ╔═╡ 784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
html"""<style>.dont-panic{ display: none }</style>"""


# ╔═╡ 25ab8427-22f1-4b68-b4a8-bae98af47fe9
#=╠═╡
function plotcv(cvresult, celldata)
    !isdefined(Main, :PlutoRunner) && return
    fig = Figure(size = (700, 300))
    ax = Axis(fig[1, 1])
    j_we = -hcat(cvresult.j_we...)[1, :] * celldata.F
    lines!(ax, cvresult.voltages, j_we / ufac"mA/cm^2")
    #lines!(ax,cvresult.voltages, -jx*edata.F/ufac"mA/cm^2")
    lines!(ax, cvresult.voltages, cvresult.j_cap / ufac"mA/cm^2")
    lines!(ax, cvresult.voltages, (cvresult.j_cap + j_we) / ufac"mA/cm^2")
    lines!(ax, cvresult.voltages[1:end], cvresult.j_dsp / ufac"mA/cm^2")

    c_O = [u[iO, :] |> maximum for u in cvresult.tsol.u[2:end]] / ufac"mol/dm^3"
    c_R = [u[iR, :] |> maximum for u in cvresult.tsol.u[2:end]] / ufac"mol/dm^3"

    ax2 = Axis(fig[1, 2], xlabel = "E/V", ylabel = L"c/(mol\cdot dm^{-3})")

    lines!(ax2, cvresult.voltages, c_O, color = :red, label = L"c_O(x_3)")
    lines!(
        ax2, cvresult.voltages, c_R,
        color = :blue, label = L"c_R(x_3)"
    )
    axislegend(ax2, position = :rt)
    Label(
        fig[0, 1:2],
        "ircomp=$(celldata.ircompensation) scanrate=$(sawtooth.scanrate)V/s"
    )

    return fig
end
  ╠═╡ =#

# ╔═╡ 5b8c96f4-9d64-4475-b04d-5950002bafaa
#=╠═╡
plotcv(cvresult_unc, edata_unc)
  ╠═╡ =#

# ╔═╡ 40fc09c4-8820-4208-b031-89d173d95f1a
#=╠═╡
plotcv(cvresult_pts, edata_pts)
  ╠═╡ =#

# ╔═╡ 3ef89245-a0c8-4129-bfb2-77c6cd1df8ef
#=╠═╡
plotcv(cvresult_odr, edata_odr)
  ╠═╡ =#

# ╔═╡ 2cba58da-4576-4bec-af73-af8c3581a58e
#=╠═╡
function plotresult2(cvresult, celldata)
    !isdefined(Main, :PlutoRunner) && return
    times = cvresult.tsol.t
    (; i_ref, iϕ) = celldata
    j_we = -hcat(cvresult.j_we...)[1, :] * celldata.F
    Ωdrop = (j_we + cvresult.j_cap) * R_u
    dlvolts = [ u[iϕ, 1] - u[iϕ, i_ref]    for u in cvresult.tsol.u]
    st = [sawtooth(t) for t in times]
    fig = Figure(size = (700, 400))
    Label(fig[0, 1:2], "IR compensation: $(celldata.ircompensation), scanrate=$(scanrate) V/s")
    ax1 = Axis(fig[1, 1], xlabel = "t/s", ylabel = "Δϕ/V")
    ax2 = Axis(fig[2, 1], xlabel = "t/s", ylabel = "Δϕ/V")

    lines!(ax1, times[2:end], cvresult.voltages, label = L"ϕ_0-ϕ_L")
    lines!(ax1, times, dlvolts, label = L"ϕ_0-ϕ_{DL}")
    Legend(fig[1, 2], ax1)

    lines!(
        ax2, times, st - dlvolts,
        label = L"ϕ_{sawtooth} - ϕ_{DL} "
    )
    scatter!(ax2, cvresult.times, Ωdrop, label = L"R_u\cdot (I_F+I_C)", markersize = 3, color = :red)
    Legend(fig[2, 2], ax2)
    return fig
end
  ╠═╡ =#

# ╔═╡ ce8a8eef-04d9-4c0f-8ce3-e8248dc9c623
#=╠═╡
plotresult2(cvresult_unc, edata_unc)
  ╠═╡ =#

# ╔═╡ 62d58a40-ff6c-4277-9622-fbf02dd3515e
#=╠═╡
plotresult2(cvresult_pts, edata_pts)
  ╠═╡ =#

# ╔═╡ d86e93aa-5d05-43e1-8934-d88e879d1bc8
#=╠═╡
plotresult2(cvresult_odr, edata_odr)
  ╠═╡ =#

# ╔═╡ 962d64d5-13cf-4bc7-b369-65853c74f395
function checkir(cvresult, celldata; tol = 1.0e-11)
    times = cvresult.tsol.t
    (; i_ref, iϕ) = celldata
    dlvolts = [ u[iϕ, 1] - u[iϕ, i_ref]    for u in cvresult.tsol.u]
    st = [sawtooth(t) for t in times]
    return norm(dlvolts - st, Inf) < tol
end

# ╔═╡ fc406613-fa88-4b18-8b26-53afae756cbc
@test checkir(cvresult_pts, edata_pts)

# ╔═╡ fdf4e114-a627-4941-979a-a168a63f4ff0
@test checkir(cvresult_odr, edata_odr; tol = 0.03)

# ╔═╡ Cell order:
# ╠═670585dc-bee1-4d2f-88b7-282a791badc8
# ╠═93a15e56-7ba2-4c77-951f-6c08ecef0d5d
# ╠═fb05cabb-026f-4f2e-a180-9a7b29e65451
# ╠═502e3108-72cd-40c0-b0f1-6fa922446ef9
# ╠═147576d8-e9c4-4afd-a579-0fd285b09bd7
# ╠═6ef4e7d7-f268-4547-ad46-ae01bed4e800
# ╠═96599aed-df1f-4776-8c91-e8c65110d840
# ╠═55ce733f-82ef-4c66-94e6-b948e9865385
# ╠═1896cb6f-ccf9-4914-9827-2c0ba4bf5498
# ╠═bb3d0553-79c3-423f-9ae0-63efa5a0c25c
# ╠═1043a71d-1cc8-40ad-a335-0e949ea283eb
# ╠═16490d65-5f1b-4428-ae90-42b21d6a48bd
# ╠═76cd84dd-505e-4ee4-9ef7-eafaf122eaa9
# ╠═1089d336-2dce-4cb0-a6de-feb7008f30c9
# ╠═9977fef4-3ba8-4456-a01a-3d174c434c21
# ╠═9d7921c9-1b19-4399-aab0-e43dbc974117
# ╠═b43f5473-1624-4fb6-bbe6-c23a1322841c
# ╠═7c64dcde-9530-4494-bf07-dce3e9e6d247
# ╠═2efc764c-f7a0-4afe-9bca-4bcb32232510
# ╠═5b8c96f4-9d64-4475-b04d-5950002bafaa
# ╠═40fc09c4-8820-4208-b031-89d173d95f1a
# ╠═3ef89245-a0c8-4129-bfb2-77c6cd1df8ef
# ╟─023b176b-6174-4b87-b716-6435062ebbb0
# ╟─f68da385-91a5-4878-8e71-11d380fa7d9b
# ╠═9095f5dd-4e2f-4c2d-89c3-e7aaf2c872d3
# ╠═ef27ab55-3502-4747-960e-be25240c5df8
# ╠═ce8a8eef-04d9-4c0f-8ce3-e8248dc9c623
# ╠═62d58a40-ff6c-4277-9622-fbf02dd3515e
# ╠═d86e93aa-5d05-43e1-8934-d88e879d1bc8
# ╠═28c3bed8-f225-4405-bbdd-90951ef6211f
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╠═25ab8427-22f1-4b68-b4a8-bae98af47fe9
# ╠═2cba58da-4576-4bec-af73-af8c3581a58e
# ╠═962d64d5-13cf-4bc7-b369-65853c74f395
# ╠═fc406613-fa88-4b18-8b26-53afae756cbc
# ╠═fdf4e114-a627-4941-979a-a168a63f4ff0
