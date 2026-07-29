### A Pluto.jl notebook ###
# v1.0.1

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
# This cell is deactivated in file, so it is not executed during CI.
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

# ╔═╡ 9a2c3973-3d57-4f50-9067-43c505a8bf97
PlutoUI.TableOfContents()

# ╔═╡ da275cc8-c35d-4b82-8c23-2f44cf6455e7
md"""
# Implementing IR compensation
"""

# ╔═╡ 349f30a6-2eb9-410a-a7c6-0deb265c03d2
md"""
## General data
"""

# ╔═╡ a7ac6c61-28a1-4699-834c-25cf84504a12
md"""
Data according to "Simulation of the cyclic voltammetric response of an outer-sphere redox species with inclusion of electrical double layer structure and ohmic potential drop" by Levey et al, PCCP 2023, [https://doi.org/10.1039/D3CP00098B](https://doi.org/10.1039/D3CP00098B).
"""

# ╔═╡ fb05cabb-026f-4f2e-a180-9a7b29e65451
begin
    scanrate = 1ufac"V/s"
    vmin = -0.5ufac"V"
    vmax = 0.1ufac"V"
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


# ╔═╡ ccd44b6c-8e15-48b9-bdfe-15ca7dc3c2e4
begin
    ε1 = 6.0 * ph"VacuumElectricPermittivity"
    ε2 = 30.0 * ph"VacuumElectricPermittivity"
    x2 = 0.29ufac"nm"
    x3 = 0.59ufac"nm"
    d1 = x2
    d2 = x3 - x2
    const C_gap_robin = 1.0 / (d1 / ε1 + d2 / ε2)
    const C_gap_dirichlet = 1.0e15
    C_gap_robin / ufac"μF/cm^2"
end

# ╔═╡ ee617289-890a-476d-b1c5-1d1267e1cda8
md"""
Dielectric decrement: this is in the moment an experiemental feature due to lack of an implementation of a better model.
"""

# ╔═╡ 147576d8-e9c4-4afd-a579-0fd285b09bd7
function ε_dec(x)
    if x[1] < 1.0e-9
        return 0.01
    else
        return 1.0
    end
end

# ╔═╡ 06833006-7050-498e-8316-612523a6697a
md"""
Electrolyte data for calculations without IR compensation
"""

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
    # Position of "pseudo reference electrode" for pseudopotentiostat approach
    xref = [10.0 * ufac"nm"],
    # Dielectric decrement
    ε_dec,
    # IR compensation
    ircompensation = NoIRCompensation()
    # alternative activity coefficients
    # actcoeff! = pbi_gamma!,
)

# ╔═╡ 1896cb6f-ccf9-4914-9827-2c0ba4bf5498
@test iselectroneutral(edata_unc.c_bulk, edata_unc)

# ╔═╡ 01b1639b-7928-4912-9e4b-e76305663ed5
md"""
Electrolyte data for calculations with IR compensation by the "pseudopotentionstat" approach.
"""

# ╔═╡ b43f5473-1624-4fb6-bbe6-c23a1322841c
edata_pts = copy(edata_unc, ircompensation = PseudoPotentiostat())

# ╔═╡ c541518f-116a-45d9-a6ad-04e2e3e1854d
md"""
Estimate of uncompensated resistance based on bulk concentrations
"""

# ╔═╡ bb3d0553-79c3-423f-9ae0-63efa5a0c25c
Ru = L / LiquidElectrolytes.conductivity(edata_unc, edata_unc.c_bulk)

# ╔═╡ 3f3e21f6-1d1f-4e8d-a378-3530f040781c
md"""
Electrolyte data for calculations with IR compensation by ohmic drop compensation via the estimate of the uncompensated resistance.
"""

# ╔═╡ d26f2c0f-d935-439a-9b30-bf47a40a7bc6
md"""
Sawtooth voltage control:
"""

# ╔═╡ 16490d65-5f1b-4428-ae90-42b21d6a48bd
sawtooth = SawTooth(; scanrate, vmin, vmax)

# ╔═╡ 9478563d-b382-420e-be99-58d525ab8240
md"""
### Boundary Conditions
"""

# ╔═╡ 82518a2b-04eb-4f64-8ab5-f46580a76bf1
md"""
Redox rection function using the rate expressions according to Landstorfer et al.:
"""

# ╔═╡ 76cd84dd-505e-4ee4-9ef7-eafaf122eaa9
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

# ╔═╡ 2fbb39ef-be14-4d2c-9958-f0f5132ff183
function redoxreaction_robin(y, u, bnode, data)
    (; iϕ, F, RT, ϕ_we, ircompensation, iϕ_we) = data
    ϕ_PET = u[iϕ]
    ϕ_L = 0
    # E is the applied electrode potential (vs. the reference electrode).
    if isactive(ircompensation)
        E = u[iϕ_we]
    else
        E = ϕ_we
    end
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

# ╔═╡ c219c417-d61a-4460-aa2f-37870e376c56
function redoxreaction(y, u, bnode, data)
    if data.C_gap < C_gap_dirichlet
        redoxreaction_robin(y, u, bnode, data)
    else
        redoxreaction_dirichlet(y, u, bnode, data)
    end
    return nothing
end

# ╔═╡ 591bfd86-fe3e-46cb-9c10-781848fb94ee
edata_odr = copy(
    edata_unc,
    ircompensation = OhmicDropEstimation(;
        factor = 0.95,
        Ru = Ru,
        species = 1,
        ne = 1,
        # Redox reaction function to be evaluated during ohmic drop compensation
        redoxreaction,
    )
)

# ╔═╡ 675ce578-8028-4af7-a236-5279a7dadc7b
md"""
Half cell boundary condition. In order to accomodate the different types of IR compensation, we need to be careful when to call evaluate the redox reaction, and when to set the working electrode voltage.
"""

# ╔═╡ 1089d336-2dce-4cb0-a6de-feb7008f30c9
function halfcellbc(f, u, bnode, data)
    (; Γ_we, Γ_bulk, ϕ_we, iϕ, ircompensation) = data
    bulkbcondition(f, u, bnode, data; region = Γ_bulk)
    if bnode.region == Γ_we
        if !isactive(ircompensation)
            potentialbcondition!(f, u, bnode, data, ϕ_we)
        end
        # With IR compensation, working electrode voltage is set by a generic
        # operator which adds the necessary compensation value to ϕ_we

        if !isa(ircompensation, OhmicDropEstimation)
            redoxreaction(f, u, bnode, data)
        end
        # Ohmic drop compensation needs to evaluate the faradaic current, so the
        # reaction expression is invoked in the generic operator

    end
    return nothing
end

# ╔═╡ 9331c262-ce81-4c1f-997a-7e078c3a88f7
md"""
## Grid and  simulation function
"""

# ╔═╡ 1043a71d-1cc8-40ad-a335-0e949ea283eb
begin
    hmin = 0.05ufac"nm"
    hmax = L / 10
    X = geomspace(0, L, hmin, hmax)
    grid = simplexgrid(X)
end

# ╔═╡ 9977fef4-3ba8-4456-a01a-3d174c434c21
function simulate(
        grid, celldata;
        damp_initial = 0.5,
        tol_round = 1.0e-9,
        max_round = 3,
        C_gap = C_gap_dirichlet,
        kwargs...
    )
    LiquidElectrolytes.trace!(10)
    xcelldata = copy(celldata; C_gap)
    pnpcell = PNPSystem(
        grid;
        bcondition = halfcellbc,
        celldata = xcelldata
    )
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
    return pnpcell, cvresult, xcelldata
end

# ╔═╡ 97670da1-940a-46e8-9388-c904618066d3
md"""
## Dirichlet boundary conditions
"""

# ╔═╡ 9d7921c9-1b19-4399-aab0-e43dbc974117
cell_unc_dirichlet, cvresult_unc_dirichlet, celldata_unc_dirichlet = simulate(grid, edata_unc, C_gap = C_gap_dirichlet)

# ╔═╡ 7c64dcde-9530-4494-bf07-dce3e9e6d247
cell_pts_dirichlet, cvresult_pts_dirichlet, celldata_pts_dirichlet = simulate(grid, edata_pts, C_gap = C_gap_dirichlet)

# ╔═╡ 2efc764c-f7a0-4afe-9bca-4bcb32232510
cell_odr_dirichlet, cvresult_odr_dirichlet, celldata_odr_dirichlet = simulate(grid, edata_odr, C_gap = C_gap_dirichlet)

# ╔═╡ 80f3ba36-30d3-4759-b92a-9fb87c63dd23
md"""
### CVs and surface concentrations 
"""

# ╔═╡ 166bb0c9-17f7-4c43-8fc6-3bd3d4afd42c
md"""
### Voltages and concentrations
"""

# ╔═╡ 023b176b-6174-4b87-b716-6435062ebbb0
md"""
t/s: $@bind t_dirichlet PlutoUI.Slider(range(cvresult_unc_dirichlet.tsol.t[begin], cvresult_unc_dirichlet.tsol.t[end], length = 1001), show_value = true, default=cvresult_unc_dirichlet.tsol.t[end]/5)
"""

# ╔═╡ 61e0320d-44d3-40ec-8f3c-58ffc1239669
md"""
### Time dependent voltages 
"""

# ╔═╡ 08ef3f36-40c4-415f-a3e6-e6d9fc574d54
md"""
### Tests
"""

# ╔═╡ baa733de-ba00-4124-9169-e2e959c7946d
md"""
## Robin boundary conditions
"""

# ╔═╡ 392b4681-5b78-42f5-ad88-d5efb73578c5
cell_unc_robin, cvresult_unc_robin, celldata_unc_robin = simulate(grid, edata_unc, C_gap = C_gap_robin)

# ╔═╡ 2ad28a94-96f5-4ce6-8b38-b6983bcba8d3
cell_pts_robin, cvresult_pts_robin, celldata_pts_robin = simulate(grid, edata_pts, C_gap = C_gap_robin)

# ╔═╡ 9ec13079-54d8-4e7b-9e21-2bb97ec57685
cell_odr_robin, cvresult_odr_robin, celldata_odr_robin = simulate(grid, edata_odr, C_gap = C_gap_robin)

# ╔═╡ b4bbe022-5df7-43f0-8f3e-325bc908c484
md"""
### CVs and surface concentrations
"""

# ╔═╡ 5f8a5f0a-8f39-4d4f-b9a2-1d0f10b9aec4
md"""
### Voltages and concentrations
"""

# ╔═╡ 481705b5-3833-4269-8f87-bef3e489d00b
md"""
t/s: $@bind t_robin PlutoUI.Slider(range(cvresult_unc_robin.tsol.t[begin], cvresult_unc_robin.tsol.t[end], length = 1001), show_value = true, default=cvresult_unc_robin.tsol.t[end]/5)
"""


# ╔═╡ ea4200b5-1443-4068-96d5-271f48da3707
md"""
### Time dependent voltages
"""

# ╔═╡ 2a800901-54bf-4d85-8133-04cf3bf7a032
md"""
### Tests
"""

# ╔═╡ edf4bfb7-62cc-4ee5-8a8a-a05416d478e7
md"""
## Appendix: plot functions
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
xplotsol(grid, cvresult_unc_dirichlet.tsol(t_dirichlet), celldata_unc_dirichlet; xcut = 1ufac"nm")
  ╠═╡ =#

# ╔═╡ 9095f5dd-4e2f-4c2d-89c3-e7aaf2c872d3
#=╠═╡
xplotsol(grid, cvresult_pts_dirichlet.tsol(t_dirichlet), celldata_pts_dirichlet; xcut = 1ufac"nm")
  ╠═╡ =#

# ╔═╡ ef27ab55-3502-4747-960e-be25240c5df8
#=╠═╡
xplotsol(grid, cvresult_odr_dirichlet.tsol(t_dirichlet), celldata_odr_dirichlet; xcut = 1ufac"nm")
  ╠═╡ =#

# ╔═╡ 9c7ce60a-3abe-44c2-ad6f-4784b5a83a3a
#=╠═╡
xplotsol(grid, cvresult_unc_robin.tsol(t_robin), celldata_unc_robin; xcut = 1ufac"nm")
  ╠═╡ =#

# ╔═╡ 20927491-e3f9-4e18-a92f-1abc803e55c3
#=╠═╡
xplotsol(grid, cvresult_pts_robin.tsol(t_robin), celldata_pts_robin; xcut = 1ufac"nm")
  ╠═╡ =#

# ╔═╡ 9296e3fa-95b3-40e2-b3ff-08a2e2cb7cbb
#=╠═╡
xplotsol(grid, cvresult_odr_robin.tsol(t_robin), celldata_odr_robin; xcut = 1ufac"nm")
  ╠═╡ =#

# ╔═╡ 784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
html"""<style>.dont-panic{ display: none }</style>"""


# ╔═╡ 25ab8427-22f1-4b68-b4a8-bae98af47fe9
#=╠═╡
function plotcv(cvresult, celldata)
    !isdefined(Main, :PlutoRunner) && return
    fig = Figure(size = (700, 300))
    ax = Axis(fig[1, 1], title="CV")
    j_we = -hcat(cvresult.j_we...)[1, :] * celldata.F
    lines!(ax, cvresult.voltages, j_we / ufac"mA/cm^2", label=L"I_F")
    lines!(ax, cvresult.voltages, cvresult.j_cap / ufac"mA/cm^2", label=L"I_C")
    lines!(ax, cvresult.voltages, (cvresult.j_cap + j_we) / ufac"mA/cm^2",label=L"I_F+I_C")

	axislegend(ax, position = :lt)

    c_O = [u[iO, :] |> maximum for u in cvresult.tsol.u[2:end]] / ufac"mol/dm^3"
    c_R = [u[iR, :] |> maximum for u in cvresult.tsol.u[2:end]] / ufac"mol/dm^3"

    ax2 = Axis(fig[1, 2], xlabel = "E/V", ylabel = L"c/(mol\cdot dm^{-3})",
			  title= "Max concentration in DL")

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
plotcv(cvresult_unc_dirichlet, celldata_unc_dirichlet)
  ╠═╡ =#

# ╔═╡ 40fc09c4-8820-4208-b031-89d173d95f1a
#=╠═╡
plotcv(cvresult_pts_dirichlet, celldata_pts_dirichlet)
  ╠═╡ =#

# ╔═╡ 3ef89245-a0c8-4129-bfb2-77c6cd1df8ef
#=╠═╡
plotcv(cvresult_odr_dirichlet, celldata_odr_dirichlet)
  ╠═╡ =#

# ╔═╡ 5082c976-51e2-4d11-aa89-de7e26fbd602
#=╠═╡
plotcv(cvresult_unc_robin, celldata_unc_robin)
  ╠═╡ =#

# ╔═╡ cec0dd3f-6dbc-4c93-ae90-b02753a7d33b
#=╠═╡
plotcv(cvresult_pts_robin, celldata_pts_robin)
  ╠═╡ =#

# ╔═╡ 3a9ce9cf-a94c-4837-bbb5-585f92de7551
#=╠═╡
plotcv(cvresult_odr_robin, celldata_odr_robin)
  ╠═╡ =#

# ╔═╡ 2cba58da-4576-4bec-af73-af8c3581a58e
#=╠═╡
function plotresult2(cvresult, celldata)
    !isdefined(Main, :PlutoRunner) && return
    times = cvresult.tsol.t
    (; iref, iϕ) = celldata
    j_we = hcat(cvresult.j_reaction...)[1, :] * celldata.F
    Ωdrop = (j_we + cvresult.j_cap) * Ru
    st = cvresult.sawtooth
    fig = Figure(size = (700, 400))
    Label(fig[0, 1:2], "IR compensation: $(celldata.ircompensation), scanrate=$(scanrate) V/s")
    ax1 = Axis(fig[1, 1], xlabel = "t/s", ylabel = "Δϕ/V")
    ax2 = Axis(fig[2, 1], xlabel = "t/s", ylabel = "Δϕ/V")

    lines!(ax1, cvresult.times, cvresult.voltages, label = L"ϕ_0-ϕ_L")
    lines!(ax1, cvresult.times, cvresult.dlvoltages, label = L"ϕ_0-ϕ_{DL}")
    Legend(fig[1, 2], ax1)

    lines!(ax2, cvresult.times, st - cvresult.dlvoltages, label = L"ϕ_{sawtooth} - ϕ_{DL} ")
    scatter!(ax2, cvresult.times, Ωdrop, label = L"R_u\cdot (I_F+I_C)", markersize = 3, color = :red)
    Legend(fig[2, 2], ax2)
    return fig
end
  ╠═╡ =#

# ╔═╡ ce8a8eef-04d9-4c0f-8ce3-e8248dc9c623
#=╠═╡
plotresult2(cvresult_unc_dirichlet, celldata_unc_dirichlet)
  ╠═╡ =#

# ╔═╡ 62d58a40-ff6c-4277-9622-fbf02dd3515e
#=╠═╡
plotresult2(cvresult_pts_dirichlet, celldata_pts_dirichlet)
  ╠═╡ =#

# ╔═╡ d86e93aa-5d05-43e1-8934-d88e879d1bc8
#=╠═╡
plotresult2(cvresult_odr_dirichlet, celldata_odr_dirichlet)
  ╠═╡ =#

# ╔═╡ c58759ee-2a2f-42bc-ac66-4f84b01ffbab
#=╠═╡
plotresult2(cvresult_unc_robin, celldata_unc_robin)
  ╠═╡ =#

# ╔═╡ 207144e2-ee6b-441d-82f7-dbbe699dcc27
#=╠═╡
plotresult2(cvresult_pts_robin, celldata_pts_robin)
  ╠═╡ =#

# ╔═╡ 95d0d8fa-984b-4588-b3ab-97c64b697343
#=╠═╡
plotresult2(cvresult_odr_robin, celldata_odr_robin)
  ╠═╡ =#

# ╔═╡ 962d64d5-13cf-4bc7-b369-65853c74f395
function checkir(cvresult, celldata; tol = 1.0e-11)
    times = cvresult.tsol.t
    (; iref, iϕ) = celldata
    return norm(cvresult.dlvoltages - cvresult.sawtooth, Inf) < tol
end

# ╔═╡ fc406613-fa88-4b18-8b26-53afae756cbc
@test checkir(cvresult_pts_dirichlet, celldata_pts_dirichlet)

# ╔═╡ fdf4e114-a627-4941-979a-a168a63f4ff0
@test checkir(cvresult_odr_dirichlet, celldata_odr_dirichlet; tol = 0.03)

# ╔═╡ c8a176a6-030d-4a8c-86bb-2c8ed5e64612
@test checkir(cvresult_pts_robin, celldata_pts_robin)

# ╔═╡ e30667e7-7ef5-40af-851a-e189cd78fd6b
@test checkir(cvresult_odr_robin, celldata_odr_robin, tol = 0.03)

# ╔═╡ Cell order:
# ╠═670585dc-bee1-4d2f-88b7-282a791badc8
# ╠═93a15e56-7ba2-4c77-951f-6c08ecef0d5d
# ╠═9a2c3973-3d57-4f50-9067-43c505a8bf97
# ╟─da275cc8-c35d-4b82-8c23-2f44cf6455e7
# ╟─349f30a6-2eb9-410a-a7c6-0deb265c03d2
# ╟─a7ac6c61-28a1-4699-834c-25cf84504a12
# ╠═fb05cabb-026f-4f2e-a180-9a7b29e65451
# ╠═502e3108-72cd-40c0-b0f1-6fa922446ef9
# ╠═ccd44b6c-8e15-48b9-bdfe-15ca7dc3c2e4
# ╟─ee617289-890a-476d-b1c5-1d1267e1cda8
# ╠═147576d8-e9c4-4afd-a579-0fd285b09bd7
# ╟─06833006-7050-498e-8316-612523a6697a
# ╠═55ce733f-82ef-4c66-94e6-b948e9865385
# ╠═1896cb6f-ccf9-4914-9827-2c0ba4bf5498
# ╟─01b1639b-7928-4912-9e4b-e76305663ed5
# ╠═b43f5473-1624-4fb6-bbe6-c23a1322841c
# ╟─c541518f-116a-45d9-a6ad-04e2e3e1854d
# ╠═bb3d0553-79c3-423f-9ae0-63efa5a0c25c
# ╟─3f3e21f6-1d1f-4e8d-a378-3530f040781c
# ╠═591bfd86-fe3e-46cb-9c10-781848fb94ee
# ╟─d26f2c0f-d935-439a-9b30-bf47a40a7bc6
# ╠═16490d65-5f1b-4428-ae90-42b21d6a48bd
# ╟─9478563d-b382-420e-be99-58d525ab8240
# ╟─82518a2b-04eb-4f64-8ab5-f46580a76bf1
# ╠═76cd84dd-505e-4ee4-9ef7-eafaf122eaa9
# ╠═2fbb39ef-be14-4d2c-9958-f0f5132ff183
# ╠═c219c417-d61a-4460-aa2f-37870e376c56
# ╟─675ce578-8028-4af7-a236-5279a7dadc7b
# ╠═1089d336-2dce-4cb0-a6de-feb7008f30c9
# ╟─9331c262-ce81-4c1f-997a-7e078c3a88f7
# ╠═1043a71d-1cc8-40ad-a335-0e949ea283eb
# ╠═9977fef4-3ba8-4456-a01a-3d174c434c21
# ╟─97670da1-940a-46e8-9388-c904618066d3
# ╠═9d7921c9-1b19-4399-aab0-e43dbc974117
# ╠═7c64dcde-9530-4494-bf07-dce3e9e6d247
# ╠═2efc764c-f7a0-4afe-9bca-4bcb32232510
# ╟─80f3ba36-30d3-4759-b92a-9fb87c63dd23
# ╠═5b8c96f4-9d64-4475-b04d-5950002bafaa
# ╠═40fc09c4-8820-4208-b031-89d173d95f1a
# ╠═3ef89245-a0c8-4129-bfb2-77c6cd1df8ef
# ╟─166bb0c9-17f7-4c43-8fc6-3bd3d4afd42c
# ╟─023b176b-6174-4b87-b716-6435062ebbb0
# ╠═f68da385-91a5-4878-8e71-11d380fa7d9b
# ╠═9095f5dd-4e2f-4c2d-89c3-e7aaf2c872d3
# ╠═ef27ab55-3502-4747-960e-be25240c5df8
# ╟─61e0320d-44d3-40ec-8f3c-58ffc1239669
# ╠═ce8a8eef-04d9-4c0f-8ce3-e8248dc9c623
# ╠═62d58a40-ff6c-4277-9622-fbf02dd3515e
# ╠═d86e93aa-5d05-43e1-8934-d88e879d1bc8
# ╟─08ef3f36-40c4-415f-a3e6-e6d9fc574d54
# ╠═fc406613-fa88-4b18-8b26-53afae756cbc
# ╠═fdf4e114-a627-4941-979a-a168a63f4ff0
# ╟─baa733de-ba00-4124-9169-e2e959c7946d
# ╠═392b4681-5b78-42f5-ad88-d5efb73578c5
# ╠═2ad28a94-96f5-4ce6-8b38-b6983bcba8d3
# ╠═9ec13079-54d8-4e7b-9e21-2bb97ec57685
# ╟─b4bbe022-5df7-43f0-8f3e-325bc908c484
# ╠═5082c976-51e2-4d11-aa89-de7e26fbd602
# ╠═cec0dd3f-6dbc-4c93-ae90-b02753a7d33b
# ╠═3a9ce9cf-a94c-4837-bbb5-585f92de7551
# ╟─5f8a5f0a-8f39-4d4f-b9a2-1d0f10b9aec4
# ╟─481705b5-3833-4269-8f87-bef3e489d00b
# ╠═9c7ce60a-3abe-44c2-ad6f-4784b5a83a3a
# ╠═20927491-e3f9-4e18-a92f-1abc803e55c3
# ╠═9296e3fa-95b3-40e2-b3ff-08a2e2cb7cbb
# ╟─ea4200b5-1443-4068-96d5-271f48da3707
# ╠═c58759ee-2a2f-42bc-ac66-4f84b01ffbab
# ╠═207144e2-ee6b-441d-82f7-dbbe699dcc27
# ╠═95d0d8fa-984b-4588-b3ab-97c64b697343
# ╟─2a800901-54bf-4d85-8133-04cf3bf7a032
# ╠═c8a176a6-030d-4a8c-86bb-2c8ed5e64612
# ╠═e30667e7-7ef5-40af-851a-e189cd78fd6b
# ╟─edf4bfb7-62cc-4ee5-8a8a-a05416d478e7
# ╠═28c3bed8-f225-4405-bbdd-90951ef6211f
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╠═25ab8427-22f1-4b68-b4a8-bae98af47fe9
# ╠═2cba58da-4576-4bec-af73-af8c3581a58e
# ╠═962d64d5-13cf-4bc7-b369-65853c74f395
