### A Pluto.jl notebook ###
# v0.20.3

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

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin
    import Pkg # hide
    Pkg.activate(joinpath(@__DIR__, "..", "docs")) # hide
    using Revise

    using PlutoUI, HypertextLiteral
    using LiquidElectrolytes
    using LessUnitful
    using ExtendableGrids
    using VoronoiFVM
    using GridVisualize
    using StaticArrays
    using Interpolations
    if isdefined(Main, :PlutoRunner)
        using CairoMakie
        CairoMakie.activate!(type = "svg")
        default_plotter!(CairoMakie)
    end
end

# ╔═╡ 78685b0c-b27a-450a-9bae-5948c578d688
pkgdir(LiquidElectrolytes)

# ╔═╡ 84e05551-7d51-4b2c-88f2-b186ad6a244a
md"""
## Setup
"""

# ╔═╡ 16b9af50-11c5-4bdf-b5d8-8d2d9331b5e9
md"""
### Units
"""

# ╔═╡ 515baffb-2e22-401d-aacd-15971dd4365e
begin
    LessUnitful.@phconstants R N_A e
    @unitfactors nm cm μF mol dm s V
    const F = N_A * e
end;

# ╔═╡ bc896b9c-03a9-4da0-be80-f096766cb039
md"""
### Data
"""

# ╔═╡ 50ff0114-f9b1-4bd0-8190-153fef07ef3a
begin
    const vmin = -1 * V
    const vmax = 1 * V
    const vdelta = 0.1 * V
    const molarity = 0.1
    const nref = 0
    const κ = 10.0
    const scheme = :μex
    const R0 = 4.0e-15mol / (cm^2 * s)
    const epsreg = 1.0e-20


    const Δg = 0.0
    const β = 0.5
    const ϕ_we = 0.0
    const ihplus = 1
    const iso4 = 2
    const io2 = 3
    const z = [1, -2, 0]
    const allκ = [κ, κ, 0]

end

# ╔═╡ 231ebfc2-51b6-4027-8be2-75891dfed7e0
md"""
### Control
"""

# ╔═╡ a4e0379d-f7f6-4b61-bf38-5eb17f67505a
solver_control =
    (max_round = 3, tol_round = 1.0e-10, reltol = 1.0e-7, tol_mono = 1.0e-7, verbose = "e")

# ╔═╡ 970389b5-d2c1-4992-9978-aca1ccd3d2fc
md"""
### Reaction description
"""

# ╔═╡ 5482615c-cb5b-4c44-99c4-19416eabae7f
md"""
```math
  4H^+ + O_2 + 4e^- \leftrightharpoons  2H_2O
```
"""

# ╔═╡ 136cf1aa-75c6-4aaa-b93b-801391ec800c
md"""
In the following reaction function, the  balance with the solvent is fulfilled automatically via the incompressibility constraint. Any material which is removed via the boundary reaction is automatically replaced by solvent with the corresponding volume. So "ignoring the solvent" here is correct.
"""

# ╔═╡ b916eb92-8ba8-49aa-bd7c-1bfc91e813d4
function halfcellbc(f, u, bnode, data)
    bulkbcondition(f, u, bnode, data)
    (; iϕ, eneutral, ϕ_we, Γ_we, RT) = data

    return if bnode.region == Γ_we
        f .= 0.0
        if !data.eneutral
            boundary_dirichlet!(
                f,
                u,
                bnode;
                species = iϕ,
                region = data.Γ_we,
                value = data.ϕ_we,
            )
        end
        μh2o, μ = chemical_potentials!(MVector{4, eltype(u)}(undef), u, data)
        A =
            (4 * μ[ihplus] + μ[io2] - 2μh2o + Δg + eneutral * F * (u[iϕ] - data.ϕ_we)) /
            (RT)
        r = rrate(R0, β, A)
        f[ihplus] -= 4 * r
        f[io2] -= r
    end
end


# ╔═╡ 392a648c-12a2-41fb-b3f6-aa3cfe3cbcd7
grid = let
    hmin = 1.0e-1 * nm * 2.0^(-nref)
    hmax = 1.0 * nm * 2.0^(-nref)
    L = 20.0 * nm
    X = geomspace(0, L, hmin, hmax)
    simplexgrid(X)
end

# ╔═╡ 420a82e0-5fc2-47ea-8916-d88910655d50
md"""
## Nernst-Planck halfcell
"""

# ╔═╡ 60d85155-9fa7-4740-9769-212ceef1918b
begin
    celldata = ElectrolyteData(;
        nc = 3,
        z,
        κ = allκ,
        Γ_we = 1,
        Γ_bulk = 2,
        eneutral = false,
        scheme,
        epsreg,
    )

    c_bulk = celldata.c_bulk
    c_bulk[io2] = 0.001 * mol / dm^3
    c_bulk[iso4] = molarity * mol / dm^3
    c_bulk[ihplus] = 2.0 * molarity * mol / dm^3


end

# ╔═╡ bd3b54ce-de22-4f8b-aa27-7f56c5dca140
let
    c0, barc = c0_barc(c_bulk, celldata)
    @info c0, barc
    v0 = celldata.v0
    v = celldata.v
    κx = celldata.κ
    cv = c0 * v0
    for i in 1:celldata.nc
        cv += c_bulk[i] * (v[i] + κx[i] * v0)
    end
    cv
end

# ╔═╡ 8683cf81-4b53-4099-a668-37c3641802dd
@assert isapprox(celldata.c_bulk' * celldata.z, 0, atol = 1.0e-12)

# ╔═╡ 1d23b6b5-88cf-4200-b612-82dffcc5cca7
cell = PNPSystem(grid; bcondition = halfcellbc, celldata)

# ╔═╡ 763c393c-e0c8-447e-92e4-f2a5f0de2a30
begin

    result = LiquidElectrolytes.ivsweep(
        cell;
        store_solutions = true,
        voltages = vmin:vdelta:vmax,
        solver_control...,
    )
    tsol = voltages_solutions(result)
    for it in 1:length(tsol.t)
        tsol.u[it][io2, :] /= mol / dm^3
        tsol.u[it][ihplus, :] /= mol / dm^3
        tsol.u[it][iso4, :] /= mol / dm^3
    end
    volts = result.voltages
end

# ╔═╡ 2717060b-9e75-439f-8aa8-614f0930a163
md"""
### Discussion
"""

# ╔═╡ 5bc4f11f-24c6-4af8-a554-1b5771f1f2b0
function currents_h2o(result; kwargs...)
    return -4κ * currents(result, ihplus; kwargs...) - (κ + 2) * currents(result, io2; kwargs...)
end

# ╔═╡ 4b966743-9e82-46e7-843a-c8d33d3cb2e4
length(volts)

# ╔═╡ b81676e8-dcec-49fd-b350-f26ac61243ec
let

    vis = GridVisualizer(;
        size = (600, 400),
        tilte = "IV Curve",
        xlabel = "Φ_WE/V",
        ylabel = "I",
        legend = :lb,
    )
    scalarplot!(
        vis,
        volts,
        currents(result, ihplus; electrode = :bulk),
        linestyle = "--",
        label = "H+, bulk",
        color = :red,
    )
    scalarplot!(
        vis,
        volts,
        currents(result, ihplus),
        color = :red,
        clear = false,
        linestyle = :solid,
        label = "H+, we",
    )

    scalarplot!(
        vis,
        volts,
        currents(result, io2; electrode = :bulk),
        linestyle = "--",
        label = "O2, bulk",
        color = :green,
        clear = false,
    )
    scalarplot!(
        vis,
        volts,
        currents(result, io2),
        color = :green,
        clear = false,
        linestyle = :solid,
        label = "O2, we",
    )

    scalarplot!(
        vis,
        volts,
        currents(result, io2, electrode = :bulk),
        linestyle = "--",
        label = "O2, bulk",
        color = :green,
        clear = false,
    )
    scalarplot!(
        vis,
        volts,
        currents(result, io2),
        color = :green,
        clear = false,
        linestyle = :solid,
        label = "O2, we",
    )


    scalarplot!(
        vis,
        volts,
        currents_h2o(result; electrode = :bulk) / 100,
        linestyle = "--",
        label = "H2O/100, bulk",
        color = :blue,
        clear = false,
    )
    scalarplot!(
        vis,
        volts,
        currents_h2o(result) / 100,
        color = :blue,
        clear = false,
        linestyle = :solid,
        label = "H2O/100, we",
    )


    reveal(vis)
end

# ╔═╡ d1de762f-5f84-4bfb-bb77-2c5c61b5411e
ix = ihplus

# ╔═╡ ab97b97f-a7ad-4fc5-aa35-8876e8d79eea
let
    jwe = currents(result, ix, electrode = :we)
    jbulk = currents(result, ix, electrode = :bulk)
    scalarplot(
        volts,
        (jwe - jbulk),
        size = (600, 200),
    )
end

# ╔═╡ d7b10140-7db7-4be0-88c3-53ba1f203310
@bind vshow PlutoUI.Slider(range(vmin, vmax, length = 101), show_value = true)

# ╔═╡ 9226027b-725d-446e-bc14-dd335a60ec09
let
    vinter = linear_interpolation(volts, currents(result, io2))
    sol = tsol(vshow)
    c0 = solventconcentration(sol * (mol / dm^3), celldata) / (mol / dm^3)
    check =
        celldata.v0 * c0 +
        celldata.v[io2] * sol[io2, :] +
        (celldata.v[ihplus] + κ * celldata.v0) * sol[ihplus, :] +
        (celldata.v[iso4] + κ * celldata.v0) * sol[iso4, :]
    @info extrema(check)
    title = "Φ_we=$(vshow), I=$(round(vinter(vshow), sigdigits = 3))"
    vis = GridVisualizer(;
        size = (600, 250),
        yscale = :log,
        limits = (1.0e-6, 100),
        legend = :rt,
        title,
    )
    scalarplot!(vis, grid, sol[io2, :]; color = :green, label = "O_2")
    scalarplot!(vis, grid, sol[iso4, :]; color = :gray, clear = false, label = "SO4--")
    scalarplot!(vis, grid, sol[ihplus, :]; color = :red, clear = false, label = "H+")
    scalarplot!(vis, grid, c0; color = :blue, clear = false, label = "H2O")
    reveal(vis)
end

# ╔═╡ 201e1cda-0090-40cc-94a2-2d547cc1351e
md"""
## Alternative cell with explicit solvent
"""

# ╔═╡ d4b327a7-4a6f-4aad-994c-419365691b1c
md"""
This does not work and is not correct: there is no binary diffusion coefficient
between ih2o and the other solutes. 
"""

# ╔═╡ f6d0dc7c-ecea-4915-aaa6-cd20a749fafa
begin
    ih2o = 4
    xz = [1, -2, 0, 0]
    xκ_all = [κ, κ, 0, 0]
end

# ╔═╡ 1e33f0a5-676c-40ab-89ab-b410bcdd9b39
function xhalfcellbc(f, u, bnode, data)
    bulkbcondition(f, u, bnode, data)
    (; iϕ, eneutral, ϕ_we, Γ_we, RT) = data

    return if bnode.region == Γ_we
        f .= 0.0
        if !eneutral
            boundary_dirichlet!(
                f,
                u,
                bnode;
                species = iϕ,
                region = data.Γ_we,
                value = data.ϕ_we,
            )
        end
        μsolv, μ = @inline chemical_potentials!(MVector{5, eltype(u)}(undef), u, data)
        A =
            (
            4 * μ[ihplus] + μ[io2] - 2μ[ih2o] +
                Δg +
                eneutral * 4 * F * (u[iϕ] - data.ϕ_we)
        ) / (RT)
        r = rrate(R0, β, A)
        f[ihplus] -= 4 * r
        f[io2] -= r
        f[ih2o] += 2r + 4κ * r # shedding of the solvation shell of ihplus
    end
end


# ╔═╡ 8a83ef93-415b-4fe9-97cb-6e381382059e
begin
    xcelldata = ElectrolyteData(;
        nc = 4,
        z = xz,
        κ = xκ_all,
        Γ_we = 1,
        Γ_bulk = 2,
        eneutral = false,
        scheme,
        epsreg,
    )

    xc_bulk = xcelldata.c_bulk

    c_slv = 1.0e-6 * mol / dm^3
    xc_bulk[io2] = 0.001 * mol / dm^3
    xc_bulk[iso4] = molarity * mol / dm^3
    xc_bulk[ihplus] = 2.0 * molarity * mol / dm^3

    xc_bulk[ih2o] =
        55.4 * mol / dm^3 - xc_bulk[io2] - (1 + κ) * xc_bulk[iso4] -
        (1 + κ) * xc_bulk[ihplus] - c_slv


    @assert isapprox(xcelldata.c_bulk' * xcelldata.z, 0, atol = 1.0e-12)

    xcell = PNPSystem(grid; bcondition = xhalfcellbc, celldata = xcelldata)


end

# ╔═╡ 6e3cf2f4-9d4e-4611-8702-9143d85c3861
let
    @info xc_bulk
    c0, barc = c0_barc(xc_bulk, xcelldata)
    @info c0, barc, xc_bulk[ih2o]
    v0 = xcelldata.v0
    v = xcelldata.v
    κx = xcelldata.κ
    cv = c0 * v0
    for i in 1:xcelldata.nc
        cv += xc_bulk[i] * (v[i] + κx[i] * v0)
    end
    cv
end

# ╔═╡ 4c537d80-e6e2-4fc2-b388-9e277c5b67b2
begin
    #    xvolts, xj_we, xj_bulk, xsols

    xresult = LiquidElectrolytes.ivsweep(
        xcell;
        voltages = vmin:vdelta:vmax,
        store_solutions = true,
        solver_control...,
    )
    xtsol = voltages_solutions(xresult)

    for it in 1:length(xtsol.t)
        xtsol.u[it][io2, :] /= mol / dm^3
        xtsol.u[it][ihplus, :] /= mol / dm^3
        xtsol.u[it][iso4, :] /= mol / dm^3
        xtsol.u[it][ih2o, :] /= mol / dm^3
    end
end

# ╔═╡ 2738089b-8fb1-4cd7-8f3f-87a5086a13a7
@bind xvshow PlutoUI.Slider(range(vmin, vmax, length = 101), show_value = true)

# ╔═╡ c1098a16-c62c-4eeb-b9be-bda21f6b2bca
let
    sol = xtsol(xvshow)
    vinter = linear_interpolation(xresult.voltages, currents(xresult, io2))
    c0 = solventconcentration(sol * (mol / dm^3), xcelldata) / (mol / dm^3)
    title = "Φ_we=$(vshow), I=$(round(vinter(vshow), sigdigits = 3))"

    vis = GridVisualizer(;
        size = (600, 250),
        yscale = :log,
        limits = (1.0e-6, 100),
        legend = :rt,
        title,
    )

    check =
        xcelldata.v0 * c0 +
        xcelldata.v[io2] * sol[io2, :] +
        (xcelldata.v[ihplus] + κ * xcelldata.v0) * sol[ihplus, :] +
        (xcelldata.v[iso4] + κ * xcelldata.v0) * sol[iso4, :] +
        (xcelldata.v[ih2o]) * sol[ih2o, :]
    @info extrema(check)

    scalarplot!(vis, grid, sol[io2, :]; color = :green, label = "O_2")
    scalarplot!(vis, grid, sol[iso4, :]; color = :gray, clear = false, label = "SO4--")
    scalarplot!(vis, grid, sol[ihplus, :]; color = :red, clear = false, label = "H+")
    scalarplot!(vis, grid, sol[ih2o, :]; color = :blue, clear = false, label = "H2O")
    scalarplot!(vis, grid, c0; color = :magenta, clear = false, label = "H2O")
    reveal(vis)
end

# ╔═╡ e74e788f-9f84-44a8-9846-17e53ce3b3be
@bind yvshow PlutoUI.Slider(range(vmin, vmax, length = 101), show_value = true)

# ╔═╡ 8d0910d4-bfbb-4659-ab6e-9ff0e696722d
let
    vis = GridVisualizer(; size = (600, 200), legend = :rt)
    j_bulk = currents(result, io2)
    x_bulk = currents(xresult, io2)
    sol = tsol(yvshow)
    c0 = solventconcentration(sol * (mol / dm^3), celldata) / (mol / dm^3)
    #    j = linear_interpolation(volts, j_bulk)(yvshow)
    #    jh2o = -j[ihplus] - j[iso4] - j[io2]
    #    @info j[1:3], jh2o
    #    @info linear_interpolation(xvolts, xj_bulk)(yvshow)[1:4]
    xsol = xtsol(yvshow)
    scalarplot!(vis, grid, xsol[ih2o, :], color = :green, label = "explicit")
    scalarplot!(vis, grid, c0, color = :red, label = "implicit", clear = false)
    reveal(vis)
end

# ╔═╡ f9b4d4dc-7def-409f-b40a-f4eba1163741
TableOfContents()

# ╔═╡ 7a93e9a8-8a2d-4b11-84ef-691706c0eb0f
begin
    hrule() = html"""<hr>"""
    highlight(mdstring, color) =
    htl"""<blockquote style="padding: 10px; background-color: $(color);">$(mdstring)</blockquote>"""

    macro important_str(s)
        return :(highlight(Markdown.parse($s), "#ffcccc"))
    end
    macro definition_str(s)
        return :(highlight(Markdown.parse($s), "#ccccff"))
    end
    macro statement_str(s)
        return :(highlight(Markdown.parse($s), "#ccffcc"))
    end


    html"""
            <style>
             h1{background-color:#dddddd;  padding: 10px;}
             h2{background-color:#e7e7e7;  padding: 10px;}
             h3{background-color:#eeeeee;  padding: 10px;}
             h4{background-color:#f7f7f7;  padding: 10px;}
            
    	     pluto-log-dot-sizer  { max-width: 655px;}
             pluto-log-dot.Stdout { background: #002000;
    	                            color: #10f080;
                                    border: 6px solid #b7b7b7;
                                    min-width: 18em;
                                    max-height: 300px;
                                    width: 675px;
                                        overflow: auto;
     	                           }
    	
        </style>
    """
end

# ╔═╡ 5beb3a0d-e57a-4aea-b7a0-59b8ce9ff5ce
hrule()


# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╠═78685b0c-b27a-450a-9bae-5948c578d688
# ╟─84e05551-7d51-4b2c-88f2-b186ad6a244a
# ╟─16b9af50-11c5-4bdf-b5d8-8d2d9331b5e9
# ╠═515baffb-2e22-401d-aacd-15971dd4365e
# ╟─bc896b9c-03a9-4da0-be80-f096766cb039
# ╠═50ff0114-f9b1-4bd0-8190-153fef07ef3a
# ╟─231ebfc2-51b6-4027-8be2-75891dfed7e0
# ╠═a4e0379d-f7f6-4b61-bf38-5eb17f67505a
# ╟─970389b5-d2c1-4992-9978-aca1ccd3d2fc
# ╠═5482615c-cb5b-4c44-99c4-19416eabae7f
# ╠═136cf1aa-75c6-4aaa-b93b-801391ec800c
# ╠═b916eb92-8ba8-49aa-bd7c-1bfc91e813d4
# ╠═392a648c-12a2-41fb-b3f6-aa3cfe3cbcd7
# ╟─420a82e0-5fc2-47ea-8916-d88910655d50
# ╠═60d85155-9fa7-4740-9769-212ceef1918b
# ╠═bd3b54ce-de22-4f8b-aa27-7f56c5dca140
# ╠═8683cf81-4b53-4099-a668-37c3641802dd
# ╠═1d23b6b5-88cf-4200-b612-82dffcc5cca7
# ╠═763c393c-e0c8-447e-92e4-f2a5f0de2a30
# ╟─2717060b-9e75-439f-8aa8-614f0930a163
# ╠═5bc4f11f-24c6-4af8-a554-1b5771f1f2b0
# ╠═4b966743-9e82-46e7-843a-c8d33d3cb2e4
# ╠═b81676e8-dcec-49fd-b350-f26ac61243ec
# ╠═d1de762f-5f84-4bfb-bb77-2c5c61b5411e
# ╠═ab97b97f-a7ad-4fc5-aa35-8876e8d79eea
# ╠═d7b10140-7db7-4be0-88c3-53ba1f203310
# ╠═9226027b-725d-446e-bc14-dd335a60ec09
# ╟─201e1cda-0090-40cc-94a2-2d547cc1351e
# ╠═d4b327a7-4a6f-4aad-994c-419365691b1c
# ╠═f6d0dc7c-ecea-4915-aaa6-cd20a749fafa
# ╠═1e33f0a5-676c-40ab-89ab-b410bcdd9b39
# ╠═8a83ef93-415b-4fe9-97cb-6e381382059e
# ╠═6e3cf2f4-9d4e-4611-8702-9143d85c3861
# ╠═4c537d80-e6e2-4fc2-b388-9e277c5b67b2
# ╠═2738089b-8fb1-4cd7-8f3f-87a5086a13a7
# ╠═c1098a16-c62c-4eeb-b9be-bda21f6b2bca
# ╟─e74e788f-9f84-44a8-9846-17e53ce3b3be
# ╟─8d0910d4-bfbb-4659-ab6e-9ff0e696722d
# ╟─5beb3a0d-e57a-4aea-b7a0-59b8ce9ff5ce
# ╟─f9b4d4dc-7def-409f-b40a-f4eba1163741
# ╟─7a93e9a8-8a2d-4b11-84ef-691706c0eb0f
