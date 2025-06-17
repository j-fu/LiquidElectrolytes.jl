### A Pluto.jl notebook ###
# v0.20.9

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

# ╔═╡ a70cef7d-2a2f-4155-bdf3-fec9df94c63f
begin
    import Pkg # hide
    Pkg.activate(joinpath(@__DIR__, "..", "docs")) # hide
    using Revise # hide
    using PlutoUI, HypertextLiteral
    using LiquidElectrolytes
    using Printf
    using LessUnitful
    using ExtendableGrids
    using VoronoiFVM
    using GridVisualize
    using LinearSolve, ExtendableSparse
    using StaticArrays
    using Interpolations
    using LaTeXStrings
    using DrWatson
    using Test
    using DoubleFloats

    #VoronoiFVM.set_logging!(false)
    if isdefined(Main, :PlutoRunner)
        using CairoMakie, Colors

        default_plotter!(CairoMakie)
        CairoMakie.activate!(type = "svg")
        CairoMakie.set_theme!(CairoMakie.theme_latexfonts(), fontsize = 12)
    end
end


# ╔═╡ baa3e08e-5d64-4c8f-9f6d-5fdb40e97bc5
begin
    using HypertextLiteral: @htl_str, @htl
    using UUIDs: uuid1
end

# ╔═╡ 954648a2-b1fa-4934-a9f7-914dc77101eb
TableOfContents(aside = false)

# ╔═╡ 1e7e4a5d-3ea0-4a3f-a3b3-cc5c041cbfe6
md"""
## Data
"""

# ╔═╡ f15c4585-a2eb-4278-aedf-0ebfde6f40fd
md"""
```math
	O + ne^- \leftrightharpoons R
```
"""

# ╔═╡ a201c02d-8663-4d8b-9e17-f780676261c6
begin
    @phconstants N_A e R ε_0
    const F = N_A * e
    @unitfactors cm μF mol dm s mA A nm C
end;

# ╔═╡ bccf3c34-0e09-44a9-b195-55c7f82d416d
md"""
We have  ``\int_0^\infty be^{-bx}= 1``
"""

# ╔═╡ b1a6fb55-e9c1-40c6-ba79-c93b760f8fe5
md"""
## Result plots
"""

# ╔═╡ 90f3f622-be77-4107-b8f5-87712a8dd560
md"""
## Summary
"""

# ╔═╡ 7b62d217-aa2a-4fff-a96d-f75e201612a7
md"""
- With negative applied potential, O which is more positive can take in the available electrons and react to R. Vice versa, with positive potential, we "suck" the electrons from the more negative R, leaving the more positive O behind.

- Double-layer-wise, the attraction to the electrode can go into the opposite direction... 
- We have the consistency with the electroneutral case, there, we see the reaction direction well
- Can we calculate effective diffusivities in the double layer ?
- There are two effects limiting the current in the double layer
  - Displacement by the background electrolyte
  - Electrostatic repulsion
- Further verification:
  - Monotonicity of electrochemical potential
  - Zero electrochem. potential with  R=0
- Errors introduced by interpolation between solutions are comparably large and may result in non-monotone electrochemical potentials
"""

# ╔═╡ 93ba3568-fe21-4371-beaa-e4d05cb45797
md"""
## Implementation + Solution
"""

# ╔═╡ cfbf88e4-db45-4653-bfcf-ce39b5ec86e1
function tunnelbc(f, u, bnode, data)
    (; nc, Γ_we, Γ_bulk, ϕ_we, ip, iϕ, v, v0, RT) = data
    bulkbcondition(f, u, bnode, data; region = Γ_bulk)
    return if bnode.region == Γ_we
        boundary_dirichlet!(f, u, bnode; species = iϕ, region = Γ_we, value = ϕ_we)
    end
end

# ╔═╡ 3bf6f385-20f4-49b1-a791-ff54ef87b575
md"""
## Experimental Code
"""

# ╔═╡ 40567761-e66d-4fa9-93f1-4ffc38725488
function ismonotone(u, i)
    ∇u = u[i, 2:end] - u[i, 1:(end - 1)]
    return all(x -> x > 0, ∇u) || all(x -> x < 0, ∇u)
end

# ╔═╡ 96ee98b8-22c3-4b58-bd30-45aec74d95ae
function flux(u, data)
    (; nc) = data
    nn = size(u, 2)
    flx = zeros(size(u, 1), nn - 1)
    for i in 1:(nn - 1)
        @views LiquidElectrolytes.pnpflux(flx[:, i], u[:, [i + 1, i]], (index = i,), data)
    end
    return flx
end


# ╔═╡ 8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
html"""<hr>"""

# ╔═╡ 5998bd78-63bf-4c7c-8392-38e6342ee961
md"""
## Appendix: Technicalities
"""

# ╔═╡ 5944b4ff-a75c-4487-ba69-2a28297d4d83
pdir(args...) = projectdir("..", args...)

# ╔═╡ 5c4510c7-e721-497e-8c81-f50f850a93af
function zstr(z)
    if z == -2
        return "2-"
    elseif z == -1
        return "-"
    elseif z == 0
        return ""
    elseif z == 1
        return "+"
    elseif z == 2
        return "2+"
    elseif z == 3
        return "3+"
    else
        error("zstr($(z)) not implemented")
    end
end

# ╔═╡ 784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
html"""<style>.dont-panic{ display: none }</style>"""

# ╔═╡ afe4745f-f9f1-4e23-8735-cbec6fb79c41
begin
    function floataside(text::Markdown.MD; top = 1)
        uuid = uuid1()
        return @htl(
            """
            		<style>


            		@media (min-width: calc(700px + 30px + 300px)) {
            			aside.plutoui-aside-wrapper-$(uuid) {

            	color: var(--pluto-output-color);
            	position:fixed;
            	right: 1rem;
            	top: $(top)px;
            	width: 400px;
            	padding: 10px;
            	border: 3px solid rgba(0, 0, 0, 0.15);
            	border-radius: 10px;
            	box-shadow: 0 0 11px 0px #00000010;
            	/* That is, viewport minus top minus Live Docs */
            	max-height: calc(100vh - 5rem - 56px);
            	overflow: auto;
            	z-index: 40;
            	background-color: var(--main-bg-color);
            	transition: transform 300ms cubic-bezier(0.18, 0.89, 0.45, 1.12);

            			}
            			aside.plutoui-aside-wrapper > div {
            #				width: 300px;
            			}
            		}
            		</style>

            		<aside class="plutoui-aside-wrapper-$(uuid)">
            		<div>
            		$(text)
            		</div>
            		</aside>

            		"""
        )
    end
    floataside(stuff; kwargs...) = floataside(md"""$(stuff)"""; kwargs...)
end;


# ╔═╡ efd397bd-9b4d-402f-ad42-cc0f19dc48b6

floataside(
    @bind guidata confirm(
        #! format: off
        PlutoUI.combine() do Child
md"""
  __User Data__	``\quad O + ne^- \leftrightharpoons R``
 - ``z_R``: $(Child("zR", NumberField(-2:2;default=-1))) 
   ``\quad n:`` $(Child("n", NumberField(0:2;default=1)))
   ``\quad κ:`` $(Child("κ", NumberField(0:10;default=10)))
 - ``M/(mol/L)`` $(Child("fgmol",NumberField(0.1:0.1:4;default=0.1)))
   ``M_{bg}/(mol/L)`` $(Child("bgmol",NumberField(0.1:0.1:4;default=2)))
 - Double64: $(Child("double64",CheckBox()))
 - tunnel: $(Child("tunnel", CheckBox()))
   ``\quad β/cm^{-1}``: $(Child("β", TextField(10;default="1.0e8")))
"""
        end,
        #! format: on
        label = "Submit"
    );
    top = 40
)


# ╔═╡ 08bd3963-eb07-460c-aa92-862d82be4e70
begin
    verbose = ""
    nref = 0
    voltages = (-1:0.025:1) * ufac"V"
    dlcap = false
    const molR = guidata.fgmol
    const molO = guidata.fgmol
    const molS = guidata.bgmol
    const iO = 1
    const iR = 2
    const iS = 3
    const iX = 4
    const iϕ = 5
    const zR = guidata.zR
    const zO = zR + guidata.n
    const zS = -1
    const zX = 1

    const molX = molS - zR * molR - zO * molO
    const z = [zO, zR, zS, zX]
    const c_bulk = [molO, molR, molS, molX] * mol / dm^3
    @info c_bulk' * z
    @assert abs(c_bulk' * z < 1.0e-12)
    const κ = Float64(guidata.κ)
    const R0 = 1.0e-10 * mol / (cm^2 * s)
    const Δg = 0.0
    const β = 0.5
    const βt = parse(Float64, guidata.β) / cm
    etunnel(x) = βt * exp(-βt * x)

    pnpdata = ElectrolyteData(;
        nc = 4,
        z,
        c_bulk,
        κ = fill(κ, 4),
        Γ_we = 1,
        Γ_bulk = 2,
    )


    pnpdata
end

# ╔═╡ c30965b0-63da-48c4-8d63-63469a4c8fe2
begin
    hmin = 1.0e-2 * nm * 2.0^(-nref)
    L = 20.0 * nm
    hmax = L / 20 * 2.0^(-nref)
    X = geomspace(0, L, hmin, hmax)
    grid = simplexgrid(X)
end

# ╔═╡ 605520f0-8ff1-4cd8-b634-c2a0d5d38675
gridplot(grid, size = (600, 150))

# ╔═╡ 0791e740-a47c-4c3e-9d2c-b0c5daf7a403
function halfcellbc(f, u, bnode, data)
    (; nc, Γ_we, Γ_bulk, ϕ_we, ip, iϕ, v, v0, RT, κ) = data
    bulkbcondition(f, u, bnode, data; region = Γ_bulk)
    if bnode.region == Γ_we
        if !data.eneutral
            boundary_dirichlet!(f, u, bnode; species = iϕ, region = Γ_we, value = ϕ_we)
        end
        c0, barc = c0_barc(u, data)
        μR = chemical_potential(u[iR], barc, u[ip], v[iR] + κ[iR] * v0, data)
        μO = chemical_potential(u[iO], barc, u[ip], v[iO] + κ[iO] * v0, data)
        A = (μR - μO + Δg + data.eneutral * (zR - zO) * F * (u[iϕ] - ϕ_we)) / RT
        r = rrate(R0, β, A)
        f[iR] -= r
        f[iO] += r
    end
    return nothing
end


# ╔═╡ f38923b3-b25e-4a2c-a502-91931e522019
function tunnelreaction(f, u, node, data)
    (; nc, Γ_we, Γ_bulk, ϕ_we, ip, iϕ, v, v0, RT) = data

    c0, barc = c0_barc(u, data)
    x = node[1]
    g = etunnel(x)
    μR = chemical_potential(u[iR], barc, u[ip], v[iR] + κ * v0, data)
    μO = chemical_potential(u[iO], barc, u[ip], v[iO] + κ * v0, data)
    A = (μR - μO) / RT
    r = rrate(g * R0, β, A)
    f[iR] -= r
    f[iO] += r
    return nothing
end


# ╔═╡ 258f970b-5fe8-46b6-9c81-f84524e3e6b0
ZR = zstr(zR)

# ╔═╡ e6150334-31f8-450d-8a5d-da6e17ab81c1
ZO = zstr(zO)

# ╔═╡ ea14cd1b-5a64-4911-80b2-fde0306120fe
valuetype = guidata.double64 ? Double64 : Float64

# ╔═╡ 793265dd-bfeb-4083-a2fa-80e74f06fdb2
function sweep(; eneutral = true, tunnel = false, bikerman = true)
    celldata = deepcopy(pnpdata)
    celldata.eneutral = eneutral
    if tunnel
        pnpcell = PNPSystem(grid; bcondition = tunnelbc, reaction = tunnelreaction, celldata)
    else
        pnpcell = PNPSystem(grid; bcondition = halfcellbc, celldata, valuetype)
    end
    result = ivsweep(
        pnpcell; voltages, store_solutions = true, verbose = "en",
        maxiters = 25,
        Δu_opt = 1.0e-3,
        Δp_min = 1.0e-4,
        abstol = 1.0e-12,
        reltol = 1.0e-6,
        tol_round = 1.0e-8,
        max_round = 5,
        damp_initial = 0.5
    )
    currs = (zO - zR) * LiquidElectrolytes.currents(result, iR, electrode = :bulk) #zO*LiquidElectrolytes.currents(result,iO, electrode=:bulk)
    sol = LiquidElectrolytes.voltages_solutions(result)
    q = chargedensity(sol, celldata)
    return result, currs, sol, q, celldata
end

# ╔═╡ 9c3e818f-5828-44e3-83a3-4986914feb85
pnpresult, pnpcurrs, pnpsol, pnpcharge, pnpd = sweep(eneutral = false, tunnel = false);

# ╔═╡ 52579e12-9658-4cfc-b014-841166245bc9
nnpresult, nnpcurrs, nnpsol, nnpcharge = sweep(eneutral = true, tunnel = false);

# ╔═╡ d6965320-b46a-4c07-8e73-dd9369690506
if guidata.tunnel
    tresult, tcurrs, tsol = sweep(eneutral = false, tunnel = true)
else
    tresult, tcurrs, tsol = (nothing, nothing, nothing)
    tcharge = nothing
end

# ╔═╡ c26ad7b4-d6f7-4ad1-844c-78c0aa7d1b3d
begin
    mylog = RLog(valuetype)
    mylog10(x) = mylog(x) / log(10)
    Makie.defined_interval(::typeof(mylog10)) = Makie.defined_interval(log10)
    Makie.defaultlimits(::T) where {T <: typeof(mylog10)} = Makie.defaultlimits(log10)
end


# ╔═╡ 142bbf36-3616-4c57-81ba-010816861419
floataside(
    @bind plotscale confirm(
        PlutoUI.combine() do Child
            md"""
            __Plot scaling__

            S: $(Child("S",TextField(8;default="1")))
            M: $(Child("M",TextField(8;default="1.0e2")))
            L: $(Child("L",TextField(8;default="1.0e5")))
            """
        end,
        label = "Submit"
    );
    top = 250
)

# ╔═╡ 07987612-58df-4a49-a929-55ce191681a1
floataside(
    md"""
    __Voltage:__ $(@bind iv PlutoUI.Slider(1:length(pnpresult.voltages), show_value=false))
    """, top = 380
)

# ╔═╡ 2847be38-9480-48f0-b4a9-b781541c2cdd
begin
    v = pnpresult.voltages[iv]
    floataside(
        md"""
        Applied Voltage: v=$(round(v, sigdigits=4))V
        """, top = 450
    )
end

# ╔═╡ af5c26e0-026e-4ddd-98c7-94bf4e534256
if isdefined(Main, :PlutoRunner)
    let
        function yplot(fg, imax; xlabel = "")
            ax = Axis(fg; ylabel = L"I/(mA/cm^2)", xlabel)
            ylims!(ax, [-imax, imax])
            # Plotter.scatter!(ax,pnpresult.voltages[1:12:end], -pnpcurrs[1:12:end]/(mA/cm^2),
            #          color = :red, markershape = :utriangle, markersize = 3)

            lines!(
                ax,
                nnpresult.voltages, -nnpcurrs / (mA / cm^2),
                color = :green,
                linewidth = 2
            )
            scatterlines!(
                ax,
                nnpresult.voltages[1:10:end], -nnpcurrs[1:10:end] / (mA / cm^2),
                color = :green, label = "NNP",
                linewidth = 2
            )
            lines!(
                ax,
                pnpresult.voltages, -pnpcurrs / (mA / cm^2), color = :red, label = "PNP", linewidth = 2
            )
            lines!(ax, [v, v], [-imax, imax], color = :gray)

            if guidata.tunnel
                lines!(ax, tresult.voltages, -tcurrs / (mA / cm^2), color = :blue, label = "PNP-tunnel")
            end
            return ax
        end
        fig = Figure(size = (650, 400), fontsize = 12)
        title = guidata.tunnel ? L"z_R=%$(zR), n=%$(guidata.n), β=%$(guidata.β)/cm" : L"z_R=%$(zR), n=%$(guidata.n)"
        ax1 = yplot(fig[1, :], parse(Float64, plotscale.S))
        ax2 = yplot(fig[2, :], parse(Float64, plotscale.M))
        ax3 = yplot(fig[3, :], parse(Float64, plotscale.L); xlabel = L"U/V")
        Label(fig[0, :], title, tellwidth = false, fontsize = 16)
        Legend(
            fig[1:3, 2], ax2; labelsize = 10,
            backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.5)
        )

        if false
            if guidata.tunnel
                CairoMakie.save(
                    savename(
                        pdir("redox_tunnel"),
                        (zR = zR, n = guidata.n, beta = guidata.β),
                        "png"
                    ), fig
                )
            else
                CairoMakie.save(savename(pdir("redox"), (zR = zR, n = n), "png"), fig)
            end
        end
        fig
    end
end

# ╔═╡ ae0af916-986f-42cf-8a7a-215171ee5a55
begin
    tu = nothing
    if guidata.tunnel
        tu = tsol(v)
    end
    if !isnothing(tu)
        tu = max.(1.0e-15, tu / (mol / dm^3))
        μtu = LiquidElectrolytes.electrochemical_potentials(tsol(v), pnpd)
    end
    uu = pnpsol[:, :, iv]
    ru = uu / (mol / dm^3)
    nu = nnpsol(v) / (mol / dm^3)
    #    _, μu = LiquidElectrolytes.chemical_potentials(pnpsol(v), pnpd)
    μu = LiquidElectrolytes.electrochemical_potentials(pnpsol(v), pnpd)
end

# ╔═╡ 6c7099a1-dfad-4fd4-b940-bb60748f58fa
# ╠═╡ show_logs = false
if isdefined(Main, :PlutoRunner)
    let
        XX = (X[2:end] / nm)
        if v < 0
            dir = "O\\to R,"
        elseif v ≈ 0
            dir = ""
        else
            dir = "R\\to O,"
        end

        title = guidata.tunnel ?
            L"U=%$(round(v,digits=4))V,\, %$(dir)\, z_R=%$(zR),\, z_O=%$(zO),\, β=%$(guidata.β)/cm" :
            L"U=%$(round(v,digits=4))V,\, %$(dir)\, z_R=%$(zR),\, z_O=%$(zO)"


        fig = Figure(size = (650, 400))
        ax1 = Axis(
            fig[1, 1];
            #       xlabel=L"x/nm",
            ylabel = L"c/(mol/dm^3)",
            xscale = mylog10,
            yscale = log10
        )
        xlims!(ax1, XX[1], L / nm)
        ylims!(ax1, 1.0e-10, 1.0e2)

        lines!(ax1, XX, ru[iO, 2:end], color = :red, linestyle = :solid, label = L"O^{%$(ZO)}, pnp")
        lines!(ax1, XX, ru[iR, 2:end], color = :blue, linestyle = :solid, label = L"R^{%$(ZR)}, pnp")
        lines!(ax1, XX, ru[iS, 2:end], color = :lightgreen, linestyle = :solid, label = L"S^-, pnp")
        lines!(ax1, XX, ru[iX, 2:end], color = :orange, linestyle = :solid, label = L"X^+, pnp")
        lines!(ax1, XX, nu[iO, 2:end], color = :red, linestyle = :dash, label = L"O^{%$(ZO)}, nnp")
        lines!(ax1, XX, nu[iR, 2:end], color = :blue, linestyle = :dash, label = L"R^{%$(ZR)}, nnp")
        lines!(ax1, XX, abs.(pnpcharge(v)[1, 2:end] / (e / nm^3)), color = :black, linestyle = :solid, label = L"q")

        if !isnothing(tu)
            lines!(ax1, XX, tu[iO, 2:end], color = :red, linestyle = :dashdot, linewidth = 3, label = L"O^{%$(ZO)}, pnp-t")
            lines!(ax1, XX, tu[iR, 2:end], color = :blue, linestyle = :dashdot, linewidth = 3, label = L"R^{%$(ZR)}, pnp-t")
            lines!(
                ax1, XX, etunnel.(X[2:end]) / etunnel(0),
                color = :gray, linestyle = :solid
            )
        end
        Legend(
            fig[1, 2], ax1; labelsize = 10,
            backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.5)
        )

        ################################################################################
        ax2 = Axis(
            fig[2, 1];
            #    yscale=log10,
            xscale = log10,
            xlabel = L"x/nm",
            ylabel = L"μ/V",
        )
        #limits=(-2.5e4, 1.0e5),
        xlims!(1.0e-2, L / nm)
        lines!(ax2, XX, μu[iO, 2:end], color = :red, linestyle = :solid, label = L"O^{%$(ZO)}, pnp")
        lines!(ax2, XX, μu[iR, 2:end], color = :blue, linestyle = :solid, label = L"R^{%$(ZR)}, pnp")
        if !isnothing(tu)
            lines!(ax2, XX, μtu[iO, 2:end], color = :red, linestyle = :dashdot, label = L"O^{%$(ZO)}, pnp-t", linewidth = 3)
            lines!(ax2, XX, μtu[iR, 2:end], color = :blue, linestyle = :dashdot, label = L"R^{%$(ZR)}, pnp-t", linewidth = 3)
        end
        #lines!(ax2,XX,uu[iϕ,2:end], color=:green,linestyle=:solid, label=L"\phi, pnp"),
        #lines!(ax2,XX,μu[iS,2:end], color=:lightgreen,linestyle=:solid, label=L"S^-, pnp", clear=false)
        #lines!(ax2,XX,μu[iX,2:end], color=:orange,linestyle=:solid, label=L"X^+, pnp")
        Legend(
            fig[2, 2], ax2; labelsize = 10,
            backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.5)
        )
        Label(fig[0, :], title, tellwidth = false, fontsize = 16)

        fig
    end
end

# ╔═╡ 33f0cc8b-099a-4963-864f-84a6e937b921
cx = sum(uu[1:4, :], dims = 1) / (mol / dm^3)

# ╔═╡ e99c4c5f-474a-4eb7-9599-0e51d60fb08c
[c0_barc(uu[:, i], pnpdata)[2] for i in 1:size(uu, 2)] / (mol / dm^3)

# ╔═╡ 42b2d17d-95cd-4e40-8a5f-a17e35319e05
μu

# ╔═╡ c6ac1bfd-1720-4c6b-9cf3-2b1f68012394
ismonotone(μu, 1) && ismonotone(μu, 2)

# ╔═╡ 7bd39e66-089a-4a63-af1f-5aa90b7427a1
μu[iO, :] .- μu[iO, end]

# ╔═╡ 0f55903e-7529-4f63-8de2-1d6e772e6bca
begin
    flx = pnpdata.F * flux(ru, pnpd)
    for i in 1:size(flx, 2)
        for j in 1:size(flx, 1)
            flx[j, i] /= (X[i + 1] - X[i])
        end
    end
    abs.((flx[iR, :] + flx[iO, :]) ./ flx[iO, :])
end

# ╔═╡ a9d4ea79-8e78-462d-a09c-4a798e5da725
function xflux(i)
    dμ = μu[i, 1:(end - 1)] - μu[i, 2:end]
    ac = pnpdata.F * pnpdata.D[i] * (ru[i, 1:(end - 1)] + ru[i, 2:end]) / 2
    return ac .* dμ ./ (X[2:end] - X[1:(end - 1)])
end

# ╔═╡ 999a18c7-1ff0-425b-aa92-239847373bb3
abs.(xflux(iR) + xflux(iO)) ./ abs.(xflux(iR))

# ╔═╡ Cell order:
# ╠═a70cef7d-2a2f-4155-bdf3-fec9df94c63f
# ╟─954648a2-b1fa-4934-a9f7-914dc77101eb
# ╟─1e7e4a5d-3ea0-4a3f-a3b3-cc5c041cbfe6
# ╠═f15c4585-a2eb-4278-aedf-0ebfde6f40fd
# ╠═a201c02d-8663-4d8b-9e17-f780676261c6
# ╠═c30965b0-63da-48c4-8d63-63469a4c8fe2
# ╟─605520f0-8ff1-4cd8-b634-c2a0d5d38675
# ╠═08bd3963-eb07-460c-aa92-862d82be4e70
# ╟─bccf3c34-0e09-44a9-b195-55c7f82d416d
# ╟─b1a6fb55-e9c1-40c6-ba79-c93b760f8fe5
# ╟─6c7099a1-dfad-4fd4-b940-bb60748f58fa
# ╟─af5c26e0-026e-4ddd-98c7-94bf4e534256
# ╟─90f3f622-be77-4107-b8f5-87712a8dd560
# ╠═7b62d217-aa2a-4fff-a96d-f75e201612a7
# ╟─93ba3568-fe21-4371-beaa-e4d05cb45797
# ╠═0791e740-a47c-4c3e-9d2c-b0c5daf7a403
# ╠═cfbf88e4-db45-4653-bfcf-ce39b5ec86e1
# ╠═f38923b3-b25e-4a2c-a502-91931e522019
# ╠═793265dd-bfeb-4083-a2fa-80e74f06fdb2
# ╠═9c3e818f-5828-44e3-83a3-4986914feb85
# ╠═52579e12-9658-4cfc-b014-841166245bc9
# ╠═d6965320-b46a-4c07-8e73-dd9369690506
# ╠═ae0af916-986f-42cf-8a7a-215171ee5a55
# ╟─3bf6f385-20f4-49b1-a791-ff54ef87b575
# ╠═40567761-e66d-4fa9-93f1-4ffc38725488
# ╟─33f0cc8b-099a-4963-864f-84a6e937b921
# ╠═e99c4c5f-474a-4eb7-9599-0e51d60fb08c
# ╠═42b2d17d-95cd-4e40-8a5f-a17e35319e05
# ╠═c6ac1bfd-1720-4c6b-9cf3-2b1f68012394
# ╠═7bd39e66-089a-4a63-af1f-5aa90b7427a1
# ╠═96ee98b8-22c3-4b58-bd30-45aec74d95ae
# ╠═0f55903e-7529-4f63-8de2-1d6e772e6bca
# ╠═a9d4ea79-8e78-462d-a09c-4a798e5da725
# ╠═999a18c7-1ff0-425b-aa92-239847373bb3
# ╟─8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
# ╟─5998bd78-63bf-4c7c-8392-38e6342ee961
# ╠═baa3e08e-5d64-4c8f-9f6d-5fdb40e97bc5
# ╠═5944b4ff-a75c-4487-ba69-2a28297d4d83
# ╠═c26ad7b4-d6f7-4ad1-844c-78c0aa7d1b3d
# ╠═efd397bd-9b4d-402f-ad42-cc0f19dc48b6
# ╠═5c4510c7-e721-497e-8c81-f50f850a93af
# ╠═258f970b-5fe8-46b6-9c81-f84524e3e6b0
# ╠═e6150334-31f8-450d-8a5d-da6e17ab81c1
# ╠═ea14cd1b-5a64-4911-80b2-fde0306120fe
# ╠═142bbf36-3616-4c57-81ba-010816861419
# ╠═07987612-58df-4a49-a929-55ce191681a1
# ╠═784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╠═2847be38-9480-48f0-b4a9-b781541c2cdd
# ╠═afe4745f-f9f1-4e23-8735-cbec6fb79c41
