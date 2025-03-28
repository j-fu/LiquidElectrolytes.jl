### A Pluto.jl notebook ###
# v0.20.5

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

# ╔═╡ 43b83ec3-3c85-4e2c-9a7d-e102ebc29619
begin
    @phconstants N_A e R ε_0
    const F = N_A * e
    @unitfactors cm μF mol dm s mA A nm C
end;


# ╔═╡ e73cfcbb-c3c4-49e2-a127-a505654aa8b4
function zstr(z::Int)
    if z == -1
        return "-"
    elseif z < -1
        return "$(-z)-"
    elseif z == 0
        return ""
    elseif z == 1
        return "+"
    elseif z > 1
        return "$(z)+"
    end
end


# ╔═╡ 8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
html"""<hr>"""

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


# ╔═╡ bb03204a-cb5a-4595-ace6-8b52a9440494
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
 - scanrate/``(V/s)``: $(Child("scanrate", TextField(6;default="0.1")))
   nperiods: $(Child("nperiods", NumberField(1:10;default=1)))
 - ``L/nm``: $(Child("L", TextField(10;default="20")))
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

# ╔═╡ 2496ac21-e027-49d2-a78d-964193d808e0
begin
    verbose = ""
    nref = 0
    sawtooth = SawTooth(
        scanrate = parse(Float64, guidata.scanrate),
        vmin = -0.5, vmax = 0.5
    )
    const nperiods = guidata.nperiods
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
    scheme = :μex
    #scheme = :act
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
        scheme
    )

    @assert iselectroneutral(c_bulk, pnpdata)
    @assert isincompressible(c_bulk, pnpdata)
    pnpdata
end

# ╔═╡ f0745600-be37-49d3-92fe-231bcc2e7013
function halfcellbc(f, u, bnode, data)
    (; nc, Γ_we, Γ_bulk, ip, iϕ, v, v0, RT, κ) = data
    bulkbcondition(f, u, bnode, data; region = Γ_bulk)
    ϕ_we = sawtooth(bnode.time)
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

# ╔═╡ 430f4195-7693-4682-8e81-4ef2d218a2f9
begin
    hmin = 1.0e-2 * nm * 2.0^(-nref)
    L = parse(Float64, guidata.L) * nm
    hmax = L / 20 * 2.0^(-nref)
    X = geomspace(0, L, hmin, hmax)
    grid = simplexgrid(X)
end

# ╔═╡ b218dbfc-c627-489c-8f27-8ebb7a721a44
gridplot(grid, size = (600, 150))

# ╔═╡ 483eefd0-a96d-4c32-a939-2517f08a6b3d
valuetype = guidata.double64 ? Double64 : Float64

# ╔═╡ 43cc60e3-8a4d-4fac-9454-7e599243e058
begin
    mylog = RLog(valuetype)
    mylog10(x) = mylog(x) / log(10)
    Makie.defined_interval(::typeof(mylog10)) = Makie.defined_interval(log10)
    Makie.defaultlimits(::T) where {T <: typeof(mylog10)} = Makie.defaultlimits(log10)

end

# ╔═╡ 76dab1ef-520c-41ea-b939-80e7212b27a9
function sweep(pnpdata; eneutral = true, tunnel = false, bikerman = true)
    celldata = deepcopy(pnpdata)
    celldata.eneutral = eneutral
    pnpcell = PNPSystem(grid; bcondition = halfcellbc, celldata, valuetype)
    return result = cvsweep(
        pnpcell;
        voltages = sawtooth,
        nperiods,
        store_solutions = true,
        verbose = "en"
    )

end

# ╔═╡ 02b9f166-e56d-42ad-aa9e-6c7b084860bc
pnpresult = sweep(pnpdata; eneutral = false, tunnel = false)

# ╔═╡ 9d9d5b00-3453-4a62-8b0d-f671140b7a11
let
    fig = Figure()
    ax = Axis(fig[1, 1], yscale = log10)
    T = pnpresult.times
    #lines!(ax,T, voltages.(T))
    lines!(ax, T[2:end], T[2:end] - T[1:(end - 1)])
    fig
end

# ╔═╡ c7d0a14b-b24f-4dfa-90d9-65f06d86ed50
nnpresult = sweep(pnpdata; eneutral = true, tunnel = false)

# ╔═╡ 5a692843-9dc7-4104-a0aa-cff229b1fabb
let
    fig = Figure(size = (650, 400))
    ax = Axis(fig[1, 1])
    lines!(
        ax, voltages(pnpresult), currents(pnpresult, iO),
        color = RGBf.(range(0.1, 1, length(voltages(pnpresult))), 0.0, 0.0)
    )
    lines!(
        ax, voltages(pnpresult), currents(nnpresult, iO),
        color = RGBf.(0.0, range(0.1, 1, length(voltages(nnpresult))), 0.0)
    )
    #ylims!(-0.0001, 0.0001)
    fig
end

# ╔═╡ fe0ee3a5-3825-4ab7-b0f8-b1ca8c698683
floataside(
    md"""
    __Time:__ $(@bind it PlutoUI.Slider(1:length(pnpresult.tsol.t)-1, show_value=false))
    """, top = 380
)


# ╔═╡ 57a31265-0012-43c4-87b1-ff1f2f4d0c29
let
    ZR = zstr(zR)
    ZO = zstr(zO)

    XX = (X[2:end] / nm)
    uu = pnpresult.tsol[:, :, it + 1]
    ru = uu / (mol / dm^3)

    fig = Figure(size = (650, 400))
    ax1 = Axis(
        fig[1, 1];
        #       xlabel=L"x/nm",
        ylabel = L"c/(mol/dm^3)",
        xscale = log10,
        yscale = log10,
        title = "V=$(pnpresult.voltages[it])"
    )
    xlims!(ax1, XX[1], L / nm)
    ylims!(ax1, 1.0e-10, 1.0e2)

    lines!(ax1, XX, ru[iO, 2:end], color = :red, linestyle = :solid, label = L"O^{%$(ZO)}, pnp")
    lines!(ax1, XX, ru[iR, 2:end], color = :blue, linestyle = :solid, label = L"R^{%$(ZR)}, pnp")
    lines!(ax1, XX, ru[iS, 2:end], color = :lightgreen, linestyle = :solid, label = L"S^-, pnp")
    lines!(ax1, XX, ru[iX, 2:end], color = :orange, linestyle = :solid, label = L"X^+, pnp")

    Legend(
        fig[1, 2], ax1; labelsize = 10,
        backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.5)

    )

    fig
end

# ╔═╡ Cell order:
# ╠═a70cef7d-2a2f-4155-bdf3-fec9df94c63f
# ╠═43b83ec3-3c85-4e2c-9a7d-e102ebc29619
# ╠═2496ac21-e027-49d2-a78d-964193d808e0
# ╠═43cc60e3-8a4d-4fac-9454-7e599243e058
# ╠═430f4195-7693-4682-8e81-4ef2d218a2f9
# ╠═b218dbfc-c627-489c-8f27-8ebb7a721a44
# ╠═9d9d5b00-3453-4a62-8b0d-f671140b7a11
# ╠═f0745600-be37-49d3-92fe-231bcc2e7013
# ╠═483eefd0-a96d-4c32-a939-2517f08a6b3d
# ╠═76dab1ef-520c-41ea-b939-80e7212b27a9
# ╟─5a692843-9dc7-4104-a0aa-cff229b1fabb
# ╟─57a31265-0012-43c4-87b1-ff1f2f4d0c29
# ╠═02b9f166-e56d-42ad-aa9e-6c7b084860bc
# ╠═c7d0a14b-b24f-4dfa-90d9-65f06d86ed50
# ╠═bb03204a-cb5a-4595-ace6-8b52a9440494
# ╠═fe0ee3a5-3825-4ab7-b0f8-b1ca8c698683
# ╠═e73cfcbb-c3c4-49e2-a127-a505654aa8b4
# ╟─8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
# ╟─baa3e08e-5d64-4c8f-9f6d-5fdb40e97bc5
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╟─afe4745f-f9f1-4e23-8735-cbec6fb79c41
