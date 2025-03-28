### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 4a215269-da12-431b-9e56-beec46112998
# ╠═╡ skip_as_script = true
#=╠═╡
begin
   using Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "docs")) # hide
    using Revise
	    using GridVisualize

   	import CairoMakie	
	    using CairoMakie: lines!
  	 	default_plotter!(CairoMakie)
 		CairoMakie.activate!(type="png")
  
end
  ╠═╡ =#

# ╔═╡ a70cef7d-2a2f-4155-bdf3-fec9df94c63f
begin
    using LiquidElectrolytes
    using Printf
    using ExtendableGrids
    using LinearAlgebra
    using PlutoUI
    using ExtendableFEM
    using VoronoiFVM
    using LessUnitful
    using Test

end

# ╔═╡ baa3e08e-5d64-4c8f-9f6d-5fdb40e97bc5
begin
    using HypertextLiteral: @htl_str, @htl
    using UUIDs: uuid1
end

# ╔═╡ 8db5d72c-34a9-411e-ade5-5f48b99f9231
begin
    @phconstants e N_A
    const F = e * N_A
    @unitfactors mol m dm cm nm μA s mV V mPa

end

# ╔═╡ 8af5f374-d379-44a6-b567-ca34f6a2149e
begin
    const L = 50 * nm
    const W = 5 * nm
    const nref = 0
    const Δp = 1.0
    const μ = 0.89 * mPa * s
    const cyl = false
    nr = 8 * 2^nref + 1
    nl = Int((nr - 1) * L / W) + 1
    Y = range(0, L, length = nl)
    flowX = range(0, W, length = nr)
    pnpX = geomspace(0, W, 2 * W / nr, 0.25 * W / nr)
    pnpgrid = simplexgrid(pnpX, Y)
    flowgrid = simplexgrid(flowX, Y)
    flowgrid = pnpgrid
    xgrid = simplexgrid(1.0e9 * 2.5 * pnpX, 1.0e9 * Y)
    pnpgrid[CellEdges]
    pnpgrid
end

# ╔═╡ f49f1197-257b-4fdb-ad80-5236c3aa1ee7
#=╠═╡
let
    if isdefined(Main,:PlutoRunner)
        vis = GridVisualizer(layout = (1, 2), size = (650, 300), linewidth = 0.1)
        gridplot!(vis[1, 1], flowgrid, title = "flowgrid", aspect = 0.5)
        gridplot!(vis[1, 2], pnpgrid, title = "pnpgrid", aspect = 0.5)
        reveal(vis)
    end
end
  ╠═╡ =#

# ╔═╡ 154c142a-c2b2-4229-80ec-9d7cf924c04c
begin
    const cbulk = 1.0 * mol / dm^3
    const Δϕ = 0.5 * V
    const σ = -15 * μA * s / cm^2
    const σ_embed = [σ]
    embed(data, λ) = σ_embed[1] = σ * λ
end

# ╔═╡ 21d159e0-1639-451e-9054-2729a5942b5f
function pnpbcond(y, u, bnode, data)
    (; iϕ, nc) = data
    boundary_neumann!(y, u, bnode, species = iϕ, region = 2, value = σ_embed[1])
    boundary_dirichlet!(y, u, bnode, species = iϕ, region = 1, value = -Δϕ / 2)
    boundary_dirichlet!(y, u, bnode, species = iϕ, region = 3, value = Δϕ / 2)
    for i in 1:nc
        boundary_dirichlet!(y, u, bnode, species = i, region = 1, value = cbulk)
        boundary_dirichlet!(y, u, bnode, species = i, region = 3, value = cbulk)
    end
    return
end

# ╔═╡ 3c6bac8c-1286-4ef2-8462-5c055fe38e9d
function flowbcond(flowsolver)
    assign_operator!(
        flowsolver,
        HomogeneousBoundaryData(
            LiquidElectrolytes.velocity_unknown(flowsolver);
            regions = [4], mask = (1, 0, 1)
        )
    )
    return assign_operator!(
        flowsolver,
        HomogeneousBoundaryData(LiquidElectrolytes.velocity_unknown(flowsolver); regions = [2])
    )
end

# ╔═╡ 280f5233-621e-4d9d-95f2-8c97a2c343fc
default_pnpdata = ElectrolyteData();

# ╔═╡ eb84642a-a26d-4b26-8530-35a499a179a8
function PNPData(; molar_volume = default_pnpdata.v0, scheme = :μex, κ = 20)
    pnpdata = ElectrolyteData()
    pnpdata.edgevelocity = zeros(num_edges(pnpgrid))
    pnpdata.solvepressure = false
    pnpdata.scheme = scheme
    pnpdata.κ .= κ
    pnpdata.D .= 1.0e-9 * m^2 / s
    pnpdata.pscale = 1
    pnpdata.v .= molar_volume
    return pnpdata
end

# ╔═╡ 9aaa7a1f-23ea-466b-9d54-6193879ae9fa
gcdata = PNPData(molar_volume = 0, κ = 0)

# ╔═╡ f78d5edc-9e67-4d93-87f2-a08a9fc3b845
dgldata = PNPData()

# ╔═╡ d176f9d9-4f50-4d99-97b7-bbe9d1dd1396
begin
    gcsolver = LiquidElectrolytes.PNPStokesSolver(;
        flowgrid,
        pnpgrid,
        μ,
        velospace = H1BR,
        pnpdata = gcdata,
        pnpbcond,
        pnpreaction = (y, u, bnode, data) -> nothing,
        flowbcond
    )
end;

# ╔═╡ 11006383-62e1-4491-a1c9-a28b99407c0b
begin
    dglsolver = LiquidElectrolytes.PNPStokesSolver(;
        flowgrid,
        pnpgrid,
        μ,
        velospace = H1BR,
        pnpdata = dgldata,
        pnpbcond,
        pnpreaction = (y, u, bnode, data) -> nothing,
        flowbcond
    )
end;

# ╔═╡ bb157989-88a1-42f8-b292-c4d1c436db03
gcpnpsol, gcflowsol = solve(
    gcsolver;
    verbose = "n",
    niter = 20,
    damp0 = 0.1,
    damp_initial = 1,
    tol_round = 1.0e-9,
    max_round = 3,
    embed,
    nembed = 10
)

# ╔═╡ 49a89a4b-d8ab-4bdb-9603-45dce344fc00
dglpnpsol, dglflowsol = solve(
    dglsolver;
    verbose = "n",
    niter = 20,
    damp0 = 0.1,
    damp_initial = 1,
    tol_round = 1.0e-9,
    max_round = 3,
    embed,
    nembed = 10
)

# ╔═╡ 666d12aa-5fcb-49d5-a2ce-b245274af693
#=╠═╡
if isdefined(Main,:PlutoRunner)
    LiquidElectrolytes.flowplot(gcflowsol, gcsolver.flowsolver;
                                Plotter = CairoMakie,
                                vscale = 0.5,
                                aspect = 0.25,
                                rasterpoints = 30)
end
  ╠═╡ =#

# ╔═╡ 804642e9-a290-4c37-b500-3e11345c1ad5
#=╠═╡
if isdefined(Main,:PlutoRunner)
    LiquidElectrolytes.flowplot(dglflowsol, dglsolver.flowsolver;
                                Plotter = CairoMakie,
                                vscale = 0.5,
                                aspect = 0.25,
                                rasterpoints = 30)
end
  ╠═╡ =#

# ╔═╡ 00d1e5e7-f53f-4aa4-b06f-d283b9cde21a
#=╠═╡
let
    if isdefined(Main,:PlutoRunner)
        vis = GridVisualizer(; layout = (2, 2), size = (650, 600))
        scalarplot!(vis[1, 1], xgrid, gcpnpsol[1, :] / cbulk)
        scalarplot!(vis[1, 2], xgrid, gcpnpsol[2, :] / cbulk)
        scalarplot!(vis[2, 1], xgrid, gcpnpsol[3, :], colormap = :bwr, limits = (-1, 1))
        scalarplot!(vis[2, 2], xgrid, gcpnpsol[4, :])
        reveal(vis)
    end
end
  ╠═╡ =#

# ╔═╡ a91c5744-e744-4c25-a748-30d7155e9b63
#=╠═╡
let
    if isdefined(Main,:PlutoRunner)
        vis = GridVisualizer(; layout = (2, 2), size = (650, 600))
        scalarplot!(vis[1, 1], xgrid, dglpnpsol[1, :] / cbulk)
        scalarplot!(vis[1, 2], xgrid, dglpnpsol[2, :] / cbulk)
        scalarplot!(vis[2, 1], xgrid, dglpnpsol[3, :], colormap = :bwr, limits = (-1, 1))
        scalarplot!(vis[2, 2], xgrid, dglpnpsol[4, :])
        reveal(vis)
    end
end
  ╠═╡ =#

# ╔═╡ fe4a9591-0842-42ac-9c77-76f6125d6688
function midvelo(flowsol, pnpsol, solver)
    nodevelo = LiquidElectrolytes.node_velocity(flowsol, solver.flowsolver)
    nx = length(pnpX)
    iymid = (length(Y) - 1) ÷ 2 + 1
    midrange = ((iymid - 1) * nx + 1):(iymid * nx)
    c_p = pnpsol[1, midrange]
    c_m = pnpsol[2, midrange]
    ϕ = -pnpsol[3, midrange]
    p = pnpsol[4, midrange]
    u = nodevelo[2, midrange]
    return c_p, c_m, ϕ, p, u
end

# ╔═╡ 98a70371-71b0-4416-8395-6a03eb6873c0
#=╠═╡
if isdefined(Main, :PlutoRunner)
	fig=CairoMakie.Figure(size=(650,300))
	c_p, c_m, ϕ, p, u = midvelo(gcflowsol,gcpnpsol, gcsolver)
	ax1=CairoMakie.Axis(fig[1,1],
		yaxisposition=:left,
		xticks=[0,1,2,3,4,5],
		ylabel="u/(m/s)",
		xlabel="x/nm",
		title="Zero ion size"
		)
	CairoMakie.ylims!(ax1,0,0.8)
	CairoMakie.xlims!(ax1,0,5)
	ax2=CairoMakie.Axis(fig[1,1],
		yaxisposition=:right,
		ylabel="c/(mol/L); ϕ/(10mV)"
		)
	CairoMakie.ylims!(ax2,0,10)
	CairoMakie.xlims!(ax2,0,5)

	data=[
lines!(ax1, pnpX / nm, u, color=:green,linewidth = 2)
 lines!(ax2, pnpX / nm, c_p / (mol / dm^3), color=:red, linewidth = 2)
 lines!(ax2, pnpX / nm, c_m / (mol / dm^3), color=:blue, linewidth = 2)
lines!(ax2, pnpX / nm, ϕ / (10mV), color=:black, linewidth = 2)
	]
CairoMakie.axislegend(ax1, data, ["u_z","c^+","c^-","Φ"], position=:ct, backgroundcolor = :transparent)

	
		c_p, c_m, ϕ, p, u = midvelo(dglflowsol,dglpnpsol, dglsolver)
	ax1=CairoMakie.Axis(fig[1,2],
		yaxisposition=:left,
		xticks=[0,1,2,3,4,5],
		ylabel="u/(m/s)",
		xlabel="x/nm",
		title="Finite ion size\n Solvation"
		)
	CairoMakie.ylims!(ax1,0,0.8)
	CairoMakie.xlims!(ax1,0,5)
	ax2=CairoMakie.Axis(fig[1,2],
		yaxisposition=:right,
		ylabel="c/(mol/L); ϕ/(10mV)"
		)
	CairoMakie.ylims!(ax2,0,10)
	CairoMakie.xlims!(ax2,0,5)

	data=[
lines!(ax1, pnpX / nm, u, color=:green, linewidth = 2),
 lines!(ax2, pnpX / nm, c_p / (mol / dm^3), color=:red,linewidth = 2),
 lines!(ax2, pnpX / nm, c_m / (mol / dm^3), color=:blue, linewidth = 2),
lines!(ax2, pnpX / nm, ϕ / (10mV), color=:black, linewidth = 2)
	]
CairoMakie.axislegend(ax1, data, ["u_z","c^+","c^-","Φ"], position=:ct, backgroundcolor = :transparent)
	fig
	
end
  ╠═╡ =#

# ╔═╡ 43b89443-58ad-4321-8db7-ead924a4b17b
begin
    gc_c_p, gc_c_m, gc_ϕ, gc_p, gc_u = midvelo(gcflowsol, gcpnpsol, gcsolver)
    dgl_c_p, dgl_c_m, dgl_ϕ, dgl_p, dgl_u = midvelo(dglflowsol, dglpnpsol, dglsolver)
    @test gc_u[1] < dgl_u[1]
    @test gc_c_p[end] > dgl_c_p[end]
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


# ╔═╡ Cell order:
# ╠═4a215269-da12-431b-9e56-beec46112998
# ╠═a70cef7d-2a2f-4155-bdf3-fec9df94c63f
# ╠═8db5d72c-34a9-411e-ade5-5f48b99f9231
# ╠═8af5f374-d379-44a6-b567-ca34f6a2149e
# ╠═f49f1197-257b-4fdb-ad80-5236c3aa1ee7
# ╠═154c142a-c2b2-4229-80ec-9d7cf924c04c
# ╠═21d159e0-1639-451e-9054-2729a5942b5f
# ╠═3c6bac8c-1286-4ef2-8462-5c055fe38e9d
# ╠═280f5233-621e-4d9d-95f2-8c97a2c343fc
# ╠═eb84642a-a26d-4b26-8530-35a499a179a8
# ╠═9aaa7a1f-23ea-466b-9d54-6193879ae9fa
# ╠═f78d5edc-9e67-4d93-87f2-a08a9fc3b845
# ╠═d176f9d9-4f50-4d99-97b7-bbe9d1dd1396
# ╠═11006383-62e1-4491-a1c9-a28b99407c0b
# ╠═bb157989-88a1-42f8-b292-c4d1c436db03
# ╠═49a89a4b-d8ab-4bdb-9603-45dce344fc00
# ╠═666d12aa-5fcb-49d5-a2ce-b245274af693
# ╠═804642e9-a290-4c37-b500-3e11345c1ad5
# ╠═00d1e5e7-f53f-4aa4-b06f-d283b9cde21a
# ╠═a91c5744-e744-4c25-a748-30d7155e9b63
# ╠═fe4a9591-0842-42ac-9c77-76f6125d6688
# ╠═98a70371-71b0-4416-8395-6a03eb6873c0
# ╠═43b89443-58ad-4321-8db7-ead924a4b17b
# ╟─8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
# ╟─baa3e08e-5d64-4c8f-9f6d-5fdb40e97bc5
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╟─afe4745f-f9f1-4e23-8735-cbec6fb79c41
