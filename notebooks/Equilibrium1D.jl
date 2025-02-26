### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 14f8c67a-759a-4646-811c-01d03e3cf726
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "docs")) # hide
    using Revise
    using HypertextLiteral
    using PlutoUI
    using ExtendableGrids
    using GridVisualize
    using NLsolve
    using LiquidElectrolytes
	using LiquidElectrolytes: set_molarity!, apply_voltage!, update_derived!, calc_cmol, calc_QBL, calc_φ, calc_p, ysum
    using VoronoiFVM: VoronoiFVM, unknowns, history, solve, boundary_dirichlet!
    using LessUnitful
    using Colors
    if isdefined(Main,:PlutoRunner)
	using CairoMakie	
	default_plotter!(CairoMakie)
	CairoMakie.activate!(type="png")
    end
end

# ╔═╡ 6c4eac8c-e243-4e0b-9daf-2ef2d35fe873
md"""
# Equilibrium1D.jl
"""

# ╔═╡ 0f2a3eb0-9818-4bf8-9dde-ba28ea3dd2f5
htl""" Description of the equilibrium problem: 
<a href="./open?path=$(joinpath(@__DIR__,\"..\",\"src\",\"equilibrium.jl\"))" target="_blank"> here </a>"""

# ╔═╡ f36552fd-affd-44e0-83b5-3401459f0560
TableOfContents(title="")

# ╔═╡ 5bce5b30-b5fd-4e13-8b20-8218edaa6c60
md"""
Results should coincide with Fuhrmann, CPC 2015
"""

# ╔═╡ b43533f6-a948-418c-8539-2d54aa8e5943
begin
	const V=ufac"V"
	const eV=ufac"eV"
	const nm=ufac"nm"
	const cm=ufac"cm"
    const μF=ufac"μF"
    const iA=1
    const iC=2
end


# ╔═╡ 5391c130-6706-4f06-8335-6474875fe210
md"""
### Geometry
"""

# ╔═╡ 7062de01-4986-45eb-8bcf-73ff23445dc2
Vmax=2*V

# ╔═╡ 6cd423b6-a6a7-485c-8531-99a51e19de14
L=20nm

# ╔═╡ 0254a47d-33cf-4af7-bb72-c584b9bf98b6
hmin=0.001*nm

# ╔═╡ 7a123ca0-d407-4dc0-9c01-f90ac53b690d
hmax=1*nm

# ╔═╡ 2fdef862-536a-4401-8b5c-bc7e80b6a224
X=ExtendableGrids.geomspace(0,L,hmin,hmax)

# ╔═╡ ac75d0a2-df48-48be-af40-dbb7b7670bf2
grid=ExtendableGrids.simplexgrid(X)

# ╔═╡ 634e25da-8636-413c-9f29-b3b0f3418e69
const solvation=10

# ╔═╡ 1263e426-c510-47d4-a0c3-6023cf98c11f
md"""
### Solution for various applied voltages
"""

# ╔═╡ 4f57312e-1b08-43a3-b038-d30bfc62e754
begin
    data=EquilibriumData();
	data.χ=78; 
	data.κ.=solvation; 
	set_molarity!(data,0.01)
	update_derived!(data)
end

# ╔═╡ d0776266-c208-4237-99c8-1364527eaeb1
dlcap0(data)/(μF/cm^2)

# ╔═╡ 681a3a96-28a8-4168-b263-02a360da0ca8
sys_sy=create_equilibrium_system(grid,data)

# ╔═╡ 08fbabbf-51b4-4398-b183-cea0bfeda76e
sys_pp=create_equilibrium_pp_system(grid,data,Γ_bulk=2)

# ╔═╡ 1ab27251-0999-49a7-a970-29e70d8fd800
inival=unknowns(sys_sy,inival=0);

# ╔═╡ 7899d9b8-0edc-443f-a5a9-96a01187ff74
md"""
Change applied voltage. $(@bind voltage PlutoUI.Slider(-Vmax:0.05:Vmax,show_value=true,default=0))
"""

# ╔═╡ 9545c38c-35d4-4adf-bd80-82af84e93564
begin
	apply_voltage!(sys_sy,voltage)
	apply_voltage!(sys_pp,voltage)
	sol_sy=solve(sys_sy,inival=inival,log=true)
	sol_pp=solve(sys_pp,inival=inival,log=true)
	inival.=sol_sy
end;

# ╔═╡ b737de17-593e-40b9-b147-41516225530d
history(sol_sy),history(sol_pp)

# ╔═╡ 11e94e40-fbcf-4d03-ab86-09cc2e75c5c4
cmol=calc_cmol(sol_pp,sys_pp)

# ╔═╡ 10c87e6d-4485-4d17-b7f2-8ead685e17f4
q=calc_QBL(sol_pp,sys_pp)

# ╔═╡ 1765d933-c6a4-4302-b19e-98d27954c0b4
function plot(sol,sys)
	vis=GridVisualizer(resolution=(600,200),legend=:rt,limits=(-60,60))
	
	scalarplot!(vis,grid,calc_φ(sol,sys),xlimits=(0,5*nm),color=:green,clear=true,label="φ/V")
	scalarplot!(vis,grid,cmol[iA,:],color=:blue,clear=false,label="c_A/(mol/L)")
	scalarplot!(vis,grid,cmol[iC,:],color=:red,clear=false,label="c_C/(mol/L)")
	scalarplot!(vis,grid,clear=false,calc_p(sol,sys)/ufac"GPa",color=:magenta,label="p/GPa",
	title="φ_0=$(voltage)",ylabel="")
	reveal(vis)
end;

# ╔═╡ dcec8754-282f-42c5-9a27-39c07f735af7
plot(sol_sy,sys_sy)

# ╔═╡ 0ba63b59-bfb5-4bdf-9590-210fced0684c
ysum(sys_pp,sol_pp)

# ╔═╡ a95d311c-7f31-456b-974f-526948eef169
scalarplot(grid,ysum(sys_pp,sol_pp),xlimits=(0,5nm),show=true,
	title="Molar fraction sum constraint for pressure Poisson at $(voltage)V",size=(600,200))

# ╔═╡ 768bc478-339f-4bb5-8736-e0377d219744
md"""
### Double layer capacitance
"""

# ╔═╡ c20f9bde-5ff5-47fe-b543-758159f8add5
molarities=[0.001,0.01,0.1,1]

# ╔═╡ d176f826-b8ec-4bdf-b75f-55f96e596a34
function capsplot(sys,molarities)
	vis=GridVisualizer(resolution=(500,300),legend=:rt,clear=true,xlabel="φ/V",ylabel="C_dl/(μF/cm^2)")
	hmol=1/length(molarities)
	for imol=1:length(molarities)
		c=RGB(1-imol*hmol,0,imol*hmol)
		t=@elapsed volts,caps=dlcapsweep_equi(sys,vmax=1V,nsteps=100,molarity=molarities[imol])
		cdl0=dlcap0(sys.physics.data)
		@info "elapsed=$(t)"
	    scalarplot!(vis,volts,caps/(μF/cm^2),
			color=c,clear=false,label="$(molarities[imol])M",markershape=:none)
		scalarplot!(vis,[0],[cdl0]/(μF/cm^2),
			clear=false,markershape=:circle,label="")
	end
#	save(plotsdir("1DResults.pdf"),vis)
	reveal(vis)
end

# ╔═╡ c158042f-9782-448d-9387-5c219ab3958f
md"""
#### Calculation via molar fraction sum constraint
"""

# ╔═╡ 97c852a9-a14a-41d6-8c95-3cea2d760f4a
 capsplot(sys_sy,molarities)

# ╔═╡ 00b889de-3a0c-4e90-92ff-ef9e3a9a69d3
md"""
#### Calculation with pressure Poisson
"""

# ╔═╡ 0e01553e-4d78-48b3-aa1c-e4d01d99d2ce
capsplot(sys_pp,molarities)

# ╔═╡ aec30344-5f87-4c48-bc57-20a3cd41feeb
md"""
### Calculation with modified Poisson-Boltzmann
"""

# ╔═╡ 758f5a7c-a928-4fa2-9366-3973fa161b48
function bcondition(f,u,bnode,data)
	(;Γ_we,Γ_bulk,ϕ_we) = data
       	iϕ,ip=1,2 
	## Dirichlet ϕ=ϕ_we at Γ_we
	boundary_dirichlet!(f,u,bnode,species=iϕ,region=Γ_we,value=ϕ_we)
	boundary_dirichlet!(f,u,bnode,species=iϕ,region=Γ_bulk,value=data.ϕ_bulk)
	boundary_dirichlet!(f,u,bnode,species=ip,region=Γ_bulk,value=data.p_bulk)
    end

# ╔═╡ 7e08d884-62ac-40f7-a705-d617cbec3c9e
function create_pbsystem()
    celldata = ElectrolyteData(;
                               nc = 2,
                               Γ_we = 1,
                               Γ_bulk = 2)
	
	celldata.κ.=solvation
	pbsys=PBSystem(grid;celldata,bcondition)
end


# ╔═╡ 9bda8daa-91f8-40fe-b467-2fd6c9fcf0a2
pbsys=create_pbsystem()

# ╔═╡ 75532ae4-15c0-4a7f-8694-1b419e0657ec
result=dlcapsweep(pbsys;inival=unknowns(pbsys,inival=0.0),iϕ=1,voltages=-1:0.01:1)

# ╔═╡ 72118f4d-95ae-4b8d-9a7a-af36190ea521
function pbcapsplot(sys,molarities)
	vis=GridVisualizer(resolution=(500,300),legend=:rt,clear=true,xlabel="φ/V",ylabel="C_dl/(μF/cm^2)")
	hmol=1/length(molarities)
	voltages=range(-1.0,1,length=101)
	for imol=1:length(molarities)
		c=RGB(1-imol*hmol,0,imol*hmol)
		t=@elapsed result=dlcapsweep(pbsys;inival=unknowns(pbsys,inival=0),iϕ=1,voltages,molarity=molarities[imol]*ufac"mol/dm^3")
		(;voltages,dlcaps)=result
		cdl0=dlcap0(sys.physics.data)
		@info "elapsed=$(t)"
	    scalarplot!(vis,voltages,dlcaps/(μF/cm^2),
			color=c,clear=false,label="$(molarities[imol])M",markershape=:none)
		scalarplot!(vis,[0],[cdl0]/(μF/cm^2),
			clear=false,markershape=:circle,label="")
	end
#	save(plotsdir("1DResults.pdf"),vis)
	reveal(vis)
end

# ╔═╡ f4493da0-3f79-45ae-995f-7222ce4ba892
pbcapsplot(pbsys,molarities)


# ╔═╡ Cell order:
# ╟─6c4eac8c-e243-4e0b-9daf-2ef2d35fe873
# ╟─0f2a3eb0-9818-4bf8-9dde-ba28ea3dd2f5
# ╠═14f8c67a-759a-4646-811c-01d03e3cf726
# ╟─f36552fd-affd-44e0-83b5-3401459f0560
# ╟─5bce5b30-b5fd-4e13-8b20-8218edaa6c60
# ╠═b43533f6-a948-418c-8539-2d54aa8e5943
# ╟─5391c130-6706-4f06-8335-6474875fe210
# ╠═7062de01-4986-45eb-8bcf-73ff23445dc2
# ╠═6cd423b6-a6a7-485c-8531-99a51e19de14
# ╠═0254a47d-33cf-4af7-bb72-c584b9bf98b6
# ╠═7a123ca0-d407-4dc0-9c01-f90ac53b690d
# ╠═2fdef862-536a-4401-8b5c-bc7e80b6a224
# ╠═ac75d0a2-df48-48be-af40-dbb7b7670bf2
# ╠═634e25da-8636-413c-9f29-b3b0f3418e69
# ╟─1263e426-c510-47d4-a0c3-6023cf98c11f
# ╠═4f57312e-1b08-43a3-b038-d30bfc62e754
# ╠═d0776266-c208-4237-99c8-1364527eaeb1
# ╠═681a3a96-28a8-4168-b263-02a360da0ca8
# ╠═08fbabbf-51b4-4398-b183-cea0bfeda76e
# ╠═1ab27251-0999-49a7-a970-29e70d8fd800
# ╠═9545c38c-35d4-4adf-bd80-82af84e93564
# ╠═b737de17-593e-40b9-b147-41516225530d
# ╠═11e94e40-fbcf-4d03-ab86-09cc2e75c5c4
# ╠═10c87e6d-4485-4d17-b7f2-8ead685e17f4
# ╟─7899d9b8-0edc-443f-a5a9-96a01187ff74
# ╟─1765d933-c6a4-4302-b19e-98d27954c0b4
# ╟─dcec8754-282f-42c5-9a27-39c07f735af7
# ╠═0ba63b59-bfb5-4bdf-9590-210fced0684c
# ╠═a95d311c-7f31-456b-974f-526948eef169
# ╟─768bc478-339f-4bb5-8736-e0377d219744
# ╠═c20f9bde-5ff5-47fe-b543-758159f8add5
# ╠═d176f826-b8ec-4bdf-b75f-55f96e596a34
# ╟─c158042f-9782-448d-9387-5c219ab3958f
# ╠═97c852a9-a14a-41d6-8c95-3cea2d760f4a
# ╟─00b889de-3a0c-4e90-92ff-ef9e3a9a69d3
# ╠═0e01553e-4d78-48b3-aa1c-e4d01d99d2ce
# ╟─aec30344-5f87-4c48-bc57-20a3cd41feeb
# ╠═758f5a7c-a928-4fa2-9366-3973fa161b48
# ╠═7e08d884-62ac-40f7-a705-d617cbec3c9e
# ╠═9bda8daa-91f8-40fe-b467-2fd6c9fcf0a2
# ╠═75532ae4-15c0-4a7f-8694-1b419e0657ec
# ╠═72118f4d-95ae-4b8d-9a7a-af36190ea521
# ╠═f4493da0-3f79-45ae-995f-7222ce4ba892
