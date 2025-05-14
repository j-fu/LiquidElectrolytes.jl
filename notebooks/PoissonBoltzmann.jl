### A Pluto.jl notebook ###
# v0.20.8

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

# ╔═╡ 18ce65e8-078d-4dec-bc4d-a05ce77acfd3
# ╠═╡ skip_as_script = true
#=╠═╡
begin
    import Pkg # hide
    Pkg.activate(joinpath(@__DIR__, "..", "docs")) # hide
    using Revise # hide
	using GridVisualize
		using CairoMakie	
		default_plotter!(CairoMakie)
		CairoMakie.activate!(type="svg")
end
  ╠═╡ =#

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin
    using PlutoUI, HypertextLiteral, UUIDs
    using LinearAlgebra
    using Interpolations
    using VoronoiFVM, ExtendableGrids
    using LiquidElectrolytes
    using LessUnitful, Unitful
end

# ╔═╡ 25a19579-9cba-4416-87cb-83d00e5926f3
md"""
# PoissonBoltzmann.jl
"""

# ╔═╡ 1afa972d-74cf-47ce-86cf-8751a9e85218
md"""
Under development, may be removed, though. Same as VoronoiFVMPoisson, but with the intention do the same things with LiquidElectrolytes.jl
"""

# ╔═╡ b181a82b-46e7-4dae-b0b1-6886dd40886a
md"""
## Constants, Units & Parameters
"""

# ╔═╡ 5a5e3586-2009-4381-ae9d-a0bd34539b65
md"""
See LessUnitful.jl  `ph` and `ufac` and Unitful.jl for `u`
"""

# ╔═╡ b6338637-1f3c-4410-8e82-3afdaf603656
begin
    const ε = 78.49
    const ε_0 = ph"ε_0"
    const F = ph"N_A" * ph"e"
    const K = ufac"K"
    const nm = ufac"nm"
    const dm = ufac"dm"
    const V = ufac"V"
    const mol = ufac"mol"
end;

# ╔═╡ 371df98d-9f55-4244-8243-e6b597cc8d84
md"""
Parameters:
"""

# ╔═╡ 10be79ba-9ab9-4e6f-8568-2de93ea77964
begin
    const molarity = 0.01
    const c_bulk = [molarity, molarity] * mol / dm^3 # bulk ion molar concentrations
    const z = [1, -1] # species charge numbers
    const c̄ = 55.508mol / dm^3 # solvent (water) molar concentration
    const nref = 0 # grid refinement level
    const T = 298.15K # temperature
    const L = 20.0nm # computational domain size
    const hmin = 1.0e-1 * nm * 2.0^(-nref) # grid size at working electrode
    const hmax = 1.0 * nm * 2.0^(-nref) # grid size at bulk

end;

# ╔═╡ 7da889ce-9c6b-4abc-b19d-6311aacc32b1
X = geomspace(0, L, hmin, hmax)

# ╔═╡ a2a8132a-d54d-48a7-aa02-ff83593f0e16
grid = simplexgrid(X)

# ╔═╡ db023dac-5ce0-4126-bfd0-9036ccec11f4
#=╠═╡
gridplot(grid,size=(600,200))
  ╠═╡ =#

# ╔═╡ 834c81c8-9035-449f-b331-73bbee87756c
md"""
Check for bulk electroneutrality:
"""

# ╔═╡ a4a4c6a4-ea90-4d23-b13a-2020790b2889
md"""
## Classical Poisson-Boltzmann

```math
  -\nabla \cdot εε_0 \nabla \phi = F\sum_{i=1}^n z_ic_0 \exp\left(z_i\phi \frac{F}{RT}\right)
```
"""

# ╔═╡ 41715397-020c-4505-a61a-4f2910318423
begin
	function  pbo_gamma(γ, c,p, electrolyte)
     	for i=1:electrolyte.nc
		   γ[i]=1
    	end
	    return nothing
    end
    pbo_electrolyte = ElectrolyteData(;γ! = pbo_gamma, c_bulk)
 end

# ╔═╡ 7f209979-0776-40ab-9b2b-3b2d0145fa2c
iselectroneutral(c_bulk, pbo_electrolyte)

# ╔═╡ 80a42b07-d6b1-41d6-bf1d-66a01eeb5019
md"""
Derived data:
"""

# ╔═╡ c57e82bb-79d8-4553-b31d-d1448c38c649
md"""
Debye length= $(debyelength(pbo_electrolyte) |> x->round(x,sigdigits=5) |>u"nm")
"""

# ╔═╡ e2179a6c-cc1c-4850-89d4-d3f87fd5e6ee
md"""
Double layer capacitance at zero voltage for symmetric binary electrolyte = 
$(dlcap0(pbo_electrolyte) |> x->round(x,sigdigits=5) |>u"μF/cm^2")
"""

# ╔═╡ 66a988e0-9388-4bbe-8d92-9177dfca8ac4
function pb_bcondition(f, u, bnode, data)
    (; Γ_we, Γ_bulk, ϕ_we, iϕ, ip) = data
    ## Dirichlet ϕ=ϕ_we at Γ_we
    boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = ϕ_we)
    boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_bulk, value = data.ϕ_bulk)
    return boundary_dirichlet!(f, u, bnode, species = ip, region = Γ_bulk, value = data.p_bulk)

end


# ╔═╡ 91a661ee-9c0c-496e-b6fb-4bb692cb7255
pbo_system = PBSystem(grid; celldata = pbo_electrolyte, bcondition = pb_bcondition)

# ╔═╡ 240b9017-bb20-4e0f-b372-3248a1cf0a9f
pbo_result = dlcapsweep(
    pbo_system; voltages = range(0, 1, length = 101),
    store_solutions = true
)


# ╔═╡ 7ce95a6a-dd5f-4c2a-810b-6483f02ce661
pbo_volts = voltages(pbo_result)

# ╔═╡ 7a014d29-535b-4389-987c-cbfcde5fd8a4
pbo_caps = pbo_result.dlcaps

# ╔═╡ c4af7d22-52ad-495e-b1db-0cd129173c60
pbo_sols = voltages_solutions(pbo_result);

# ╔═╡ ad804cc3-93cd-421b-8314-a7d11e051e3d
@bind pbo_v PlutoUI.Slider(range(pbo_volts[begin], pbo_volts[end], length = 201), default = 0.05, show_value = true)

# ╔═╡ db5b6820-5a53-465c-b380-66756fd722a6
md"""
## "Poisson-Bikerman"
"""

# ╔═╡ dbc500d6-99e8-40d8-9b60-dcccf1889ce8
begin
	function  pbi_gamma(γ, c,p, electrolyte)
		sumc=zero(eltype(c))
		for i=1:electrolyte.nc
			sumc+=c[i]
		end
	    g= 1.0/(1.0 - sumc*electrolyte.v0)
     	for i=1:electrolyte.nc
		   γ[i]=g
    	end
	    return nothing
    end
   
    pbi_electrolyte = ElectrolyteData(c_bulk=pbo_electrolyte.c_bulk, 
									  γ! =pbi_gamma
									 )
	pbi_electrolyte.v .= pbi_electrolyte.v0
    pbi_electrolyte
end

# ╔═╡ 5ccc639b-0fe9-4bd0-a672-4a9d8910d0aa
pbi_system = PBSystem(grid; celldata = pbi_electrolyte, bcondition = pb_bcondition)

# ╔═╡ 7bcad22f-4ad8-4f54-b920-4e5fbe3f2104
pbi_result = dlcapsweep(
    pbi_system; voltages = range(0, 1, length = 101),
    verbose = "",
    store_solutions = true,
	damp_initial=0.1
)


# ╔═╡ 00ea2ca4-5051-4f38-aeb3-122a9664ba39
pbi_volts = voltages(pbi_result)

# ╔═╡ 930aa7cc-47da-434a-b185-98283b3baa0d
pbi_caps = pbi_result.dlcaps

# ╔═╡ f45d7707-1626-42e7-8211-8cdcc4cccbea
pbi_sols = voltages_solutions(pbi_result);

# ╔═╡ 26266ce9-1fd5-4dec-9ea7-2a79c724685d
@bind pbi_v PlutoUI.Slider(range(pbi_volts[begin], pbi_volts[end], length = 201), default = 0.05, show_value = true)

# ╔═╡ bd6749ab-8b75-4a45-8ac4-0650ba4a903c
md"""
## "Poisson-Bikerman" #2 
"""

# ╔═╡ b31af4a2-716c-4c17-b283-0103630fdf38
begin
   
    pbi2_electrolyte = ElectrolyteData(c_bulk=pbo_electrolyte.c_bulk)
	pbi2_electrolyte.v .= pbi_electrolyte.v0
	pbi2_electrolyte.κ .= 0
	pbi2_electrolyte.M .= pbi_electrolyte.M0
    update_derived!(pbi2_electrolyte)
end

# ╔═╡ 56016697-5ce9-4e0a-b14a-bb5d276a358c
pbi2_system = PBSystem(grid; celldata = pbi2_electrolyte, bcondition = pb_bcondition)

# ╔═╡ d261c10c-6003-45b9-91ce-7b7f9acbdc8e
pbi2_result = dlcapsweep(
    pbi2_system; voltages = range(0, 1, length = 101),
    verbose = "",
    store_solutions = true,
	damp_initial=0.1
)


# ╔═╡ d99c6463-9ae2-479a-8c6f-ef59af17fffa
pbi2_volts = voltages(pbi_result)

# ╔═╡ cb8b5bdf-148e-4834-b423-79c221a43f7b
pbi2_caps = pbi2_result.dlcaps

# ╔═╡ 29140461-4e06-4435-87ff-81061dc80d4c
@bind pbi2_v PlutoUI.Slider(range(pbi2_volts[begin], pbi2_volts[end], length = 201), default = 0.05, show_value = true)

# ╔═╡ d64cf39a-5f91-4d3f-ad82-1aef00a0a290
pbi2_sols = voltages_solutions(pbi2_result);

# ╔═╡ e373824d-6e0d-44c8-8c8a-192cccf10298
#=╠═╡
function plotcdl(pbvolts,pbcaps,pbv)
	vc=linear_interpolation(pbvolts,pbcaps)
	cdl=round(vc(pbv)/ufac"μF/cm^2",sigdigits=4)
	vis=GridVisualizer(size=(500,200),xlabel="Δϕ",ylabel="C_dl",legend=:lt)
	scalarplot!(vis,pbvolts,pbcaps/ufac"μF/cm^2")
	scalarplot!(vis,[pbv],[cdl],markershape=:circle,clear=false,markersize=20,
	label="$(cdl)")
	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ 811fca9b-bac0-4003-a0f9-bfefcbfbfa30
#=╠═╡
plotcdl(pbo_volts,pbo_caps,pbo_v)
  ╠═╡ =#

# ╔═╡ eb7d0765-d5e7-4ef9-916d-764c5aca9822
#=╠═╡
plotcdl(pbi_volts,pbi_caps,pbi_v)
  ╠═╡ =#

# ╔═╡ fc55ef51-b7f4-4640-b8e5-ad7b4f8a411c
#=╠═╡
plotcdl(pbi2_volts,pbi2_caps,pbi2_v)
  ╠═╡ =#

# ╔═╡ e480b378-7afa-4906-9e8a-70eac7712b5e
#=╠═╡
function plotsols(pbsols,pbv)
    vis=GridVisualizer(size=(700,200),layout=(1,2),xlabel="x/nm")
	sol=pbsols(pbv)
	ϕ=sol[3,:]
	cp=sol[1,:]/ufac"mol/dm^3"
	cm=sol[2,:]/ufac"mol/dm^3"
	scalarplot!(vis[1,1],X/nm,ϕ,color=:green,ylabel="ϕ/V")
	scalarplot!(vis[1,2],X/nm,cp,color=:red,ylabel="c/ (mol/L)")
	scalarplot!(vis[1,2],X/nm,cm,color=:blue,clear=false)
	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ 65a25bd3-c420-4da3-b808-0e1888c6b4ba
#=╠═╡
plotsols(pbo_sols,pbo_v)
  ╠═╡ =#

# ╔═╡ 5d15d24d-6317-4ad1-a53e-5af5c5bcf28a
#=╠═╡
plotsols(pbi_sols,pbi_v)
  ╠═╡ =#

# ╔═╡ f29de74d-66dd-4d95-9d6a-c33bfb8361c9
#=╠═╡
plotsols(pbi2_sols,pbi2_v)
  ╠═╡ =#

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
        	/* Headers */
             h1{background-color:#dddddd;  padding: 10px;}
             h2{background-color:#e7e7e7;  padding: 10px;}
             h3{background-color:#eeeeee;  padding: 10px;}
             h4{background-color:#f7f7f7;  padding: 10px;}

    		/* "Terminal"  */
    	     pluto-log-dot-sizer  { max-width: 655px;}
             pluto-log-dot.Stdout { background: #002000;
    	                            color: #10f080;
                                    border: 6px solid #b7b7b7;
                                    min-width: 18em;
                                    max-height: 300px;
                                    width: 675px;
                                    overflow: auto;
     	                           }
            /* Standard cell width etc*/
    		main {
       			flex: 1;
    		    max-width: calc(700px + 25px + 6px); /* 700px + both paddings */
        		padding-top: 0px;
        		padding-bottom: 4rem;
        		padding-left: 25px;
        		padding-right: 6px;
        		align-content: center;
        		width: 100%;
    			}

           /* Cell width for slides*/
    		xmain {
    			margin: 0 auto;
    			max-width: 750px;
        		padding-left: max(20px, 3%);
        		padding-right: max(20px, 3%);
    	        }

    	
    	
        </style>
    """
end

# ╔═╡ 9b5f389f-b105-4610-bab7-f79305fedc31
hrule()

# ╔═╡ 5beb3a0d-e57a-4aea-b7a0-59b8ce9ff5ce
hrule()

# ╔═╡ 86e1fffe-cf70-4d93-8131-153f0e98c94b
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

# ╔═╡ 4652157f-3288-45b6-877e-17ed8dc7d18e
function slidemodeswitch()
    uuid = uuid1()
    return html"""
    	<script id=$(uuid)>

        const right = document.querySelector('button.changeslide.next')
        const left = document.querySelector('button.changeslide.prev')

        let fullScreen = false

        const func = (e) => {
            if (e.key == "F10") {
                e.preventDefault()
                window.present()
                if (fullScreen) {
                    document.exitFullscreen().then(() => fullScreen = false)
                } else {

                    document.documentElement.requestFullscreen().then(() => fullScreen = true)
                }
            }
            if (document.body.classList.contains('presentation')) {
             
                if (e.target.tagName == "TEXTAREA") return
                if (e.key == "PageUp") {
                    e.preventDefault()
                    left.click()
                    return
             }

                if (e.key == "PageDown") {
                    e.preventDefault()
                    right.click()
                    return
                }
                if (e.key == "Escape") {
                    window.present()
                    fullScreen = false
                document.exitFullscreen().catch(() => {return})
                }
            }
        }

        document.addEventListener('keydown',func)

        invalidation.then(() => {document.removeEventListener('keydown',func)})
    </script>
    """
end;


# ╔═╡ Cell order:
# ╟─25a19579-9cba-4416-87cb-83d00e5926f3
# ╟─1afa972d-74cf-47ce-86cf-8751a9e85218
# ╠═18ce65e8-078d-4dec-bc4d-a05ce77acfd3
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─b181a82b-46e7-4dae-b0b1-6886dd40886a
# ╟─5a5e3586-2009-4381-ae9d-a0bd34539b65
# ╠═b6338637-1f3c-4410-8e82-3afdaf603656
# ╟─371df98d-9f55-4244-8243-e6b597cc8d84
# ╠═10be79ba-9ab9-4e6f-8568-2de93ea77964
# ╠═7da889ce-9c6b-4abc-b19d-6311aacc32b1
# ╠═a2a8132a-d54d-48a7-aa02-ff83593f0e16
# ╠═db023dac-5ce0-4126-bfd0-9036ccec11f4
# ╟─834c81c8-9035-449f-b331-73bbee87756c
# ╟─a4a4c6a4-ea90-4d23-b13a-2020790b2889
# ╠═41715397-020c-4505-a61a-4f2910318423
# ╠═7f209979-0776-40ab-9b2b-3b2d0145fa2c
# ╟─80a42b07-d6b1-41d6-bf1d-66a01eeb5019
# ╟─c57e82bb-79d8-4553-b31d-d1448c38c649
# ╟─e2179a6c-cc1c-4850-89d4-d3f87fd5e6ee
# ╠═66a988e0-9388-4bbe-8d92-9177dfca8ac4
# ╠═91a661ee-9c0c-496e-b6fb-4bb692cb7255
# ╠═240b9017-bb20-4e0f-b372-3248a1cf0a9f
# ╠═7ce95a6a-dd5f-4c2a-810b-6483f02ce661
# ╠═7a014d29-535b-4389-987c-cbfcde5fd8a4
# ╠═c4af7d22-52ad-495e-b1db-0cd129173c60
# ╟─ad804cc3-93cd-421b-8314-a7d11e051e3d
# ╠═811fca9b-bac0-4003-a0f9-bfefcbfbfa30
# ╠═65a25bd3-c420-4da3-b808-0e1888c6b4ba
# ╟─db5b6820-5a53-465c-b380-66756fd722a6
# ╠═dbc500d6-99e8-40d8-9b60-dcccf1889ce8
# ╠═5ccc639b-0fe9-4bd0-a672-4a9d8910d0aa
# ╠═7bcad22f-4ad8-4f54-b920-4e5fbe3f2104
# ╠═00ea2ca4-5051-4f38-aeb3-122a9664ba39
# ╠═930aa7cc-47da-434a-b185-98283b3baa0d
# ╠═f45d7707-1626-42e7-8211-8cdcc4cccbea
# ╟─26266ce9-1fd5-4dec-9ea7-2a79c724685d
# ╠═eb7d0765-d5e7-4ef9-916d-764c5aca9822
# ╠═5d15d24d-6317-4ad1-a53e-5af5c5bcf28a
# ╟─bd6749ab-8b75-4a45-8ac4-0650ba4a903c
# ╠═b31af4a2-716c-4c17-b283-0103630fdf38
# ╠═56016697-5ce9-4e0a-b14a-bb5d276a358c
# ╠═d261c10c-6003-45b9-91ce-7b7f9acbdc8e
# ╠═d99c6463-9ae2-479a-8c6f-ef59af17fffa
# ╠═cb8b5bdf-148e-4834-b423-79c221a43f7b
# ╟─29140461-4e06-4435-87ff-81061dc80d4c
# ╠═fc55ef51-b7f4-4640-b8e5-ad7b4f8a411c
# ╠═d64cf39a-5f91-4d3f-ad82-1aef00a0a290
# ╠═f29de74d-66dd-4d95-9d6a-c33bfb8361c9
# ╟─9b5f389f-b105-4610-bab7-f79305fedc31
# ╠═e373824d-6e0d-44c8-8c8a-192cccf10298
# ╠═e480b378-7afa-4906-9e8a-70eac7712b5e
# ╟─5beb3a0d-e57a-4aea-b7a0-59b8ce9ff5ce
# ╟─f9b4d4dc-7def-409f-b40a-f4eba1163741
# ╟─7a93e9a8-8a2d-4b11-84ef-691706c0eb0f
# ╟─86e1fffe-cf70-4d93-8131-153f0e98c94b
# ╟─4652157f-3288-45b6-877e-17ed8dc7d18e
