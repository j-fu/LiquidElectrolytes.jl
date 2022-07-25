var documenterSearchIndex = {"docs":
[{"location":"api/#General-remarks","page":"API","title":"General remarks","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"All physical quantities are assumed to be consistently represented through their values expressed in basic SI units (m, kg, s, A, K, mol, cd), supported by the LessUnitful.jl package built on top of Unitful.jl.","category":"page"},{"location":"api/#Electrolyte-data","page":"API","title":"Electrolyte data","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"ElectrolyteData\nAbstractElectrolyteData\ndlcap0(::ElectrolyteData)\ndebyelength(::ElectrolyteData)\nchemical_potentials!\nrrate","category":"page"},{"location":"api/#LiquidElectrolytes.ElectrolyteData","page":"API","title":"LiquidElectrolytes.ElectrolyteData","text":"mutable struct ElectrolyteData <: AbstractElectrolyteData\n\nData for electrolyte.\n\nnc::Int64\nNumber of ionic species. Default: 2\niϕ::Int64\nPotential index in species list. Default: nc + 1\nip::Int64\nPressure index in species list Default: nc + 2\nD::Vector{Float64}\nMobility coefficient Default: fill(2.0e-9 * ufac\"m^2/s\", nc)\nz::Vector{Int64}\nCharge numbers of ions Default: [(-1) ^ (i - 1) for i = 1:nc]\nM0::Float64\nMolar weight of solvent Default: 18.0153 * ufac\"g/mol\"\nM::Vector{Float64}\nMolar weight of ions Default: fill(M0, nc)\nv0::Float64\nMolar volume of solvent Default: 1 / (55.4 * ufac\"M\")\nv::Vector{Float64}\nMolar volumes of ions Default: fill(v0, nc)\nκ::Vector{Float64}\nSolvation numbers Default: fill(10.0, nc)\nc_bulk::Vector{Float64}\nBulk ion concentrations Default: fill(0.1 * ufac\"M\", nc)\nϕ_bulk::Float64\nBulk voltage Default: 0.0 * ufac\"V\"\np_bulk::Float64\nBulk pressure Default: 0.0 * ufac\"Pa\"\nΓ_bulk::Int64\nBulk boundary number Default: 1\nT::Float64\nTemperature Default: (273.15 + 25) * ufac\"K\"\nRT::Float64\nMolar gas constant scaled with temperature Default: ph\"R\" * T\nF::Float64\nFaraday constant Default: ph\"N_A*e\"\nε::Float64\nDielectric permittivity of solvent Default: 78.49\nε_0::Float64\nDielectric permittivity of vacuum Default: ph\"ε_0\"\npscale::Float64\nPressure scaling factor Default: 1.0e9\nneutralflag::Bool\nElectroneutrality switch Default: false\nscheme::Symbol\nFlux caculation scheme. This allows to choose between\n:μex (default): excess chemical potential (SEDAN) scheme, see sflux\n:act : scheme based on reciprocal activity coefficients, see aflux\n:cent : central scheme, see cflux.\nDefault: :μex\nepsreg::Float64\nRegularization parameter used in rlog  Default: 1.0e-20\n\n\n\n\n\n","category":"type"},{"location":"api/#LiquidElectrolytes.AbstractElectrolyteData","page":"API","title":"LiquidElectrolytes.AbstractElectrolyteData","text":"abstract type AbstractElectrolyteData\n\nAbstract super type for electrolytes\n\n\n\n\n\n","category":"type"},{"location":"api/#LiquidElectrolytes.dlcap0-Tuple{ElectrolyteData}","page":"API","title":"LiquidElectrolytes.dlcap0","text":"dlcap0(electrolyte)\n\nDouble layer capacitance at zero voltage for symmetric binary electrolyte.\n\nExample\n\nusing LessUnitful\nely = ElectrolyteData(c_bulk=fill(0.01ufac\"mol/dm^3\",2))\nround(dlcap0(ely),sigdigits=5) |> u\"μF/cm^2\"\n# output\n\n22.847 μF cm^-2\n\n\n\n\n\n","category":"method"},{"location":"api/#LiquidElectrolytes.debyelength-Tuple{ElectrolyteData}","page":"API","title":"LiquidElectrolytes.debyelength","text":"debyelength(electrolyte)\n\nDebye length.\n\nusing LessUnitful\nely = ElectrolyteData(c_bulk=fill(0.01ufac\"mol/dm^3\",2))\nround(debyelength(ely),sigdigits=5) |> u\"nm\"\n# output\n\n4.3018 nm\n\n\n\n\n\n","category":"method"},{"location":"api/#LiquidElectrolytes.chemical_potentials!","page":"API","title":"LiquidElectrolytes.chemical_potentials!","text":"chemical_potentials!(μ,u,electrolyte)\n\nCalculate chemical potentials from concentrations.\n\nInput:\n\nμ: memory for result (will be filled)\nu: local solution vector (concentrations, potential, pressure)\n\nReturns μ0, μ: chemical potential of solvent and chemical potentials of ions.\n\nusing LessUnitful\nely = ElectrolyteData(c_bulk=fill(0.01ufac\"mol/dm^3\",2))\nμ0,μ=chemical_potentials!([0.0,0.0],vcat(ely.c_bulk,[0,0]),ely)\nround(μ0,sigdigits=5),round.(μ,sigdigits=5)\n# output\n\n(-0.89834, [-21359.0, -21359.0])\n\n\n\n\n\n","category":"function"},{"location":"api/#LiquidElectrolytes.rrate","page":"API","title":"LiquidElectrolytes.rrate","text":"rrate(R0,β,A)\n\nReaction rate expression\n\nrrate(R0,β,A)=R0*(exp(-β*A) - exp((1-β)*A))\n\n\n\n\n\n","category":"function"},{"location":"api/#Discretization-system","page":"API","title":"Discretization system","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"PNPSystem\npnpunknowns\nelectrolytedata\nLiquidElectrolytes.bulkbcondition\nLiquidElectrolytes.solventconcentration","category":"page"},{"location":"api/#LiquidElectrolytes.PNPSystem","page":"API","title":"LiquidElectrolytes.PNPSystem","text":"PNPSystem(grid;\n         celldata=ElectrolyteData(),\n         bcondition=default_bcondition,\n         kwargs...)\n\nCreate VoronoiFVM system. Input:\n\ngrid: discretization grid\ncelldata: composite struct containing electrolyte data\nbcondition: boundary condition\nkwargs: Keyword arguments of VoronoiFVM.System\n\n\n\n\n\n","category":"function"},{"location":"api/#LiquidElectrolytes.pnpunknowns","page":"API","title":"LiquidElectrolytes.pnpunknowns","text":"pnpunknowns(sys)\n\nReturn vector of unknowns initialized with bulk data.\n\n\n\n\n\n","category":"function"},{"location":"api/#LiquidElectrolytes.electrolytedata","page":"API","title":"LiquidElectrolytes.electrolytedata","text":"electrolytedata(sys)\n\nExtract electrolyte data from system.\n\n\n\n\n\n","category":"function"},{"location":"api/#LiquidElectrolytes.bulkbcondition","page":"API","title":"LiquidElectrolytes.bulkbcondition","text":"bulkbcondition(f,u,bnode,electrolyte)\n\nBulk boundary condition for electrolyte: set potential, pressure and concentrations to bulk values.\n\n\n\n\n\n","category":"function"},{"location":"api/#LiquidElectrolytes.solventconcentration","page":"API","title":"LiquidElectrolytes.solventconcentration","text":"   solventconcentration(U::Array, electrolyte)\n\nCalculate vector of solvent concentrations from solution array.\n\n\n\n\n\n","category":"function"},{"location":"api/#Standard-calculations","page":"API","title":"Standard calculations","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"dlcapsweep\nivsweep","category":"page"},{"location":"api/#LiquidElectrolytes.dlcapsweep","page":"API","title":"LiquidElectrolytes.dlcapsweep","text":"       dlcapweep(sys;voltages=(-1:0.1:1)*ufac\"V\",\n                              δ=1.0e-4,\n                              molarity=0.1*ufac\"mol/dm^3\",\n                              solver_kwargs...)\n\nCalculate double layer capacitances for voltages given in voltages. Returns vector of voltages and vector of double layer capacitances.\n\nAssumptions:\n\nOnly one double layer in the system - close to working electrode\nSystem data has entry for voltage at working electrode ϕ_we\n1D domain\n\n\n\n\n\n","category":"function"},{"location":"api/#LiquidElectrolytes.ivsweep","page":"API","title":"LiquidElectrolytes.ivsweep","text":"   ivsweep(sys;\n              voltages=(-0.5:0.1:0.5)*ufac\"V\",\n                             ispec=1,\n                             solver_kwargs...)\n\nCalculate working electrode current corresponding to rate for species ispec for each voltage in voltages. Returns vector of voltages   and vector of currents.\n\n\n\n\n\n","category":"function"},{"location":"internal/#Electrolyte-data","page":"Internal API","title":"Electrolyte data","text":"","category":"section"},{"location":"internal/","page":"Internal API","title":"Internal API","text":"LiquidElectrolytes.charge\nLiquidElectrolytes.vrel\nLiquidElectrolytes.c0_barc\nLiquidElectrolytes.rlog\nLiquidElectrolytes.rexp\nLiquidElectrolytes.wnorm","category":"page"},{"location":"internal/#LiquidElectrolytes.charge","page":"Internal API","title":"LiquidElectrolytes.charge","text":"charge(c,electrolyte)\n\nCalculate charge from vector of concentrations\n\n\n\n\n\n","category":"function"},{"location":"internal/#LiquidElectrolytes.vrel","page":"Internal API","title":"LiquidElectrolytes.vrel","text":"vrel(ic,electrolyte)\n\nCalculate relative (wrt. solvent) molar volume of i-th species v_irel=κ_i+fracv_iv_0.\n\n\n\n\n\n","category":"function"},{"location":"internal/#LiquidElectrolytes.c0_barc","page":"Internal API","title":"LiquidElectrolytes.c0_barc","text":"c0_barc(u,electrolyte)\n\nCalculate solvent concentration c_0 and summary concentration bar c from vector of concentrations c using the incompressibility constraint (assuming κ_0=0):\n\n sum_i=0^N c_i (v_i + κ_iv_0) =1\n\nThis gives\n\n c_0v_0=1-sum_i=1^N c_i (v_i+ κ_iv_0)\n\nc_0= 1v_0 - sum_i=1^N c_i(fracv_iv_0+κ)\n\nThen we can calculate \n\n bar c= sum_i=0^N c_i\n\n\n\n\n\n","category":"function"},{"location":"internal/#LiquidElectrolytes.rlog","page":"Internal API","title":"LiquidElectrolytes.rlog","text":"rlog(u, electrolyte)\n\nCalls rlog(u;eps=electrolyte.epsreg)\n\n\n\n\n\nrlog(u; eps=1.0e-20)\n\nRegularized logarithm. Smooth linear continuation for x<eps. This means we can calculate a \"logarithm\"  of a small negative number.\n\n\n\n\n\n","category":"function"},{"location":"internal/#LiquidElectrolytes.rexp","page":"Internal API","title":"LiquidElectrolytes.rexp","text":"rexp(x;trunc=500.0)\n\nRegularized exponential. Linear continuation for x>trunc,   returns 1/rexp(-x) for x<trunc.\n\n\n\n\n\n","category":"function"},{"location":"internal/#LiquidElectrolytes.wnorm","page":"Internal API","title":"LiquidElectrolytes.wnorm","text":"wnorm(u,w,p)\n\nWeighted norm with respect to columns\n\n\n\n\n\n","category":"function"},{"location":"internal/#Finite-volume-system","page":"Internal API","title":"Finite volume system","text":"","category":"section"},{"location":"internal/","page":"Internal API","title":"Internal API","text":"LiquidElectrolytes.pnpstorage\nLiquidElectrolytes.pnpreaction\nLiquidElectrolytes.default_bcondition\nLiquidElectrolytes.pnpflux\nLiquidElectrolytes.sflux\nLiquidElectrolytes.aflux\nLiquidElectrolytes.cflux\nLiquidElectrolytes.dμex","category":"page"},{"location":"internal/#LiquidElectrolytes.pnpstorage","page":"Internal API","title":"LiquidElectrolytes.pnpstorage","text":"pnpstorage(f, u, node, electrolyte)\n\nFinite volume storage term\n\n\n\n\n\n","category":"function"},{"location":"internal/#LiquidElectrolytes.pnpreaction","page":"Internal API","title":"LiquidElectrolytes.pnpreaction","text":"pnpreaction(f, u, node, electrolyte)\n\nFinite volume reaction term\n\n\n\n\n\n","category":"function"},{"location":"internal/#LiquidElectrolytes.default_bcondition","page":"Internal API","title":"LiquidElectrolytes.default_bcondition","text":"default_bcondition(f,u,bnode,electrolyte)\n\nDefault boundary condition amounts to nothing\n\n\n\n\n\n","category":"function"},{"location":"internal/#LiquidElectrolytes.pnpflux","page":"Internal API","title":"LiquidElectrolytes.pnpflux","text":"pnpflux(f, u, edge, electrolyte)\n\nFinite volume flux. It calls either sflux, cflux or aflux.\n\n\n\n\n\n","category":"function"},{"location":"internal/#LiquidElectrolytes.sflux","page":"Internal API","title":"LiquidElectrolytes.sflux","text":"sflux(ic,dϕ,ck,cl,βk,βl,bar_ck,bar_cl,electrolyte)\n\nSedan flux,  see Gaudeul/Fuhrmann 2022\n\nAppearantly first described by Yu, Zhiping  and Dutton, Robert, SEDAN III, www-tcad.stanford.edu/tcad/programs/sedan3.html\n\nsee also the 198? Fortran code available via  https://web.archive.org/web/20210518233152/http://www-tcad.stanford.edu/tcad/programs/oldftpable.html\n\nVerification calculation is in the paper.\n\n\n\n\n\n","category":"function"},{"location":"internal/#LiquidElectrolytes.aflux","page":"Internal API","title":"LiquidElectrolytes.aflux","text":"aflux(ic,dϕ,ck,cl,βk,βl,bar_ck,bar_cl,electrolyte)\n\nFlux expression based on reciprocal activity coefficents, see Fuhrmann, CPC 2015\n\n\n\n\n\n","category":"function"},{"location":"internal/#LiquidElectrolytes.cflux","page":"Internal API","title":"LiquidElectrolytes.cflux","text":"aflux(ic,dϕ,ck,cl,βk,βl,bar_ck,bar_cl,electrolyte)\n\nFlux expression based on centrals differences, see Gaudeul/Fuhrmann 2022, Cances\n\n\n\n\n\n","category":"function"},{"location":"internal/#LiquidElectrolytes.dμex","page":"Internal API","title":"LiquidElectrolytes.dμex","text":" dμex(βk, βl, electrolyte)\n\nCalculate differences of excess chemical potentials from reciprocal activity coefficient\n\n\n\n\n\n","category":"function"},{"location":"internal/#Electrochemical-calculations","page":"Internal API","title":"Electrochemical calculations","text":"","category":"section"},{"location":"internal/","page":"Internal API","title":"Internal API","text":"LiquidElectrolytes.splitz","category":"page"},{"location":"internal/#LiquidElectrolytes.splitz","page":"Internal API","title":"LiquidElectrolytes.splitz","text":"splitz(range::AbstractRange)\n\nIf range contains zero, split it into two parts, one with values <=0 and one with values >=0. Otherwise, return the range or its reverse, such that first value always is the one with the smallest absolute value.\n\n\n\n\n\nsplitz(range::Vector)\n\nVersion of splitz(range::AbstractRange) for vectors.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"using Markdown\nMarkdown.parse(\"\"\"\n$(read(\"../../README.md\",String))\n\"\"\")","category":"page"}]
}