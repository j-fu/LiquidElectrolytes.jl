module Example220_SimpleCV
using Printf
using LessUnitful
using LiquidElectrolytes
using LiquidElectrolytes.RotatingElectrodes
using VoronoiFVM

@phconstants AvogadroConstant MolarGasConstant ElementaryCharge
const FaradayConstant = AvogadroConstant * ElementaryCharge
@unitfactors dm s cm mol K V mV

function main(;
        frequencies = [1, 4, 9, 16.0],
        Plotter = nothing,
        plotres = true,
        verbose = false,
        kreact = 1.0,
        scanrate = 10,
        num_periods = 1,
        file = "simplecv2.pdf"
    )
    Temperature = 273.15 * K + 25 * K
    FRT = FaradayConstant / (MolarGasConstant * Temperature)

    cin = 1.0 * mol / dm^3        # inlet  concentration
    D_A = 5.0e-6 * cm^2 / s       # diffusion coefficient of species A
    D_B = 5.0e-6 * cm^2 / s       # diffusion coefficient of species B

    specA = 1
    specB = 2

    eplus::Float64 = 0.0
    eminus::Float64 = 0.0

    kp::Float64 = kreact
    km::Float64 = kreact
    E0 = 0

    alpha = 0.5

    function set_voltage(deltaPhi)
        eplus = exp(-    alpha * FRT * (deltaPhi - E0))
        return eminus = exp((1.0 - alpha) * FRT * (deltaPhi - E0))
    end


    function reaction(f, u, node, data)
        return if node.region == b_disk
            r = kp * eplus * u[specA] - km * eminus * u[specB]
            f[specA] = r
            f[specB] = -r
        end
    end

    cvcontext = CVContext(
        reaction, "E7R8",
        set_voltage = set_voltage,
        num_species = 2,
        I_disk = (I) -> I[specB],
        I_ring = (I) -> -I[specB],
    )
    viscosity!(cvcontext, 0.02942 * cm^2 / s)
    diffusion!(cvcontext, specA, D_A)
    diffusion!(cvcontext, specB, D_B)
    boundary_dirichlet!(cvcontext, specA, b_in, cin)
    boundary_dirichlet!(cvcontext, specB, b_ring, 0)

    control = solvercontrol()

    control.Δt = 2.5
    control.Δt_min = 1.0e-5
    control.Δt_max = 2.5
    control.Δt_grow = 1.2
    control.Δu_opt = 1
    control.verbose = false
    control.max_lureuse = 0


    cvresults = CVResults()
    for frequency in frequencies
        add_cvresult!(
            cvresults, run_cv(
                cvcontext,
                frequency = frequency,
                phi_min = 0.5 * V,
                phi_max = -0.5 * V,
                scanrate = scanrate * mV / s,
                num_periods = num_periods,
                verbose = verbose,
                control = control
            )
        )
    end

    return if plotres
        plot_cvresults(cvresults; Plotter, file = file)
    end
end


end
