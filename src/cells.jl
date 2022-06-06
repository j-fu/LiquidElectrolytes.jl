function splitz(range::AbstractRange)
    if range[1]>=0
        return [range]
    elseif range[end]<=0
        return [reverse(range)]
    else
        [0:-step(range):range[1],0:step(range):range[end]]
    end
end

function splitz(range::Vector)
    if range[1]>=0
        return [vcat([0.0],range)]
    elseif range[end]<=0
        return [vcat([0.0],reverse(range))]
    else
        for i=1:length(range)
            if range[i]≈0.0
                return [ reverse(range[1:i]), range[i:end]]
            elseif i>1 && range[i-1]<0.0 && range[i]>0.0
                return [ vcat([0.0],reverse(range[1:i-1])), vcat([0.0],range[i:end])]
            end
        end
    end

end

function doublelayercap(sys;voltages=-1:0.1:1,δ=1.0e-4,molarity=0.1,solver_kwargs...)
    ranges=splitz(voltages)

    vplus = zeros(0)
    cdlplus = zeros(0)
    vminus = zeros(0)
    cdlminus = zeros(0)
    data=electrolytedata(sys)
    data.ϕ_we=0
    
    data.c_bulk.=molarity*mol/dm^3
    
    iϕ=data.iϕ
    inival0 = solve(sys,inival=pnpunknowns(sys))
    inival=copy(inival0)
    sol=copy(inival0)
    allprogress=sum(length,ranges)
    ϕprogress = 0
    
    function show_error(u,δ)
        @show cond(inv(Diagonal(sys.matrix))*Matrix(sys.matrix))
        @show u[1,1:5]
        @show u[2,1:5]
        @show u[3,1:5]
        @show u[4,1:5]
        @show chemical_potentials!(zeros(2),u[:,1],data)
        c0,barc= c0_barc(u[:,1],data)
        @show c0/barc, u[1,1]/barc,u[2,1]/barc
        @error "bailing out at δ=$(δ) ϕ_we=$(data.ϕ_we)V, molarity=$(molarity)"
    end

    
    control=VoronoiFVM.SolverControl(max_round=3, tol_round=1.0e-9;solver_kwargs...)
    @withprogress for range in ranges
        sol .= inival0
        success=true
        for ϕ in range 
            try
                data.ϕ_we= ϕ
                inival.=sol
                s=data.scheme
                solve!(sol, inival, sys; control)
            catch e
                println(e)
                show_error(inival,0)
                success=false
            end
#            @show extrema(sol[1,:])
            if !success
                break
            end
            Q = integrate(sys, sys.physics.reaction, sol)
            
            try
                data.ϕ_we=ϕ+δ
                inival .= sol
                solve!(sol, inival, sys; control)
            catch e
                println(e)
                show_error(inival,δ)
                success=false
            end
            if !success
                break
            end
            Qδ = integrate(sys, sys.physics.reaction, sol)
            cdl = (Qδ[iϕ] - Q[iϕ]) / δ

            if range[end]>range[1]
                push!(vplus, ϕ)
                push!(cdlplus, cdl)
            else
                push!(vminus, ϕ)
                push!(cdlminus, cdl)
            end
            ϕprogress +=1
            @logprogress ϕprogress/allprogress
        end
    end
    vcat(reverse(vminus),vplus),vcat(reverse(cdlminus),cdlplus) 
end

    




function voltagesweep(sys;voltages=-0.5:0.1:0.5,ispec=1,solver_kwargs...)
    ranges=splitz(voltages)

    factory=VoronoiFVM.TestFunctionFactory(sys)
    data=sys.physics.data
    
    tf=testfunction(factory,[data.Γ_bulk],[data.Γ_we] )


    vplus = zeros(0)
    iplus = zeros(0)
    splus = []
    vminus = zeros(0)
    iminus = zeros(0)
    sminus = []
    weights=ones(data.nc+2)
    weights[data.ip]=0
    mynorm=u->wnorm(u,weights,Inf)
    myrnorm=u->wnorm(u,weights,1)
    data=electrolytedata(sys)
    data.ϕ_we=0
    control=SolverControl(;solver_kwargs...)
    iϕ=data.iϕ
    inival0 = solve(sys;inival=pnpunknowns(sys), solver_kwargs...)
    inival=copy(inival0)
    sol=copy(inival0)
    ϕprogress=0
    function show_error(u)
        @show u[1,1:5]
        @show u[2,1:5]
        @show u[3,1:5]
        @show u[4,1:5]
        @show chemical_potentials!(zeros(size(u,1)),u[:,1],data)
        c0,barc= c0_barc(u[:,1],data)
        @show c0/barc, u[1,1]/barc,u[2,1]/barc
        @error "bailing out at ϕ_we=$(data.ϕ_we)V"
    end

    allprogress=sum(length,ranges)
    ϕprogress = 0

    
    @withprogress for range in ranges
        inival .= inival0
        for ϕ in range 
            data.ϕ_we=ϕ
            try
                solve!(sol, inival, sys;mynorm,myrnorm,control)
            catch e
                println(e)
                show_error(sol)
                break
            end
            inival .= sol
            I=-integrate(sys,sys.physics.breaction,sol; boundary=true)[:,data.Γ_we]
            if range[end]>range[1]
                push!(vplus, ϕ)
                push!(iplus, I[ispec]*F/(mA/cm^2))
                push!(splus, copy(sol))
            else
                push!(vminus, ϕ)
                push!(iminus, I[ispec]*F/(mA/cm^2))
                push!(sminus, copy(sol))
            end
            ϕprogress +=1
            @logprogress ϕprogress/allprogress
        end
    end
    popfirst!(iminus)
    popfirst!(vminus)
    popfirst!(sminus)
    vcat(reverse(vminus),vplus),vcat(reverse(iminus),iplus), vcat(reverse(sminus),splus) 
end

    
