module hpAAADO

    using RationalFunctionApproximation, GenericLinearAlgebra, LinearAlgebra, Plots, SpecialFunctions, DelimitedFiles
    using ProgressBars
    using QuadGK

    function computeZeros(r::Barycentric)
        wf = r.w_times_f
        m = length(wf)
        ty = eltype(nodes(r))
        B = diagm( [ty(0); ones(m)] )
        # Thanks to Daan Huybrechs:
        E = [0 transpose(wf); ones(m) diagm(nodes(r))]
        # super kludgy; thanks to Daan Huybrechs
        μ = 11 - 17im  # shift since 0 may well be a root
        E -= μ*B
        EB = inv(E)*B
        pp = eigvals(EB)
        large_enough = abs.(pp) .> 1e-10
        return μ .+ 1 ./ pp[large_enough]
    end

    function computePartialFractionDecomposition(poles, zeros, scale)
        coeffs = scale * ones(size(poles))
        for (idp, pole) in enumerate(poles)
            local numerator = 1
            for (idz, zero) in enumerate(zeros)
                numerator *= pole-zero
            end 
            local denominator = 1
            for (idp2, pole2) in enumerate(poles)
                if idp != idp2
                    denominator *=  pole - pole2
                end 
            end 
            coeffs[idp] *= numerator / denominator 
        end   
        return coeffs
    end 

    function sumOfExponentials(weights,poles,t)
        return weights'exp.(-poles*t)
    end 

    function computeScale(w,f) 
        return w'f / sum(w)
    end 

    function computePolesCoeffs(r)
        mypoles = real(poles(r))
        myzeros = sort(real(computeZeros(r)))
        cleanedzeros = myzeros[1:end-1]

        myw = weights(r)
        myf = values(r)
        scaling = computeScale(myw, myf)

        mycoeffs = computePartialFractionDecomposition(mypoles, real(cleanedzeros), scaling)
        return mypoles, mycoeffs
    end 

    function computeTarget(distribution, zexp, nz, QFtol, QForder, lower, upper)
        z = big.(sort(unique([0.; 10^zexp .- 10 .^range(0,zexp,nz); 10 .^range(0,zexp,nz)])))
        target = big.(zeros(size(z))) 

        println("Compute RA target")
        for idx in ProgressBar(eachindex(z[2:end]))
            local integrand = x -> distribution(x) .* z[idx+1] .^ (x)
            target[idx+1] = quadgk(integrand, big.(lower), big.(upper), atol=QFtol, order=QForder)[1]
        end 
        return z, target
    end 

    function computeKernel(distribution, texp, nt, QFtol, QForder, lower, upper)
        t = big.(10 .^(range(-10,texp,nt)))
        kernel = big.(zeros(size(t)))

        println("Compute DO kernel")
        for idx in ProgressBar(eachindex(t))
            local integrand = x -> distribution(x) ./ gamma(1-x) .* t[idx] .^ -x
            kernel[idx] = quadgk(integrand, big.(lower), big.(upper), atol=QFtol, order=QForder)[1]
        end 
        return t, kernel 
    end 

    function computeL1error(distribution, mypoles, mycoeffs, AAAtol, QForder, lower, upper)
        SOE = x -> hpAAADO.sumOfExponentials(mycoeffs, -mypoles, x)

        function time_integrand(t)
            alpha_integrand = a -> distribution(a) ./ gamma(1-a) .* t .^ -a
            return abs.(SOE.(t) - quadgk(alpha_integrand, big.(lower), big.(upper), atol = AAAtol, order=QForder)[1])
        end 
        integral, error, count = @time quadgk_count(time_integrand, big.(0.00001), big.(1), rtol=1e-3, order = QForder)
        println("integral value: ", integral)
        println("integral uncertainty: ", error)
        println("integral evaluations: ", count)
        return integral
    end 

    function displayTarget(r, z, target, filename)
        RA = x -> r(x)
        z = z[2:end]
        target = target[2:end] 
        p1 = scatter(z, target, xscale=:log10, yscale=:log10, label = "Target")
        plot!(z, RA.(z), label="RA", lw=3)
        target_error = abs.(target-RA.(z))
        p2 = plot(z[isless.(0., target_error)], target_error[isless.(0., target_error)], 
                  xscale=:log10, yscale=:log10, label = "RA-Target", lw=3)
        p = plot(p1,p2)
        display(p)
        savefig(p, filename)
    end
    
    function displayKernel(mycoeffs, mypoles, t, kernel, filename)
        SOE = x -> hpAAADO.sumOfExponentials(mycoeffs, -mypoles, x)
        p1 = scatter(t, kernel, xscale=:log10, yscale=:log10, label = "Kernel")
        plot!(t, SOE.(t), lw=3, label = "Sum of Exp")
        p2 = plot(t, abs.(kernel - SOE.(t)), xscale=:log10, yscale=:log10, label = "|Kernel-SOE|", lw=3)
        p = plot(p1,p2)
        display(p)
        savefig(p, filename)
    end 
end 