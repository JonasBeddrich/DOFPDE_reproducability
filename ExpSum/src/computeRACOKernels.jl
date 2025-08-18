using RationalFunctionApproximation, GenericLinearAlgebra, LinearAlgebra, Plots, SpecialFunctions, DelimitedFiles

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

function removeZeroFromZeros(zeros)
    return 
end 

z = big.(unique([big.(0); 
            big.(range(1e0,1e1,1000)); 
            big.(range(1e1,1e2,1000)); 
            big.(range(1e2,1e3,1000)); 
            big.(range(1e3,1e4,1000)); 
            big.(range(1e4,1e5,1000)); 
            big.(range(1e5,1e6,1000)); 
            big.(range(1e6,1e7,100)); 
            big.(range(1e7,1e8,10))]))

z = big.(range(0,1000000,10001))

limit = 7
n = 200
z = unique([big(0.); big.(10^limit .- 10 .^range(1,limit,n)); big.(10 .^range(1,limit,n))])
z2 = big.(10 .^(range(-10,0,10001)))

plot()
for alpha in [0.001,0.01,0.1,0.5,0.9,0.990,999]
    println(alpha)
    for tol in [1e-17] 
        alpha=big.(alpha)
        println(tol)
        func = x -> x .^(1-alpha)
        f = func.(z) 
        
        @time local r = aaa(z,f,tol=tol)
        
        local mypoles = poles(r)
        local myzeros = sort(real(computeZeros(r)))[1:end-1]
        
        local myw = weights(r)
        local myf = values(r)
        global mynodes = sort(nodes(r))

        # plt = scatter(z[2:end],z[2:end])
        # scatter!(mynodes[2:end], mynodes[2:end], xscale=:log10, yscale=:log10)
        # display(plt)
        # plot()
        # println(sort(mynodes))
        # println(mynodes)

        println("m", size(mynodes))

        local scaling = computeScale(myw, myf)
        local mycoeffs = computePartialFractionDecomposition(real(mypoles), real(myzeros), scaling)

        # plot(z2, x -> x.^-alpha / gamma(1-alpha), yscale=:log10)
        # plot!(z2, x -> sumOfExponentials(mycoeffs, -real(mypoles), x))

        plot!(z2, x -> abs(x.^(alpha-1) / gamma(alpha) - sumOfExponentials(mycoeffs, -real(mypoles), x)), 
        xscale=:log10, yscale=:log10, 
        # label="alpha: " * string(alpha) * " tol: " * string(tol), 
        ylim=[1e-25,1e5])

        writedlm("ExpSum/data/RA_CO/alpha_" * string(alpha) * "_m_" * string(size(mynodes)[1]) * "_tol_" * string(tol) * ".csv", 
        [real(mycoeffs), real(mypoles)], ',')
    end
end 
plot!()
