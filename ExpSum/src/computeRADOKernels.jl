include("hpAAADO.jl")
import .hpAAADO

using RationalFunctionApproximation, GenericLinearAlgebra, LinearAlgebra, Plots, SpecialFunctions, DelimitedFiles
using ProgressBars
using QuadGK

####################################################################################
# Gamma(6-a) / 120 from 0 to 2 

function distGamma6(alpha)
    return gamma(big.(6 - alpha)) / 120 
end 

function distGamma5(alpha)
    return gamma(big.(5 - alpha)) / 120 
end 

# Gamma(4-a) from 0 to 2

function distGamma4(alpha)
    return gamma(big.(4 - alpha)) 
end 

function distGamma3(alpha)
    return gamma(big.(3 - alpha)) 
end 

# Gamma(4-a) sinh(a) from 0 to 2 

function distGamma4sinh(alpha)
    return gamma(big.(4-alpha)) * sinh(alpha)
end 

function distGamma3sinh(alpha)
    return gamma(big.(3-alpha)) * sinh(1 + alpha)
end 

# c^a from 0 to 1 for different c 

function dist01TTPalpha(alpha)
    return big.(0.1)^alpha
end 

function dist05TTPalpha(alpha)
    return big.(0.5)^alpha
end 

function dist2TTPalpha(alpha)
    return big.(2)^alpha
end 

# 6 a (1-a) from 0 to 1 

function distQuadratic(alpha)
    return alpha * (1-alpha)
end 

# exp(-a) Gamma(6-a) from 0 to 2 

function distExpGamma6(alpha)
    return exp(-alpha) * gamma(big.(6-alpha))
end 

function distExpGamma5(alpha)
    return exp(-alpha-1) * gamma(big.(5-alpha))
end 

# Gamma(3-a) from 0.2 to 1.5 

function distGamma3from02to1(alpha)
    if alpha < 0.2
        return 0
    else 
        return gamma(big.(3-alpha))
    end 
end 

function distGamma2from0to05(alpha)
    if alpha > 0.5
        return 0
    else 
        return gamma(big.(2-alpha))
    end 
end 

# Bumps with support [x-0.1, x+0.1]

function distr01Bump(alpha, center)
    if center - 0.1 < alpha && alpha < center + 0.1 
        return exp(-1 / (1 - ((alpha-center)/0.1)^2)) * 22.5228
    else 
        return 0 
    end 
end 

# Bumps with support [x-0.5, x+0.5]
function distr05Bump(alpha, center)
    if center - 0.5 < alpha && alpha < center + 0.5 
        return exp(-1 / (1 - ((alpha-center)/0.5)^2)) * 4.504565
    else 
        return 0 
    end 
end 

function uniform(alpha)
    return 1 
end 

####################################################################################

QForder = 15
        
nz = 200
zexp = 8 

# for visualization only 
nt = 1000 
texp = 1 

for distribution in [distGamma6, distGamma5]
    for AAAtol in [1e-10,1e-20,1e-30] 

        println(center)

        QFtol = AAAtol * 1e-10 
        dist = distribution

        z, target = hpAAADO.computeTarget(dist, zexp, nz, QFtol, QForder, 0., 1.) 
        t, kernel = hpAAADO.computeKernel(dist, texp, nt, QFtol, QForder, 0., 1.)
        
        println("Compute rational approximation (AAA)")
        @time r = aaa(z, target, tol = AAAtol, max_degree=20)

        print("Compute partial fractions decomposition, ")
        mypoles, mycoeffs = hpAAADO.computePolesCoeffs(r) 
        println("m: ", size(mypoles))

        Datapath = "ExpSum/data/RA_DO/"
        Imagepath = "ExpSum/Images/RA_DO/"

        filename = String(Symbol(dist)) 
        filename *= "_center_" * string(center)
        filename *= "_m_" * string(size(mypoles)[1])
        filename *= "_AAAtol_" * string(AAAtol) 

        hpAAADO.displayTarget(r,z,target, Imagepath * "Target_" * filename * ".pdf")
        hpAAADO.displayKernel(mycoeffs, mypoles, t, kernel, Imagepath * "Kernel_" * filename * ".pdf")

        println("Compute L1 error")
        @time L1error = hpAAADO.computeL1error(dist, mypoles, mycoeffs, AAAtol, QForder, 0., 1.)

        writedlm(Datapath * filename * ".csv", [real(mycoeffs), real(mypoles)], ',')
        writedlm(Datapath * "L1errors/" * filename * ".csv", [L1error], ',')
        
        println()
    end 
end 