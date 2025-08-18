using Plots, DelimitedFiles, SpecialFunctions, QuadGK, CSV 

function sumOfExponentials(weights,poles,t)
    return weights'exp.(-poles*t)
end 

# filenames0001 = [
#             "ExpSum/data/RA_CO/alpha_0.001_m_4_tol_1.0e-5.csv", 
#             "ExpSum/data/RA_CO/alpha_0.001_m_13_tol_1.0e-10.csv", 
#             "ExpSum/data/RA_CO/alpha_0.001_m_23_tol_1.0e-15.csv", 
#             "ExpSum/data/RA_CO/alpha_0.001_m_33_tol_1.0e-20.csv", 
#             "ExpSum/data/RA_CO/alpha_0.001_m_43_tol_1.0e-25.csv", 
#             "ExpSum/data/RA_CO/alpha_0.001_m_52_tol_1.0e-30.csv", 
#             "ExpSum/data/RA_CO/alpha_0.001_m_61_tol_1.0e-35.csv", 
#             "ExpSum/data/RA_CO/alpha_0.001_m_71_tol_1.0e-40.csv", 
#             "ExpSum/data/RA_CO/alpha_0.001_m_80_tol_1.0e-45.csv",  
#             "ExpSum/data/RA_CO/alpha_0.001_m_89_tol_1.0e-50.csv",
#             "ExpSum/data/RA_CO/alpha_0.001_m_98_tol_1.0e-55.csv", 
#             ]

filenames001 = [
            "ExpSum/data/RA_CO/alpha_0.01_m_6_tol_1.0e-5.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_7_tol_1.0e-6.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_9_tol_1.0e-7.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_10_tol_1.0e-8.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_13_tol_1.0e-9.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_15_tol_1.0e-10.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_17_tol_1.0e-11.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_19_tol_1.0e-12.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_21_tol_1.0e-13.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_22_tol_1.0e-14.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_25_tol_1.0e-15.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_27_tol_1.0e-16.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_28_tol_1.0e-17.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_31_tol_1.0e-18.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_33_tol_1.0e-19.csv", 
            "ExpSum/data/RA_CO/alpha_0.01_m_34_tol_1.0e-20.csv", 
            # "ExpSum/data/RA_CO/alpha_0.01_m_25_tol_1.0e-15.csv", 
            # "ExpSum/data/RA_CO/alpha_0.01_m_34_tol_1.0e-20.csv", 
            # "ExpSum/data/RA_CO/alpha_0.01_m_44_tol_1.0e-25.csv", 
            # "ExpSum/data/RA_CO/alpha_0.01_m_53_tol_1.0e-30.csv", 
            # "ExpSum/data/RA_CO/alpha_0.01_m_64_tol_1.0e-35.csv", 
            # "ExpSum/data/RA_CO/alpha_0.01_m_72_tol_1.0e-40.csv", 
            # "ExpSum/data/RA_CO/alpha_0.01_m_81_tol_1.0e-45.csv",  
            # "ExpSum/data/RA_CO/alpha_0.01_m_91_tol_1.0e-50.csv",
            # "ExpSum/data/RA_CO/alpha_0.01_m_99_tol_1.0e-55.csv", 
            ]

filenames01 = [
            "ExpSum/data/RA_CO/alpha_0.1_m_8_tol_1.0e-5.csv", 
            "ExpSum/data/RA_CO/alpha_0.1_m_9_tol_1.0e-6.csv", 
            "ExpSum/data/RA_CO/alpha_0.1_m_11_tol_1.0e-7.csv", 
            "ExpSum/data/RA_CO/alpha_0.1_m_13_tol_1.0e-8.csv", 
            "ExpSum/data/RA_CO/alpha_0.1_m_16_tol_1.0e-9.csv", 
            "ExpSum/data/RA_CO/alpha_0.1_m_18_tol_1.0e-10.csv", 
            "ExpSum/data/RA_CO/alpha_0.1_m_19_tol_1.0e-11.csv", 
            "ExpSum/data/RA_CO/alpha_0.1_m_21_tol_1.0e-12.csv", 
            "ExpSum/data/RA_CO/alpha_0.1_m_23_tol_1.0e-13.csv", 
            "ExpSum/data/RA_CO/alpha_0.1_m_25_tol_1.0e-14.csv", 
            "ExpSum/data/RA_CO/alpha_0.1_m_28_tol_1.0e-15.csv", 
            "ExpSum/data/RA_CO/alpha_0.1_m_30_tol_1.0e-16.csv", 
            "ExpSum/data/RA_CO/alpha_0.1_m_31_tol_1.0e-17.csv", 
            "ExpSum/data/RA_CO/alpha_0.1_m_33_tol_1.0e-18.csv", 
            "ExpSum/data/RA_CO/alpha_0.1_m_35_tol_1.0e-19.csv",  
            "ExpSum/data/RA_CO/alpha_0.1_m_37_tol_1.0e-20.csv", 
            # "ExpSum/data/RA_CO/alpha_0.1_m_47_tol_1.0e-25.csv", 
            # "ExpSum/data/RA_CO/alpha_0.1_m_56_tol_1.0e-30.csv", 
            # "ExpSum/data/RA_CO/alpha_0.1_m_65_tol_1.0e-35.csv", 
            # "ExpSum/data/RA_CO/alpha_0.1_m_75_tol_1.0e-40.csv", 
            # "ExpSum/data/RA_CO/alpha_0.1_m_84_tol_1.0e-45.csv",  
            # "ExpSum/data/RA_CO/alpha_0.1_m_93_tol_1.0e-50.csv",
            # "ExpSum/data/RA_CO/alpha_0.1_m_101_tol_1.0e-55.csv", 
            ]

filenames05 = [
            "ExpSum/data/RA_CO/alpha_0.5_m_11_tol_1.0e-5.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_13_tol_1.0e-6.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_15_tol_1.0e-7.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_16_tol_1.0e-8.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_19_tol_1.0e-9.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_21_tol_1.0e-10.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_23_tol_1.0e-11.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_25_tol_1.0e-12.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_27_tol_1.0e-13.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_29_tol_1.0e-14.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_31_tol_1.0e-15.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_32_tol_1.0e-16.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_35_tol_1.0e-17.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_36_tol_1.0e-18.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_39_tol_1.0e-19.csv", 
            "ExpSum/data/RA_CO/alpha_0.5_m_40_tol_1.0e-20.csv", 
            # "ExpSum/data/RA_CO/alpha_0.5_m_50_tol_1.0e-25.csv", 
            # "ExpSum/data/RA_CO/alpha_0.5_m_59_tol_1.0e-30.csv", 
            # "ExpSum/data/RA_CO/alpha_0.5_m_68_tol_1.0e-35.csv", 
            # "ExpSum/data/RA_CO/alpha_0.5_m_78_tol_1.0e-40.csv", 
            # "ExpSum/data/RA_CO/alpha_0.5_m_87_tol_1.0e-45.csv",  
            # "ExpSum/data/RA_CO/alpha_0.5_m_95_tol_1.0e-50.csv",
            # "ExpSum/data/RA_CO/alpha_0.5_m_103_tol_1.0e-55.csv", 
            ]

filenames09 = [
            "ExpSum/data/RA_CO/alpha_0.9_m_13_tol_1.0e-5.csv", 
            "ExpSum/data/RA_CO/alpha_0.9_m_14_tol_1.0e-6.csv", 
            "ExpSum/data/RA_CO/alpha_0.9_m_16_tol_1.0e-7.csv", 
            "ExpSum/data/RA_CO/alpha_0.9_m_18_tol_1.0e-8.csv", 
            "ExpSum/data/RA_CO/alpha_0.9_m_20_tol_1.0e-9.csv", 
            "ExpSum/data/RA_CO/alpha_0.9_m_22_tol_1.0e-10.csv", 
            "ExpSum/data/RA_CO/alpha_0.9_m_24_tol_1.0e-11.csv", 
            "ExpSum/data/RA_CO/alpha_0.9_m_27_tol_1.0e-12.csv", 
            "ExpSum/data/RA_CO/alpha_0.9_m_29_tol_1.0e-13.csv", 
            "ExpSum/data/RA_CO/alpha_0.9_m_30_tol_1.0e-14.csv", 
            "ExpSum/data/RA_CO/alpha_0.9_m_32_tol_1.0e-15.csv", 
            "ExpSum/data/RA_CO/alpha_0.9_m_34_tol_1.0e-16.csv", 
            "ExpSum/data/RA_CO/alpha_0.9_m_36_tol_1.0e-17.csv", 
            "ExpSum/data/RA_CO/alpha_0.9_m_38_tol_1.0e-18.csv", 
            "ExpSum/data/RA_CO/alpha_0.9_m_39_tol_1.0e-19.csv",  
            "ExpSum/data/RA_CO/alpha_0.9_m_42_tol_1.0e-20.csv", 
            # "ExpSum/data/RA_CO/alpha_0.9_m_51_tol_1.0e-25.csv", 
            # "ExpSum/data/RA_CO/alpha_0.9_m_61_tol_1.0e-30.csv", 
            # "ExpSum/data/RA_CO/alpha_0.9_m_70_tol_1.0e-35.csv", 
            # "ExpSum/data/RA_CO/alpha_0.9_m_79_tol_1.0e-40.csv", 
            # "ExpSum/data/RA_CO/alpha_0.9_m_88_tol_1.0e-45.csv",  
            # "ExpSum/data/RA_CO/alpha_0.9_m_96_tol_1.0e-50.csv",
            # "ExpSum/data/RA_CO/alpha_0.9_m_104_tol_1.0e-55.csv", 
            ]

filenames099 = [
            "ExpSum/data/RA_CO/alpha_0.99_m_11_tol_1.0e-5.csv", 
            "ExpSum/data/RA_CO/alpha_0.99_m_13_tol_1.0e-6.csv", 
            "ExpSum/data/RA_CO/alpha_0.99_m_15_tol_1.0e-7.csv", 
            "ExpSum/data/RA_CO/alpha_0.99_m_17_tol_1.0e-8.csv", 
            "ExpSum/data/RA_CO/alpha_0.99_m_18_tol_1.0e-9.csv", 
            "ExpSum/data/RA_CO/alpha_0.99_m_21_tol_1.0e-10.csv", 
            "ExpSum/data/RA_CO/alpha_0.99_m_23_tol_1.0e-11.csv", 
            "ExpSum/data/RA_CO/alpha_0.99_m_24_tol_1.0e-12.csv", 
            "ExpSum/data/RA_CO/alpha_0.99_m_27_tol_1.0e-13.csv", 
            "ExpSum/data/RA_CO/alpha_0.99_m_29_tol_1.0e-14.csv", 
            "ExpSum/data/RA_CO/alpha_0.99_m_30_tol_1.0e-15.csv", 
            "ExpSum/data/RA_CO/alpha_0.99_m_32_tol_1.0e-16.csv", 
            "ExpSum/data/RA_CO/alpha_0.99_m_34_tol_1.0e-17.csv", 
            "ExpSum/data/RA_CO/alpha_0.99_m_36_tol_1.0e-18.csv", 
            "ExpSum/data/RA_CO/alpha_0.99_m_38_tol_1.0e-19.csv",  
            "ExpSum/data/RA_CO/alpha_0.99_m_41_tol_1.0e-20.csv", 
            # "ExpSum/data/RA_CO/alpha_0.99_m_50_tol_1.0e-25.csv", 
            # "ExpSum/data/RA_CO/alpha_0.99_m_59_tol_1.0e-30.csv", 
            # "ExpSum/data/RA_CO/alpha_0.99_m_68_tol_1.0e-35.csv", 
            # "ExpSum/data/RA_CO/alpha_0.99_m_77_tol_1.0e-40.csv", 
            # "ExpSum/data/RA_CO/alpha_0.99_m_86_tol_1.0e-45.csv",  
            # "ExpSum/data/RA_CO/alpha_0.99_m_95_tol_1.0e-50.csv",
            # "ExpSum/data/RA_CO/alpha_0.99_m_103_tol_1.0e-55.csv", 
            ]

# filenames0999 = [
#             "ExpSum/data/RA_CO/alpha_0.999_m_9_tol_1.0e-5.csv", 
#             "ExpSum/data/RA_CO/alpha_0.999_m_19_tol_1.0e-10.csv", 
#             "ExpSum/data/RA_CO/alpha_0.999_m_29_tol_1.0e-15.csv", 
#             "ExpSum/data/RA_CO/alpha_0.999_m_38_tol_1.0e-20.csv", 
#             "ExpSum/data/RA_CO/alpha_0.999_m_48_tol_1.0e-25.csv", 
#             "ExpSum/data/RA_CO/alpha_0.999_m_57_tol_1.0e-30.csv", 
#             "ExpSum/data/RA_CO/alpha_0.999_m_67_tol_1.0e-35.csv", 
#             "ExpSum/data/RA_CO/alpha_0.999_m_75_tol_1.0e-40.csv", 
#             "ExpSum/data/RA_CO/alpha_0.999_m_84_tol_1.0e-45.csv",  
#             "ExpSum/data/RA_CO/alpha_0.999_m_93_tol_1.0e-50.csv",
#             "ExpSum/data/RA_CO/alpha_0.999_m_101_tol_1.0e-55.csv", 
#             ]

filenamesDict = Dict([
    # (0.001, filenames0001), 
    (0.01, filenames001), 
    (0.1, filenames01), 
    (0.5, filenames05), 
    (0.9, filenames09), 
    (0.99, filenames099), 
    # (0.999, filenames0999)
])

colors = Dict([
    (0.001, :blue), 
    (0.01, :orange), 
    (0.1, :red), 
    (0.5, :green), 
    (0.9, :black), 
    (0.99, :grey), 
    (0.999, :yellow)
])

z = big.(range(1e-5,1,1001))
z_log = big.(10 .^(range(-5,0,1001)))

plt = plot()
for fn in filenamesDict
    
    L1errorDict = Dict()
    filenames = fn[2]

    global primary = true  
    local alpha 

    for filename in filenames

        tol = parse(Float64, split(split(filename, "_tol_")[2], ".csv")[1])
        alpha = parse(Float64, split(split(filename, "alpha_")[2], "_m_")[1])
        bigalpha = big(alpha)
                
        data = readdlm(filename, ',', BigFloat, '\n')
        coeffs = data[1,:]
        poles = data[2,:]
        
        integrand = x -> abs.(sumOfExponentials(coeffs, -poles, big(x)) - big(x).^(bigalpha-1) / gamma(bigalpha))
        @time L1error = quadgk(integrand, big.(0.00001), big.(1), rtol=1e-15)[1]
                
        scatter!([tol],[L1error], label=alpha, xscale=:log10, yscale=:log10, color = colors[alpha], primary=primary)
        
        L1errorDict[tol] = L1error 

        primary = false 
    end 
    CSV.write("ExpSum/data/RA_CO/L1error/L1error" * string(alpha) * ".csv", L1errorDict)
end 
display(plt)
