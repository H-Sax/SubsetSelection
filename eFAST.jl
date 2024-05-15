using Pkg
Pkg.add("Dates")
using Dates
println("Running on ", Threads.nthreads(), " threads.")
println(now())

Pkg.add("DifferentialEquations")
Pkg.add("Statistics")
Pkg.add("Distributions")
Pkg.add("QuasiMonteCarlo")
Pkg.add("JLD")
@time "Setup modules" using DifferentialEquations, Statistics, Distributions, QuasiMonteCarlo, JLD, GlobalSensitivity

function Valve(R, deltaP, open)
    dq = 0.0
    if (-open) < 0.0 
        dq =  deltaP/R
    else
        dq = 0.0
    end
    return dq

end

function ShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)
    τₑₛ = τₑₛ*τ
    τₑₚ = τₑₚ*τ
    #τ = 4/3(τₑₛ+τₑₚ)
    tᵢ = rem(t + (1 - Eshift) * τ, τ)

    Eₚ = (tᵢ <= τₑₛ) * (1 - cos(tᵢ / τₑₛ * pi)) / 2 +
         (tᵢ > τₑₛ) * (tᵢ <= τₑₚ) * (1 + cos((tᵢ - τₑₛ) / (τₑₚ - τₑₛ) * pi)) / 2 +
         (tᵢ <= τₑₚ) * 0

    E = Eₘᵢₙ + (Eₘₐₓ - Eₘᵢₙ) * Eₚ

    return E
end

function DShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)

    τₑₛ = τₑₛ*τ
    τₑₚ = τₑₚ*τ
    #τ = 4/3(τₑₛ+τₑₚ)
    tᵢ = rem(t + (1 - Eshift) * τ, τ)

    DEₚ = (tᵢ <= τₑₛ) * pi / τₑₛ * sin(tᵢ / τₑₛ * pi) / 2 +
          (tᵢ > τₑₛ) * (tᵢ <= τₑₚ) * pi / (τₑₚ - τₑₛ) * sin((τₑₛ - tᵢ) / (τₑₚ - τₑₛ) * pi) / 2
    (tᵢ <= τₑₚ) * 0
    DE = (Eₘₐₓ - Eₘᵢₙ) * DEₚ

    return DE
end

#Shi timing parameters
Eshift = 0.0
Eₘᵢₙ = 0.03

τₑₛ = 0.3
τₑₚ = 0.45 
Eₘₐₓ = 1.5
Rmv = 0.006
τ = 1.0
 
function NIK!(du, u, p, t)
    pLV, psa, psv, Vlv, Qav, Qmv, Qs = u 
    τₑₛ, τₑₚ, Rmv, Zao, Rs, Csa, Csv, Eₘₐₓ, Eₘᵢₙ = p
    # pressures (more readable names)
# the differential equations
    du[1] = (Qmv - Qav) * ShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift) + pLV / ShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift) * DShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)
    # 1 Left Ventricle
    du[2] = (Qav - Qs ) / Csa #Systemic arteries     
    du[3] = (Qs - Qmv) / Csv # Venous
    du[4] = Qmv - Qav # volume
    du[5]    = Valve(Zao, (du[1] - du[2]), u[1] - u[2])  # AV 
    du[6]   = Valve(Rmv, (du[3] - du[1]), u[3] - u[1])  # MV
    du[7]     = (du[2] - du[3]) / Rs # Systemic flow
end
##

u0 = [6.0, 6.0, 6.0, 200.0, 0.0, 0.0, 0.0]

p = [0.3, 0.45, 0.006, 0.033, 1.11, 1.13, 11.0, 1.5, 0.03]

tspan = (0, 10)

prob = ODEProblem(NIK!, u0, tspan, p)

x = LinRange(7,8,250)
@time sol = solve(prob, Vern7(),  reltol = 1e-5, abstol = 1e-5, saveat = x)

circ_time_post = function (p)
    out_global
end

###### Global Sensitiivty analysis on continous outputs ####
# LV pressure, SA pressure, LV volume 

circ_time_idxs_chunked = function (p)

    global parameterset = p

    Np::Int64 = size(p,2)

    println("Running ", Np, " cases.")

    chunksize::Int64 = Np/chunks

    println("Using ", chunks, " chunk(s) of size ", chunksize, ".")


    out = zeros(4,Np)

    for k in 1:chunks
        offset = Int((k - 1) * Np/chunks)
        startindx = Int(offset + 1)
        endindx = Int(k * Np/chunks)

        println("Starting chunk ", k, ", from ", startindx, " to ", endindx, " with offset ", offset, ".")

        pchunk = p[:,startindx:endindx]

        prob_func(prob,i,repeat) = remake(prob; u0 = [8.0, 8.0, 8.0, 265.0, 0.0, 0.0, 0.0], p=pchunk[:,i])

        ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)

        @time "    ODE Ensemble solved in " sol = solve(ensemble_prob, Vern7(),  reltol = 1e-5, abstol = 1e-5, EnsembleThreads();saveat = x,trajectories=chunksize)
        @time "    Results processed in" Threads.@threads for i in 1:chunksize
          # println(i+offset)
          out[1,i+offset] = (maximum(sol[i][4,:]) - minimum(sol[i][4,:])) #SV LV
          out[2,i+offset] = (maximum(sol[i][1,:]) - minimum(sol[i][1,:])) #PP LV
          out[3,i+offset] = (maximum(sol[i][2,:]) - minimum(sol[i][2,:])) #PP SA.P
          out[4,i+offset] = mean(sol[i][7,:]) #Mean systemic flow
        end
    end
    global out_global = out
    out
end

println("Running eFAST")
chunks::Int64 = 50000
lb = [0.255, 0.3825, 0.0051, 0.02805, 0.9435, 0.9605, 9.350, 1.275, 0.0255]
ub = [0.345, 0.5175, 0.0069, 0.03795, 1.2765, 1.2995, 12.65, 1.725, 0.0345]
@time res = gsa(circ_time_idxs_chunked,eFAST(), [[lb[i],ub[i]] for i in 1:9], samples=100000, batch=true)

#save("Nik_eFAST.jld", "data", res)


# Analysis of eFAST indices 
res = load("/home/harry/Desktop/PhD/Year 3/Li rewrite/eFAST/Nik_eFAST.jld")["data"]
# First Order
S1 = res.S1
# Total order 
ST = res.ST

# Plot the influence
using CairoMakie
CairoMakie.activate!(type = "svg")
begin
    f = Figure(size = (1100, 600), backgroundcolor = RGBf(0.98, 0.98, 0.98))
    #f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98));

    ax1 = Axis(f[1,1], xticklabelrotation = π / 3, xticklabelalign = (:right, :center), xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), yticks = (1:4, ["SV - LV", "PP - LV", "PP - SA", "Mean - Qs"]), title = "First Order Sensitivity Matrix", ylabel = "Measurements", xlabel = "Parameters")
    Label(f[1, 1, TopLeft()], "A",fontsize = 18,font = :bold,halign = :right)
    hm1 = CairoMakie.heatmap!(ax1,transpose(S1), colormap=:plasma)
    for i in 1:9, j in 1:4
        txtcolor = transpose(S1)[i, j] < -1000.0 ? :white : :black
        text!(ax1, "$(round(transpose(S1)[i,j], digits = 4))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 12)
    end


    ax2 = Axis(f[1,3], xticklabelrotation = π / 3, xticklabelalign = (:right, :center), xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), title = "Total Order Sensitivity Matrix",  xlabel = "Parameters")
    hideydecorations!(ax2)
    Label(f[1, 3, TopLeft()], "B",fontsize = 18,font = :bold,halign = :right)
    hm2 = CairoMakie.heatmap!(ax2,transpose(ST), colormap=:plasma)
    for i in 1:9, j in 1:4
        txtcolor = transpose(ST)[i, j] < -1000.0 ? :white : :black
        text!(ax2, "$(round(transpose(ST)[i,j], digits = 4))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 12)
    end
    Colorbar(f[1,2],hm2)

    gl = GridLayout(f[2, 1:3], width = Relative(0.5))
    A = transpose(ST - S1)
    ax3 = Axis(gl[1, 1], xticklabelrotation = π / 3, xticklabelalign = (:right, :center), xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), yticks = (1:4, ["SV - LV", "PP - LV", "PP - SA", "Mean - Qs"]), title = "Difference Sensitivity Matrix", ylabel = "Measurements",  xlabel = "Parameters")
    Label(gl[1, 1, TopLeft()], "C", fontsize = 18,font = :bold,halign = :right)
    hm3 = CairoMakie.heatmap!(ax3,A, colormap=:plasma)
    for i in 1:9, j in 1:4
        txtcolor = A[i, j] < -1000.0 ? :white : :black
        text!(ax3, "$(round(A[i,j], digits = 4))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 12)
    end
    Colorbar(gl[1,2],hm3)
    linkaxes!(ax1,ax2,ax3)
    f
end 


# Calculate the Orthogonality Matrices 
Orth_heat = Matrix{Float64}(undef,9,9)
for j in 1:9
    for i in 1:9
        if i==j
            Orth_heat[i,j] = 0
        else 
        Orth_heat[i,j] = sin(acos(((transpose(ST[:,i])*ST[:,j]))/(norm(ST[:,i])*norm(ST[:,j]))-1e-15))  #Slight numerical rounding error without the additional add on 
        end 
    end 
end 

#save("/home/harry/Desktop/PhD/Year 3/Li rewrite/Orthogonality/eFAST_ST.jld","data",Orth_heat)

# Calculate the average orth for each parameter 
rank = zeros(9)
for j in 1:9
    rank[j] = sum(Orth_heat[:,j])/8
end 
p=sortperm(rank,rev=true)
using CairoMakie
f = Figure();
ax = Axis(f[1,1],title="eFAST ST Sensitivity Matrix-Parameter Orthogonality Rank", xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"][p]), xlabel = "Parameters", ylabel = "Orthogonality")
CairoMakie.scatter!( rank[p])
f


# Calculate overall parameter influence Ej 
F = transpose(ST)*ST
e_deomp=eigen(F)

λ = abs.(e_deomp.values)
Q = abs.(e_deomp.vectors)

e_value_sum = sum(λ)

e = Vector{Float64}(undef,9)
for i in 1:9
    for j in 1:9
    e[i] = sum(λ[j]*Q[i,j])/e_value_sum
    end 
end 
p=sortperm(e,rev=true)
using CairoMakie
f = Figure();
ax = Axis(f[1,1],title="eFAST ST Sensitivity Matrix-Parameter Importance PCA Method", xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"][p]), xlabel = "Parameters", ylabel = "Importance")
CairoMakie.scatter!( e[p])
f

S = Array(transpose(ST))
using MAT
file = matopen("e_ST.mat", "w")
write(file, "e_ST", S)
close(file)