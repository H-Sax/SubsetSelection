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

println("Running Morris")
chunks::Int64 = 50000
lb = [0.255, 0.3825, 0.0051, 0.02805, 0.9435, 0.9605, 9.350, 1.275, 0.0255]
ub = [0.345, 0.5175, 0.0069, 0.03795, 1.2765, 1.2995, 12.65, 1.725, 0.0345]
@time res = gsa(circ_time_idxs_chunked,Morris(relative_scale = true,num_trajectory=100000),[[lb[i],ub[i]] for i in 1:9],batch = true)
#save("Nik_Morris.jld", "data", res)
using JLD, GlobalSensitivity
# Analysis of Morris indices 
res = load("/home/harry/Desktop/PhD/Year 3/Li rewrite/Morris/Nik_Morris.jld")["data"]
# Absolute mean
mu = abs.(res[:,:,1] / maximum(res[:,:,1]))
# Variances
var = abs.(res[:,:,2] / maximum(res[:,:,2]))


# Plot Mean vs Variances
using CairoMakie
CairoMakie.activate!(type = "svg")
begin
    fig = Figure(size = (620, 380))
    ax1 = Axis(fig[1,1],title = "SV - LV", xscale = log10, yscale = log10, ylabel = "Parameter Variance")
    ax2 = Axis(fig[1,2], title = "PP - LV", xscale = log10, yscale = log10)
    ax3 = Axis(fig[2,1], title = "PP - SA", xscale = log10, yscale = log10, xlabel = "Parameter Mean Influence", ylabel = "Parameter Variance")
    ax4 = Axis(fig[2,2], title = "Mean - Qs", xscale = log10, yscale = log10, xlabel = "Parameter Mean Influence")
    linkaxes!(ax1,ax2,ax3,ax4)
    hidexdecorations!(ax1, grid = false)
    hidedecorations!(ax2, grid = false)
    hideydecorations!(ax4, grid = false)
    # SV - LV Plot
    x1 = Point2f.(mu[1,1], var[1,1])
    x2 = Point2f.(mu[1,2], var[1,2])
    x3 = Point2f.(mu[1,3], var[1,3])
    x4 = Point2f.(mu[1,4], var[1,4])
    x5 = Point2f.(mu[1,5], var[1,5])
    x6 = Point2f.(mu[1,6], var[1,6])
    x7 = Point2f.(mu[1,7], var[1,7])
    x8 = Point2f.(mu[1,8], var[1,8])
    x9 = Point2f.(mu[1,9], var[1,9])
    scatter!(ax1,x1, label = L"τ_{es}", color = :black)
    scatter!(ax1,x2, label = L"τ_{ep}", color = :cyan)
    scatter!(ax1,x3, label = L"R_{mv}", color = :yellow)
    scatter!(ax1,x4, label = L"Z_{ao}", color = :hotpink2)
    scatter!(ax1,x5, label = L"R_{s}", color = :tan4)
    scatter!(ax1,x6, label = L"C_{sa}", color = :red1)
    scatter!(ax1,x7, label = L"C_{sv}", color = :silver)
    scatter!(ax1,x8, label = L"E_{max}", color = :blue3)
    scatter!(ax1,x9, label = L"E_{min}", color = :green4)
    fig[3, 1:2] = Legend(fig, ax1, orientation = :horizontal)
    # PP - LV Plot
    x1 = Point2f.(mu[2,1], var[2,1])
    x2 = Point2f.(mu[2,2], var[2,2])
    x3 = Point2f.(mu[2,3], var[2,3])
    x4 = Point2f.(mu[2,4], var[2,4])
    x5 = Point2f.(mu[2,5], var[2,5])
    x6 = Point2f.(mu[2,6], var[2,6])
    x7 = Point2f.(mu[2,7], var[2,7])
    x8 = Point2f.(mu[2,8], var[2,8])
    x9 = Point2f.(mu[2,9], var[2,9])
    scatter!(ax2,x1, label = L"τ_{es}", color = :black)
    scatter!(ax2,x2, label = L"τ_{ep}", color = :cyan)
    scatter!(ax2,x3, label = L"R_{mv}", color = :yellow)
    scatter!(ax2,x4, label = L"Z_{ao}", color = :hotpink2)
    scatter!(ax2,x5, label = L"R_{s}", color = :tan4)
    scatter!(ax2,x6, label = L"C_{sa}", color = :red1)
    scatter!(ax2,x7, label = L"C_{sv}", color = :silver)
    scatter!(ax2,x8, label = L"E_{max}", color = :blue3)
    scatter!(ax2,x9, label = L"E_{min}", color = :green4)
    # PP - SA Plot
    x1 = Point2f.(mu[3,1], var[3,1])
    x2 = Point2f.(mu[3,2], var[3,2])
    x3 = Point2f.(mu[3,3], var[3,3])
    x4 = Point2f.(mu[3,4], var[3,4])
    x5 = Point2f.(mu[3,5], var[3,5])
    x6 = Point2f.(mu[3,6], var[3,6])
    x7 = Point2f.(mu[3,7], var[3,7])
    x8 = Point2f.(mu[3,8], var[3,8])
    x9 = Point2f.(mu[3,9], var[3,9])
    scatter!(ax3,x1, label = L"τ_{es}", color = :black)
    scatter!(ax3,x2, label = L"τ_{ep}", color = :cyan)
    scatter!(ax3,x3, label = L"R_{mv}", color = :yellow)
    scatter!(ax3,x4, label = L"Z_{ao}", color = :hotpink2)
    scatter!(ax3,x5, label = L"R_{s}", color = :tan4)
    scatter!(ax3,x6, label = L"C_{sa}", color = :red1)
    scatter!(ax3,x7, label = L"C_{sv}", color = :silver)
    scatter!(ax3,x8, label = L"E_{max}", color = :blue3)
    scatter!(ax3,x9, label = L"E_{min}", color = :green4)
    # Mean - Qs Plot 
    x1 = Point2f.(mu[4,1], var[4,1])
    x2 = Point2f.(mu[4,2], var[4,2])
    x3 = Point2f.(mu[4,3], var[4,3])
    x4 = Point2f.(mu[4,4], var[4,4])
    x5 = Point2f.(mu[4,5], var[4,5])
    x6 = Point2f.(mu[4,6], var[4,6])
    x7 = Point2f.(mu[4,7], var[4,7])
    x8 = Point2f.(mu[4,8], var[4,8])
    x9 = Point2f.(mu[4,9], var[4,9])
    scatter!(ax4,x1, label = L"τ_{es}", color = :black)
    scatter!(ax4,x2, label = L"τ_{ep}", color = :cyan)
    scatter!(ax4,x3, label = L"R_{mv}", color = :yellow)
    scatter!(ax4,x4, label = L"Z_{ao}", color = :hotpink2)
    scatter!(ax4,x5, label = L"R_{s}", color = :tan4)
    scatter!(ax4,x6, label = L"C_{sa}", color = :red1)
    scatter!(ax4,x7, label = L"C_{sv}", color = :silver)
    scatter!(ax4,x8, label = L"E_{max}", color = :blue3)
    scatter!(ax4,x9, label = L"E_{min}", color = :green4)

    Label(fig[1, 1, TopLeft()], "A",fontsize = 18,font = :bold,halign = :right)
    Label(fig[2, 1, TopLeft()], "B",fontsize = 18,font = :bold,halign = :right)
    Label(fig[1, 2, TopLeft()], "C",fontsize = 18,font = :bold,halign = :right)
    Label(fig[2, 2, TopLeft()], "D",fontsize = 18,font = :bold,halign = :right)
    fig
end

# Calculate the Orthogonality Matrices 
Orth_heat = Matrix{Float64}(undef,9,9)
for j in 1:9
    for i in 1:9
        if i==j
            Orth_heat[i,j] = 0
        else 
        Orth_heat[i,j] = sin(acos(((transpose(mu[:,i])*mu[:,j]))/(norm(mu[:,i])*norm(mu[:,j]))-1e-15))  #Slight numerical rounding error without the additional add on 
        end 
    end 
end 

#save("/home/harry/Desktop/PhD/Year 3/Li rewrite/Orthogonality/Morris_mu.jld","data",Orth_heat)

# Calculate the average orth for each parameter 
rank = zeros(9)
for j in 1:9
    rank[j] = sum(Orth_heat[:,j])/8
end 
p=sortperm(rank,rev=true)
using CairoMakie
f = Figure();
ax = Axis(f[1,1],title="Morris Sensitivity Matrix-Parameter Orthogonality Rank", xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"][p]), xlabel = "Parameters", ylabel = "Orthogonality")
CairoMakie.scatter!( rank[p])
f



# Calculate overall parameter influence Ej 
F = transpose(mu)*mu
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
ax = Axis(f[1,1],title="Morris Sensitivity Matrix-Parameter Importance PCA Method", xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"][p]), xlabel = "Parameters", ylabel = "Importance")
CairoMakie.scatter!( e[p])
f

S = Array(transpose(mu))
using MAT
file = matopen("morris.mat", "w")
write(file, "morris", S)
close(file)