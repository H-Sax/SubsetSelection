using DifferentialEquations, Statistics, FiniteDiff, LinearAlgebra

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


function circ_local(p)
    _prob = ODEProblem(NIK!, [6.0, 6.0, 6.0, 200.0, 0.0, 0.0, 0.0], (0.0,10.0), p)
    newsol = solve(_prob, Vern7(),  reltol = 1e-11, abstol = 1e-11, saveat = x)
    [(maximum(newsol[4,:]) - minimum(newsol[4,:])), (maximum(newsol[1,:]) - minimum(newsol[1,:])), (maximum(newsol[2,:]) - minimum(newsol[2,:])), mean(newsol[7,:])]
end
# SV LV, PP LV, PP SA, Mean systemic flow 

## Absoloute sensititivty
s = FiniteDiff.finite_difference_jacobian(circ_local, [0.3, 0.45, 0.006, 0.033, 1.11, 1.13, 11.0, 1.5, 0.03])
## Make relative by scaling by p/y 
@time sol = solve(prob, Vern7(),  reltol = 1e-11, abstol = 1e-11, saveat = x)
measurements = [(maximum(sol[4,:]) - minimum(sol[4,:])), (maximum(sol[1,:]) - minimum(sol[1,:])), (maximum(sol[2,:]) - minimum(sol[2,:])), mean(sol[7,:])]

S = Matrix{Float64}(undef,4,9)
for i in 1:4
    for j in 1:9
        S[i,j] = (prob.p[j]/measurements[i]) * s[i,j]
    end 
end

# Relative normalised sensititivty matrix 
S = abs.(S)


#S = load("/home/harry/Desktop/PhD/Year 3/Li rewrite/Local/Local.jld")["Data"]
# Plot the normalised relative sensitivity Matrix
using CairoMakie
CairoMakie.activate!(type = "svg")
## Relative Sensitivity Matrix 
f = Figure(size = (800, 220));
# to add axis we specify a position in the figures layout as first argument f[1,1] to fill the whole plot 
ax = Axis(f[1,1], xticklabelrotation = π / 3, xticklabelalign = (:right, :center), xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), yticks = (1:4, ["SV - LV", "PP - LV", "PP - SA", "Mean - Qs"]), title = "Relative Sensitivity Matrix", ylabel = "Measurements", xlabel = "Parameters")
hm = CairoMakie.heatmap!(ax,transpose(S), colormap=:plasma)
for i in 1:9, j in 1:4
    txtcolor = transpose(S)[i, j] < -1000.0 ? :white : :black
    text!(ax, "$(round(transpose(S)[i,j], digits = 4))", position = (i, j),
        color = txtcolor, align = (:center, :center), fontsize = 12)
end
Colorbar(f[1,2],hm)
f 

# Calculate the Orthogonality Matrices 
Orth_heat = Matrix{Float64}(undef,9,9)
for j in 1:9
    for i in 1:9
        if i==j
            Orth_heat[i,j] = 0
        else 
        Orth_heat[i,j] = sin(acos(((transpose(S[:,i])*S[:,j]))/(norm(S[:,i])*norm(S[:,j]))-1e-15))  #Slight numerical rounding error without the additional add on 
        end 
    end 
end 

save("/home/harry/Desktop/PhD/Year 3/Li rewrite/Orthogonality/Local_S.jld","data",Orth_heat)

# Calculate the average orth for each parameter 
rank = zeros(9)
for j in 1:9
    rank[j] = sum(Orth_heat[:,j])/8
end 
p=sortperm(rank,rev=true)
using CairoMakie
f = Figure();
ax = Axis(f[1,1],title="Local Sensitivity Matrix-Parameter Orthogonality Rank", xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"][p]), xlabel = "Parameters", ylabel = "Orthogonality")
CairoMakie.scatter!( rank[p])
f



# Calculate overall parameter influence Ej 
F = transpose(S)*S
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
ax = Axis(f[1,1],title="Local Sensitivity Matrix-Parameter Importance PCA Method", xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"][p]), xlabel = "Parameters", ylabel = "Importance")
CairoMakie.scatter!( e[p])
f

S = Array(transpose(S))
using MAT
file = matopen("Local.mat", "w")
write(file, "local", S)
close(file)