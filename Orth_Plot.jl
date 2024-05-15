using JLD, CairoMakie
CairoMakie.activate!(type = "svg")

Local= load("path")["data"]
S1_E = load("path")["data"]
ST_E = load("path")["data"]
mu_M = load("path")["data"]
S1_S = load("path")["data"]
ST_S = load("path")["data"]

# Orth Hist Data
begin
    # Local 
    a1 = Local[2:end,1]
    a2 = Local[3:end,2]
    a3 = Local[4:end,3]
    a4 = Local[5:end,4]
    a5 = Local[6:end,5]
    a6 = Local[7:end,6]
    a7 = Local[8:end,7]
    a8 = Local[9:end,8]
    a_Local = reduce(vcat, (a1,a2,a3,a4,a5,a6,a7,a8))
    # Morris
    a1 = mu_M[2:end,1]
    a2 = mu_M[3:end,2]
    a3 = mu_M[4:end,3]
    a4 = mu_M[5:end,4]
    a5 = mu_M[6:end,5]
    a6 = mu_M[7:end,6]
    a7 = mu_M[8:end,7]
    a8 = mu_M[9:end,8]
    a_muM = reduce(vcat, (a1,a2,a3,a4,a5,a6,a7,a8))
    # eFAST S1
    a1 = S1_E[2:end,1]
    a2 = S1_E[3:end,2]
    a3 = S1_E[4:end,3]
    a4 = S1_E[5:end,4]
    a5 = S1_E[6:end,5]
    a6 = S1_E[7:end,6]
    a7 = S1_E[8:end,7]
    a8 = S1_E[9:end,8]
    a_S1E = reduce(vcat, (a1,a2,a3,a4,a5,a6,a7,a8))
    # eFAST ST
    a1 = ST_E[2:end,1]
    a2 = ST_E[3:end,2]
    a3 = ST_E[4:end,3]
    a4 = ST_E[5:end,4]
    a5 = ST_E[6:end,5]
    a6 = ST_E[7:end,6]
    a7 = ST_E[8:end,7]
    a8 = ST_E[9:end,8]
    a_STE = reduce(vcat, (a1,a2,a3,a4,a5,a6,a7,a8))
    # Sobol S1
    a1 = S1_S[2:end,1]
    a2 = S1_S[3:end,2]
    a3 = S1_S[4:end,3]
    a4 = S1_S[5:end,4]
    a5 = S1_S[6:end,5]
    a6 = S1_S[7:end,6]
    a7 = S1_S[8:end,7]
    a8 = S1_S[9:end,8]
    a_S1S = reduce(vcat, (a1,a2,a3,a4,a5,a6,a7,a8))
    # Sobol ST
    a1 = ST_S[2:end,1]
    a2 = ST_S[3:end,2]
    a3 = ST_S[4:end,3]
    a4 = ST_S[5:end,4]
    a5 = ST_S[6:end,5]
    a6 = ST_S[7:end,6]
    a7 = ST_S[8:end,7]
    a8 = ST_S[9:end,8]
    a_STS = reduce(vcat, (a1,a2,a3,a4,a5,a6,a7,a8))
end 


begin
    f = Figure(size = (1400, 1000));

    ax1 = Axis(f[1,1], xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), yticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), title = "Local Orthogonality Matrix", xlabel = "Parameters", ylabel = "Parameters")
    hm1 = CairoMakie.heatmap!(ax1,Local, colormap=:plasma)
    for i in 1:9, j in 1:9
        txtcolor = Local[i, j] < -0.0 ? :white : :black
        text!(ax1, "$(round(Local[i,j], digits = 2))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 15)
    end
    hidexdecorations!(ax1)

    ax2 = Axis(f[2,1], xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), yticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), title = "Morris Orthogonality Matrix", xlabel = "Parameters", ylabel = "Parameters")
    hm2 = CairoMakie.heatmap!(ax2,mu_M, colormap=:plasma)
    for i in 1:9, j in 1:9
        txtcolor = mu_M[i, j] < -0.0 ? :white : :black
        text!(ax2, "$(round(mu_M[i,j], digits = 2))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 15)
    end

    ax3 = Axis(f[1,2], xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), yticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), title = "eFast S1 Orthogonality Matrix", xlabel = "Parameters", ylabel = "Parameters")
    hm3 = CairoMakie.heatmap!(ax3,S1_E, colormap=:plasma)
    for i in 1:9, j in 1:9
        txtcolor = S1_E[i, j] < -0.0 ? :white : :black
        text!(ax3, "$(round(S1_E[i,j], digits = 2))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 15)
    end
    hidedecorations!(ax3)

    ax4 = Axis(f[2,2], xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), yticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), title = "eFAST ST Orthogonality Matrix", xlabel = "Parameters", ylabel = "Parameters")
    hm4 = CairoMakie.heatmap!(ax4,ST_E, colormap=:plasma)
    for i in 1:9, j in 1:9
        txtcolor = ST_E[i, j] < -0.0 ? :white : :black
        text!(ax4, "$(round(ST_E[i,j], digits = 2))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 15)
    end
    hideydecorations!(ax4)

    ax5 = Axis(f[1,3], xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), yticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), title = "Sobol S1 Orthogonality Matrix", xlabel = "Parameters", ylabel = "Parameters")
    hm5 = CairoMakie.heatmap!(ax5,S1_S, colormap=:plasma)
    for i in 1:9, j in 1:9
        txtcolor = S1_S[i, j] < -0.0 ? :white : :black
        text!(ax5, "$(round(S1_S[i,j], digits = 2))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 15)
    end
    hidedecorations!(ax5)

    ax6 = Axis(f[2,3], xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), yticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), title = "Sobol ST Orthogonality Matrix", xlabel = "Parameters", ylabel = "Parameters")
    hm6 = CairoMakie.heatmap!(ax6,ST_S, colormap=:plasma)
    for i in 1:9, j in 1:9
        txtcolor = ST_S[i, j] < -0.0 ? :white : :black
        text!(ax6, "$(round(ST_S[i,j], digits = 2))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 15)
    end
    hideydecorations!(ax6)

    ax7 = Axis(f[3,1], xticks = 0.0:0.1:1.0, title = "Local Orthogonality Histogram", xlabel = "Orthogonality Score", ylabel = "Density")
    ax8 = Axis(f[4,1], xticks = 0.0:0.1:1.0, title = "Morris Orthogonality Histogram", xlabel = "Orthogonality Score", ylabel = "Density")
    ax9 = Axis(f[3,2], xticks = 0.0:0.1:1.0, title = "eFAST S1 Orthogonality Histogram", xlabel = "Orthogonality Score", ylabel = "Density")
    ax10 = Axis(f[4,2], xticks = 0.0:0.1:1.0, title = "eFAST ST Orthogonality Histogram", xlabel = "Orthogonality Score", ylabel = "Density")
    ax11 = Axis(f[3,3], xticks = 0.0:0.1:1.0, title = "Sobol S1 Orthogonality Histogram", xlabel = "Orthogonality Score", ylabel = "Density")
    ax12 = Axis(f[4,3], xticks = 0.0:0.1:1.0, title = "Sobol ST Orthogonality Histogram", xlabel = "Orthogonality Score", ylabel = "Density")
    linkaxes!(ax7,ax8,ax9,ax10,ax11,ax12)
    hidexdecorations!(ax7, grid = false)
    hidedecorations!(ax9, grid = false)
    hideydecorations!(ax10, grid = false)
    hidedecorations!(ax11, grid = false)
    hideydecorations!(ax12, grid = false)

    hist!(ax7, a_Local, bins = 0.0:0.1:1.0, color = :purple,strokewidth = 1,strokecolor = :black , normalization = :none)
    hist!(ax8, a_muM,   bins = 0.0:0.1:1.0,   color = :purple, strokewidth = 1, strokecolor = :black , normalization = :none)
    hist!(ax9, a_S1E,   bins = 0.0:0.1:1.0,   color = :purple, strokewidth = 1, strokecolor = :black , normalization = :none)
    hist!(ax10, a_STE,  bins = 0.0:0.1:1.0,  color = :purple, strokewidth = 1, strokecolor = :black, normalization = :none)
    hist!(ax11, a_S1S,  bins = 0.0:0.1:1.0,  color = :purple, strokewidth = 1, strokecolor = :black, normalization = :none)
    hist!(ax12, a_STS,  bins = 0.0:0.1:1.0,  color = :purple, strokewidth = 1, strokecolor = :black, normalization = :none)

    Label(f[1, 1, TopLeft()], "A",fontsize = 18,font = :bold,halign = :right)
    Label(f[2, 1, TopLeft()], "D",fontsize = 18,font = :bold,halign = :right)
    Label(f[1, 2, TopLeft()], "B",fontsize = 18,font = :bold,halign = :right)
    Label(f[2, 2, TopLeft()], "E",fontsize = 18,font = :bold,halign = :right)
    Label(f[1, 3, TopLeft()], "C",fontsize = 18,font = :bold,halign = :right)
    Label(f[2, 3, TopLeft()], "F",fontsize = 18,font = :bold,halign = :right)
    Label(f[3, 1, TopLeft()], "G",fontsize = 18,font = :bold,halign = :right)
    Label(f[4, 1, TopLeft()], "J",fontsize = 18,font = :bold,halign = :right)
    Label(f[3, 2, TopLeft()], "H",fontsize = 18,font = :bold,halign = :right)
    Label(f[4, 2, TopLeft()], "K",fontsize = 18,font = :bold,halign = :right)
    Label(f[3, 3, TopLeft()], "I",fontsize = 18,font = :bold,halign = :right)
    Label(f[4, 3, TopLeft()], "L",fontsize = 18,font = :bold,halign = :right)


    CairoMakie.Colorbar(f[1:2,4],hm1, label = "Orthogonality Score", ticks = 0.0:0.2:1.0)
    f
end 













S1 = S1_S
ST = ST_S
# Plotting just Sobol Orthogonality 
using CairoMakie
CairoMakie.activate!(type = "svg")
begin
    f = Figure(size = (1100, 600), backgroundcolor = RGBf(0.98, 0.98, 0.98))
    #f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98));

    ax1 = Axis(f[1,1], xticklabelrotation = π / 3, xticklabelalign = (:right, :center), xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), yticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), title = "First Order Orthogonality Matrix", ylabel = "Parameters", xlabel = "Parameters")
    Label(f[1, 1, TopLeft()], "A",fontsize = 18,font = :bold,halign = :right)
    hm1 = CairoMakie.heatmap!(ax1,transpose(S1), colormap=:seaborn_pastel6)
    for i in 1:9, j in 1:9
        txtcolor = transpose(S1)[i, j] < -1000.0 ? :white : :black
        text!(ax1, "$(round(transpose(S1)[i,j], digits = 4))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 12)
    end


    ax2 = Axis(f[1,3], xticklabelrotation = π / 3, xticklabelalign = (:right, :center), xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), title = "Total Order Orthogonality Matrix",  xlabel = "Parameters")
    hideydecorations!(ax2)
    Label(f[1, 3, TopLeft()], "B",fontsize = 18,font = :bold,halign = :right)
    hm2 = CairoMakie.heatmap!(ax2,transpose(ST), colormap=:seaborn_pastel6)
    for i in 1:9, j in 1:9
        txtcolor = transpose(ST)[i, j] < -1000.0 ? :white : :black
        text!(ax2, "$(round(transpose(ST)[i,j], digits = 4))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 12)
    end
    Colorbar(f[1,2],hm2)

    gl = GridLayout(f[2, 1:3], width = Relative(0.5))
    A = abs.(transpose(ST - S1))
    ax3 = Axis(gl[1, 1], xticklabelrotation = π / 3, xticklabelalign = (:right, :center), xticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), yticks = (1:9, [L"\tau_{es}", L"\tau_{ep}", L"R_{mv}", L"Z_{ao}", L"R_{s}", L"C_{sa}", L"C_{sv}", L"E_{max}", L"E_{min}"]), title = "Difference Orthogonality Matrix", ylabel = "Parameters",  xlabel = "Parameters")
    Label(gl[1, 1, TopLeft()], "C", fontsize = 18,font = :bold,halign = :right)
    hm3 = CairoMakie.heatmap!(ax3,A, colormap=:seaborn_pastel6)
    for i in 1:9, j in 1:9
        txtcolor = A[i, j] < -1000.0 ? :white : :black
        text!(ax3, "$(round(A[i,j], digits = 4))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 12)
    end
    Colorbar(gl[1,2],hm3)
    linkaxes!(ax1,ax2,ax3)
    f
end 
