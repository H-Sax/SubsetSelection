using JLD, GlobalSensitivity, LinearAlgebra, DelimitedFiles
using MAT # allows reading .mat matrices files from MatLab
using Combinatorics

# read in Sensitivity Matrix from same directory
filename = "FirstPython.mat"
filepath = joinpath(@__DIR__, filename)
println("Reading Sensitivity matrix from: " * filepath)
vars = matopen(filepath)
S = read(vars, "First")
close(vars)

# (12) calculates E vector - difficulty in estimating each parameter
function effect_E(S)
    # compute covariance matrix
    X = S' * S

    # principal components; eigen decomposition
    eigens = eigen(X)
    λ, C = eigens.values, eigens.vectors
    
    # (12) calculate E measure of each parameter - multiplies each column C_j by λ, sums
    E = sum(abs.(C.*λ), dims=1) / sum(λ)

    # convert E from 1d Array to vector and return
    return vec(E)
end

# (13-16) calculates the linear dependence metric dⱼ
# where sⱼ is the sensitivity vector of chosen remaining parameter θⱼ
# and Sₚ is the set of sₖ sensitivity vectors of chosen parameters P
function dependence_dⱼ(sⱼ, Sₚ)
    # (15) left hand side matrix
    SₚᵀS = Sₚ' * Sₚ
    # (15) right hand side (using (sⱼᵀSₚ)ᵀ = Sₚᵀsⱼ )
    sⱼᵀSₚ = Sₚ' * sⱼ

    # solve (15) for a
    a = inv(SₚᵀS) * sⱼᵀSₚ

    # (13) compute minimum-distance sensitivity vector s
    s = (sum(a'.*Sₚ, dims=2))

    # (16) linear independence metric
    nom = (sⱼ' * s)[1];
    denom = norm(sⱼ) * norm(s);
    dⱼ = sin(acos(min.(nom./denom, 1))); # min used to avoid floating point errors acos(>1)

    return dⱼ
end

# Step 3) calculate minimum d over combinations of selected parameters
# where Sθ is the set of remaining parameter sensitivity vectors
# and Sₚ is the set of sₖ sensitivity vectors of chosen parameters P
function tuples_dⱼ(Sθ, Sₚ, m)
    # number of parameters chosen
    n = size(Sₚ)[2]
    # generate list of indices of Sₚ to permute
    idx = collect(1:n)
    # all n choose (m-1) combinations of indices
    combinations_Sₚ = combinations(idx, m-1)

    # initialise array of dⱼ metrics, set to 1 (max range)
    d = Array{Float64}(undef, size(Sθ)[2])
    d = fill!(d, 1)

    # TODO CHECK USE VIEWS HERE FOR TUPLES?
    for c in combinations_Sₚ
        # get dⱼ for each remaining sⱼ in Sₚ w.r.t current tuple of Sθ 'c'
        d_new = [dependence_dⱼ(sⱼ, Sₚ[:,c]) for sⱼ in eachcol(Sθ)]

        # indices of new dⱼ that are less than the current minimum
        lessers = d.>=d_new
        # update d vector with all lesser values
        d[lessers] .= d_new[lessers]
    end
    
    return d
end

# Implements Li parameter selection procedure (Section III.C)
function Li(S)
    # Determine number of measurements and parameters 
    p, m = size(S)
    # Transpose S to fit Li paper convention
    S = S'

    # stores chosen params in order of greatest rank
    Prank = Vector{Integer}(undef, p)
    # stores selected params as a boolean mask, to use or exclude parameters
    P = Vector{Bool}(undef, p)
    P = fill!(P, 0)
    # range 1 to p, to help identify actual index of param
    idx = collect(1:p)

    # Step 1; calculate overall effect of params 
    E = effect_E(S)
    # select S with highest E
    p_idx = argmax(E)
    n = 1
    
    Prank[n] = p_idx
    P[p_idx] = true

    # Steps 2 - 4; select remaining parameters
    for i in 2:p
        # CHECK IF VIEWS SLOWER? Non-contiguous in CPU
        # TODO possible improvement; use staticarray when p < 100
        selected = @view S[:,P]
        unselected = @view S[:,.!P]

        if i < m
            # Step 2 (case n < m): calculate dⱼ for each unselected parameter
            d = [dependence_dⱼ(j, selected) for j in eachcol(unselected)]
        else
            # Step 3 (case n>= m): calculate minimum d over combinations of selected parameters
            d = tuples_dⱼ(unselected, selected, m)
        end

        # Step 4: select parameter with highest identifiability index
        # calculate identifiability index for unselected params

        I = E[.!P] .* d
        # select index of greatest identifiability parameter
        max = argmax(I)
        # looks up corresponding original index from index in reduced set
        true_idx = (idx[.!P][max])
        # update ranking list and mask
        n += 1
        Prank[n] = true_idx
        P[true_idx] = true
    end

    return Prank
end

@time begin
    test = Li(S)
    println(test)
end

# used to check speed in REPL without compilation time
function run(S)
    @time begin
        test = Li(S)
    end
end