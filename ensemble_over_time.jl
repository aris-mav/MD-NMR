
## Ensemble averaging over time

F = calculateF("intra")

tp = 200 # Number of points in the t dimension

ensemble_τt = ones(size(F, 2) - tp, tp) # Initialise array
ensemble_ijt = ones(size(F, 1), tp)

for i in axes(F, 1) # For every pair of hydrogens

    for τ in 1:(size(F, 2)-tp) # For every starting time τ
        ensemble_τt[τ, :] = F[i, τ] .* F[i, τ:(τ+tp-1)] # calculate correlations and put them in a matrix
    end

    ensemble_ijt[i, :] = sum(ensemble_τt, dims=1) ./ size(ensemble_τt, 1) # average the τ dimension of that matrix

end

G = 3 / 16 * (μ₀ / 4π)^2 * ħ^2 * γ^4 * vec(sum(ensemble_ijt, dims=1) ./ size(ensemble_ijt, 1)) 

J = 2 * abs.(fft(G))

τᵣₜ = J[1] / (2 * G[1])

plot(t[1:length(G)] ./ τᵣₜ , G ./ G[1])

