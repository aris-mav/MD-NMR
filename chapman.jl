## Chapman paper method
nblocks = 3
G = mainloop(900, nblocks, "intra")
Gᵣ = 3 / 16 * (μ₀ / 4π)^2 * γ^4 * ħ^2 * real.(Base.stack(vec(mean(reshape(G, size(G)[1] ÷ nblocks, nblocks), dims=2)))[1, :])
Δω²ᵣ = 3 * Gᵣ[1]
Jᵣ = 2 * real.(fft(Gᵣ))
τᵣ = real(1 / 2 * (Jᵣ[1] / Gᵣ[1]))

G = mainloop(900, nblocks, "inter")
Gₜ = 3 / 16 * (μ₀ / 4π)^2 * γ^4 * ħ^2 * real.(Base.stack(vec(mean(reshape(G, size(G)[1] ÷ nblocks, nblocks), dims=2)))[1, :])
Δω²ₜ = 3 * Gₜ[1]
Jₜ = 2 * real.(fft(Gₜ))
τₜ = real(1 / 2 * (Jₜ[1] / Gₜ[1]))

T = 1 / (10 / 3 * Δω²ₜ * τₜ + 10 / 3 * Δω²ᵣ * τᵣ)
