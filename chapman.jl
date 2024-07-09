## Chapman paper method


nblocks = 3
G = mainloop(900, nblocks, "intra")
Gᵣ =  real.(Base.stack(vec(mean(reshape(G, size(G)[1] ÷ nblocks, nblocks), dims=2)))[1, :])
Δω²ᵣ = 3* 3 / 16 * (μ₀ / 4π)^2 * γ^4 * ħ^2 * Gᵣ[1]
Jᵣ = 2 * real.(fft(Gᵣ))
Jᵣ =  3 / 16 * (μ₀ / 4π)^2 * γ^4 * ħ^2 *Jᵣ
τᵣ = real(1 / 2 * (Jᵣ[1] / Gᵣ[1]))

G = mainloop(900, nblocks, "inter")
Gₜ = real.(Base.stack(vec(mean(reshape(G, size(G)[1] ÷ nblocks, nblocks), dims=2)))[1, :])
Δω²ₜ = 3 * 3 / 16 * (μ₀ / 4π)^2 * γ^4 * ħ^2 * Gₜ[1]
Jₜ = 2 * real.(fft(Gₜ))
Jₜ =  3 / 16 * (μ₀ / 4π)^2 * γ^4 * ħ^2 *Jₜ
τₜ = real(1 / 2 * (Jₜ[1] / Gₜ[1]))

T = 1 / (10 / 3 * Δω²ₜ * τₜ + 10 / 3 * Δω²ᵣ * τᵣ)

plot(Gₜ./Gₜ[1])


plot(Gᵣ)