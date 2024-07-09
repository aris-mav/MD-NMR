# Fmatrix = calc_F(100)
# G = calc_G(Fmatrix)

G = mainloop(900, 3, "all")

plot(t,real.(first.(G)))


## Block averaging
Gmean = Base.stack(vec(mean(reshape(G, size(G)[1] ÷ nblocks, nblocks), dims=2)))
# Gmean = Base.stack(vec(mean(reshape(G, size(G)[1]), dims=2)))

# Add prefactors (tests to match paper values)
# Gmean = Gmean * 0.0265 ./ Gmean[1, 1]
# Gmean = (13/6)*(μ₀/4pi)^2 .*Gmean

#save the plot above
# savefig("../Correlation_functions.png")

J = [2 .* real.(fft(Gmean[i, :])) for i in 1:3];

ωrange = fftfreq(length(t[1:size(Gmean)[2]]), 1 / (t[2] - t[1])) .* 1e15 * 2π; # 1e15 converts to Hz (fs to s) and 2π converts to rad/s
ω₀ = γ * 1;# rad/s
ω₀index = min(findmin(abs.(ωrange .- ω₀))[2])

# Calculate T1 and T2
T1 = 1 / (((9 / 8) * (μ₀ / 4π)^2 * γ^4 * ħ^2 * (J[2][ω₀index] + J[3][ω₀index])))
T2 = 1 / ((9 / 24) * (μ₀ / 4π)^2 * γ^4 * ħ^2 * (J[1][ω₀index] + 10J[2][ω₀index] + J[3][ω₀index]))

display("T1 = $T1 s")
display("T2 = $T2 s")



## Replicating mathematica Carles method

#G = mainloop(900, 3, "inter")
#Gmean_inter = Base.stack(vec(mean(reshape(G, size(G)[1] ÷ nblocks, nblocks), dims=2)))
#Jinter = [2 .* real.(fft(Gmean_inter[i, :])) for i in 1:3];
#plotG(Gmean_inter)
#
#G = mainloop(900, 3, "intra")
#Gmean_intra = Base.stack(vec(mean(reshape(G, size(G)[1] ÷ nblocks, nblocks), dims=2)))
#Jintra = [2 .* real.(fft(Gmean_intra[i, :])) for i in 1:3];
#plotG(Gmean_intra)
#
#T₁ = 1 / ((9 / 8 * (μ₀ / 4π)^2 * γ^4 * ħ^2 * 2 * (5.08 * (Jinter[2][1] + Jinter[3][1]) + Jintra[2][1] + Jintra[3][1])) * 10^(60))
#T₂ = 1 / ((3 / 4 * (μ₀ / 4π)^2 * γ^4 * ħ^2 * 2 * (5.08 * (3 / 8 * Jinter[1][1] + 15 / 4 * Jinter[2][1] + 3 / 8 * Jinter[3][1])) + 3 / 8 * Jintra[1][1] + 15 / 4 * Jintra[2][1] + 3 / 8 * Jintra[3][1]) * 10^(60 - 15))
