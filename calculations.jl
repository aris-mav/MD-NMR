include("./functions.jl")

## Inputs
dumpfilepath = "dump.lammpstrj"
timestep = 1.5; #femptosec
ω₀ = γ * 1 # (rad s^-1) at 1 tesla

t = time_array(dumpfilepath, timestep)

for x in ["intra", "inter"]

    F = calculateF(dumpfilepath, x);

    G_ens_av = mean(ACF.(eachrow(F))) # Ensemble average (no prefactors)

    G = prefactor * G_ens_av * 1e60 

    τ = (1 / G[1]) * trapz(t, G) # Correlation time (s)

    Δω² = delta_omega(F)
    Δω = sqrt(Δω²)/2π /1000 # KHz
    
    T = 1/ (10 / 3 * Δω² * τ )
    T1 = 1/ (J(G, t, ω₀) + 4J(G, t, 2ω₀))
    T2 = 1/ ((3/2)*J(G, t, 0) + (5/2)*J(G, t, ω₀) + J(G, t, 2ω₀))

    println("")
    println("Results for "*x*"molecular contributions")
    @show τ
    @show Δω²
    @show Δω 
    @show T
    @show T1
    @show T2
    println("")

    writedlm("./"*x*"_correlation_data.txt", [t G_ens_av])
end

