include("./functions.jl")

## Inputs
dumpfilepath = "../dump_files/dump.lammpstrj"
timestep = 1.5; #femptosec
ω₀ = γ * 1 # (rad s^-1)

t = time_array(dumpfilepath, timestep)
prefactor = (3 / 16) * (μ₀/(4π))^2 * ħ^2 * γ^4

for contributions in ["intra", "inter"]

    F = calculateF(dumpfilepath, contributions);
    G_ens_av = mean(ACF.(eachrow(F)))
    G = prefactor * G_ens_av * 1e60 ;
    τ = (1 / G[1]) * trapz(t, G)
    Δω² = 3 * prefactor * mean(mean(F.^2, dims=2)) * 1e60
    Δω = sqrt(Δω²)/2π /1000 # KHz
    T = 1 / (10 / 3 * Δω² * τ )

    T1 = 1/ (J(G, t, ω₀) + 4J(G, t, 2ω₀))
    T2 = 1/ ((3/2)*J(G, t, 0) + (5/2)*J(G, t, ω₀) + J(G, t, 2ω₀))

    println("")
    println("Results for "*contributions*"-molecular contributions")
    @show τ
    @show Δω²
    @show Δω 
    @show T1
    @show T2
    println("")
    writedlm("./"*contributions*"_correlation_data.txt", [t G_ens_av])

end

open(dumpfilepath) do io

    boxlengths = zeros(3)
    # Loop over time steps

        # Skip headers
        for _ in 1:5
            readline(io)
        end

        # Read box bounds
        for i in 1:3
            boxlengths[i] = sum(abs.(parse.(Float64, split(readline(io), ' '))))
        end
        println(boxlengths)
end
