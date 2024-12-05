include("./functions.jl")

## Inputs
dumpfilepath = "/mnt/iusers01/ceas01/f21157am/scratch/water_Nov2024/dump.lammpstrj"
timestep = 0.7; #femptosec

num_lines::Int32 = countlines(dumpfilepath)
natoms = CSV.File(dumpfilepath, skipto=4, limit=1, header=false).Column1[1]
nhydrogens::Int32 = 2 * natoms / 3
totalsteps = (floor(Int, num_lines / (natoms + 9)) ) 
#totalsteps = div(totalsteps,2)

# create time array, and convert to seconds from fs
t = collect(Float32, 0:totalsteps-1) .* (50 * timestep) .* 1e-15;

prefactor = (3 / 16) * (μ₀/(4π))^2 * ħ^2 * γ^4

# Intramolecular contributions
F = calculateF(dumpfilepath, "intra", totalsteps);
G_R_ens_av = mean(ACF.(eachrow(F)))
G_R = prefactor * G_R_ens_av * 1e60 ;
τ_R = (1 / G_R[1]) * trapz(t, G_R)
Δω²_R = 3 * prefactor * mean(mean(F.^2, dims=2)) * 1e60
Δω_R = sqrt(Δω²_R)/2π /1000 # KHz
Tr = 1 / (10 / 3 * Δω²_R * τ_R )

# intermolecular contributions
F = calculateF(dumpfilepath, "inter",totalsteps);
G_T_ens_av = mean(ACF.(eachrow(F)))
G_T = prefactor * G_T_ens_av * 1e60 ;
τ_T = (1 / G_T[1]) * trapz(t, G_T) 
Δω²_T = 3 * prefactor * mean(mean(F.^2, dims=2)) * 1e60
Δω_T = sqrt(Δω²_T)/2π /1000 # KHz
Tt = 1 / (10 / 3 * Δω²_T * τ_T )

# all contributions
F = calculateF(dumpfilepath, "all",totalsteps);
G_ens_av = mean(ACF.(eachrow(F)))
G = prefactor * G_ens_av * 1e60 ;
τ = (1 / G[1]) * trapz(t, G_T) 
Δω² = 3 * prefactor * mean(mean(F.^2, dims=2)) * 1e60
Δω = sqrt(Δω²)/2π /1000 # KHz
Tt = 1 / (10 / 3 * Δω² * τ )

# Relaxation times
T = 1 / ((10 / 3) * Δω²_R * τ_R + (10 / 3) * Δω²_T * τ_T)


J(G, t, ω) = 2 * trapz(t, (G .* cos.(ω .* t)))
ω₀ = γ * 1 # (rad s^-1)

T1_R = 1/ (J(G_R, t, ω₀) + 4J(G_R, t, 2ω₀))
T2_R = 1/ ((3/2)*J(G_R, t, 0) + (5/2)*J(G_R, t, ω₀) + J(G_R, t, 2ω₀))

T1_T = 1/ (J(G_T, t, ω₀) + 4J(G_T, t, 2ω₀))
T2_T = 1/ ((3/2)*J(G_T, t, 0) + (5/2)*J(G_T, t, ω₀) + J(G_T, t, 2ω₀))

T1 = 1/ (J(G, t, ω₀) + 4J(G, t, 2ω₀))
T2 = 1/ ((3/2)*J(G, t, 0) + (5/2)*J(G, t, ω₀) + J(G, t, 2ω₀))

println("Intramolecular")
@show τ_R
@show Δω_R 
@show T1_R
@show T2_R
println("")

println("Intermolecular")
@show τ_T
@show Δω_T 
@show T1_T
@show T2_T
println("")

println("all")
@show τ
@show Δω²
@show T1_T
@show T2_T

using DelimitedFiles
writedlm("./intra_correlation_data.txt", [t G_R])
writedlm("./inter_correlation_data.txt", [t G_T])
