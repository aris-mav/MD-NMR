include("./functions.jl")
# unicodeplots()

## Inputs
# dumpfilepath = "/lustre03/other/9847zg/dump_files/dump.lammpstrj";
dumpfilepath = "/mnt/c/Users/arism/OneDrive - The University of Manchester/Codes/MD_to_NMR/dump_files/dump.lammpstrj";
timestep = 1.5; #femptosec
nblocks = 1;
num_lines::Int32 = countlines(dumpfilepath)
natoms = CSV.File(dumpfilepath, skipto=4, limit=1, header=false).Column1[1]
nhydrogens::Int32 = 2 * natoms / 3
totalsteps = (floor(Int, num_lines / (natoms + 9)) ÷ nblocks) * nblocks

# preallocate  t (time array), and convert to seconds from fs
t = collect(Float32, 0:totalsteps-1) .* (100 * timestep) .* 1e-15

J(ω) = 2 * trapz(t, (G .* cos.(ω .* t)))

ω₀ = γ * 1 # (rad s^-1)
prefactor = (3 / 16) * (μ₀ / (4π))^2 * ħ^2 * γ^4

# Intramolecular contributions
F = calculateF(dumpfilepath, "intra",totalsteps)
G = prefactor * mean(ACF.(eachrow(F))) * 10^60
τ = (1 / G[1]) * trapz(t, G)

Δω = sqrt( 3 * G[1] )/2π

T1 = 1/ (J(ω₀) + 4J(2ω₀))
T2 = 1/ ((3/2)*J(0) + (5/2)*J(ω₀) + J(2ω₀))

# intermolecular contributions
F = calculateF(dumpfilepath, "inter",totalsteps)
G = prefactor * mean(ACF.(eachrow(F))) * 10^60
τ = (1 / G[1]) * trapz(t, G)