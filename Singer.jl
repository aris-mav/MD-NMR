include("./functions.jl")
# unicodeplots()

## Inputs
# dumpfilepath = "/lustre03/other/9847zg/dump_files/dump.lammpstrj";
#=dumpfilepath = "/mnt/c/Users/f21157am/OneDrive - The University of Manchester/Codes/MD_to_NMR/dump_files/dump.lammpstrj";=#
dumpfilepath = "C:\\Users\\f21157am\\OneDrive - The University of Manchester\\Codes\\MD_to_NMR\\dump_files\\dump.lammpstrj"
timestep = 1.5; #femptosec

num_lines::Int32 = countlines(dumpfilepath)
natoms = CSV.File(dumpfilepath, skipto=4, limit=1, header=false).Column1[1]
nhydrogens::Int32 = 2 * natoms / 3
totalsteps = (floor(Int, num_lines / (natoms + 9)) ) 

# preallocate  t (time array), and convert to seconds from fs
t = collect(Float32, 0:totalsteps-1) .* (100 * timestep) .* 1e-15

prefactor = (3 / 16) * (μ₀/(4π))^2 * ħ^2 * γ^4

# Intramolecular contributions
F = calculateF(dumpfilepath, "intra", totalsteps)
G = prefactor * mean(ACF.(eachrow(F))) * 1e60 ;
T1r = 1/ (J(ω₀) + 4J(2ω₀))

τ_R = (1 / G[1]) * trapz(t, G) # result is 2.713696293747177e-12
Δω²_R = 3 * prefactor * mean(mean(F.^2, dims=2)) * 10^60 
Δω_R = sqrt(3 * prefactor * mean(mean(F.^2, dims=2)) * 10^60 )/2π
T = 1 / (10 / 3 * Δω²_R * τ_R )

J(ω) = 2 * trapz(t, (G .* cos.(ω .* t)))
ω₀ = γ * 1 # (rad s^-1)

ω₀ = 2 * 2π * 1e6

T1r = 1/ (J(ω₀) + 4J(2ω₀))
T2r = 1/ ((3/2)*J(0) + (5/2)*J(ω₀) + J(2ω₀))


# intermolecular contributions
F = calculateF(dumpfilepath, "inter",totalsteps)
G = prefactor * mean(ACF.(eachrow(F))) * 10^60
τ_T = (1 / G[1]) * trapz(t, G) # result is 5.6951838649564116e-12 
Δω²_T = 3 * prefactor * mean(mean(F.^2, dims=2)) * 10^60

# Relaxation times
T = 1 / (10 / 3 * Δω²_R * τ_R + 10 / 3 * Δω²_T * τ_T)

