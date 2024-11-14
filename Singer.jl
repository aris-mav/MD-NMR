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
t = collect(Float32, 0:totalsteps-1) .* (100 * timestep) .* 1e-15;

prefactor = (3 / 16) * (μ₀/(4π))^2 * ħ^2 * γ^4

J(ω) = 2 * trapz(t, (G .* cos.(ω .* t)))
ω₀ = γ * 1 # (rad s^-1)

# Intramolecular contributions
F = calculateF(dumpfilepath, "intra", totalsteps);
G_R_ens_av = mean(ACF.(eachrow(F)))
G_R = prefactor * G_R_ens_av * 1e60 ;
τ_R = (1 / G[1]) * trapz(t, G_R) # result is 2.713696293747177e-12
Δω²_R = 3 * prefactor * mean(mean(F.^2, dims=2)) * 1e60
Δω_R = sqrt(Δω²_R)/2π /1000 # KHz
Tr = 1 / (10 / 3 * Δω²_R * τ_R )

# intermolecular contributions
F = calculateF(dumpfilepath, "inter",totalsteps);
G_T_ens_av = mean(ACF.(eachrow(F)))
G_T = prefactor * G_T_ens_av * 1e60 ;
τ_T = (1 / G[1]) * trapz(t, G_T) # result is 2.713696293747177e-12
Δω²_T = 3 * prefactor * mean(mean(F.^2, dims=2)) * 1e60
Δω_T = sqrt(Δω²_T)/2π /1000 # KHz
Tt = 1 / (10 / 3 * Δω²_T * τ_T )

# Relaxation times
T = 1 / ((10 / 3) * Δω²_R * τ_R + (10 / 3) * Δω²_T * τ_T)


T1 = 1/ (J(ω₀) + 4J(2ω₀))
T2 = 1/ ((3/2)*J(0) + (5/2)*J(ω₀) + J(2ω₀))

plot(t, G_T_ens_av, xlabel = "time (s)", label = "Intermolecular correlation function, no prefactors")

