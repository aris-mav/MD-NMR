include("./functions.jl")
unicodeplots()

## Inputs
dumpfilepath = "/lustre03/other/9847zg/dump_files/dump.lammpstrj";
timestep = 1.5; #femptosec
nblocks = 1;
num_lines::Int32 = countlines(dumpfilepath)
natoms = CSV.File(dumpfilepath, skipto=4, limit=1, header=false).Column1[1]
nhydrogens::Int32 = 2 * natoms / 3
totalsteps = (floor(Int, num_lines / (natoms + 9)) ÷ nblocks) * nblocks

# preallocate  t (time array)
t = collect(float32, 0:totalsteps-1) .* (csv.file(dumpfilepath, skipto=(11 + natoms), limit=1, header=false).column1 * timestep) .* 1e-15

# 
ω₀ = γ * 1; # (rad s^-1)
J(ω) =  2 * trapz(t, (G .* cos.(ω .*t)))
prefactor = (3/16) * (μ₀/4π)^2 * ħ^2 * γ^2

# Intramolecular contributions

F = calculateF(dumpfilepath,"intra")
G = prefactor * mean(ACF2.(eachrow(F))) 

τ = (1 / G[1]) * trapz(t,G)

Δω² = 3* G[1]

plot( abs.(mean(ACF.(eachrow(F))))  )
plot( abs.(mean(ACF2.(eachrow(F))))  )
plot( abs.(mean(ACF3.(eachrow(F))))  )
plot( abs.(mean(ACF4.(eachrow(F))))  )
plot( abs.(mean(autocor.(eachrow([F zeros(size(F))])))))
