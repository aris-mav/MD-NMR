include("./functions.jl")

## Inputs
dumpfilepath = "/lustre03/other/9847zg/dump_files/dump.lammpstrj";
timestep = 1; #femptosec
nblocks = 3;


# Count lines in dump
num_lines::Int32 = countlines(dumpfilepath)

# Read how many atoms are in the simulation
natoms = CSV.File(dumpfilepath, skipto=4, limit=1, header=false).Column1[1]
nhydrogens::Int32 = 2 * natoms / 3

# How many configurations are in dump (how many time steps)
# Modify nsteps so that " nsteps%nblocks=0 " (remove the remainder)
totalsteps = (floor(Int, num_lines / (natoms + 9)) ÷ nblocks) * nblocks

# Preallocate  t (time array)
t = collect(Float32, 0:totalsteps-1) .* (CSV.File(dumpfilepath, skipto=(11 + natoms), limit=1, header=false).Column1 * timestep)

# Double check prefeactors below, they should be working like this but they're not
G = mainloop(dumpfilepath, totalsteps, nblocks, "intra");
Gᵣ = 3 / 16 * (μ₀ / 4π)^2 * γ^4 * ħ^2 * real.(Base.stack(vec(mean(reshape(G, size(G)[1] ÷ nblocks, nblocks), dims=2)))[1, :]);
Δω²ᵣ = 3 * Gᵣ[1];
Jᵣ = 2 * real.(fft(Gᵣ));
τᵣ = real(1 / 2 * (Jᵣ[1] / Gᵣ[1]))
τᵣ = trapz(t[1:300]*1e-15, Gᵣ) / Gᵣ[1]

G = mainloop(dumpfilepath, totalsteps, nblocks, "inter");
Gₜ = 3 / 16 * (μ₀ / 4π)^2 * γ^4 * ħ^2 *real.(Base.stack(vec(mean(reshape(G, size(G)[1] ÷ nblocks, nblocks), dims=2)))[1, :]);
Δω²ₜ = 3 *  Gₜ[1];
Jₜ = 2 * real.(fft(Gₜ));
τₜ = real(1 / 2 * (Jₜ[1] / Gₜ[1]));

T = 1 / (10 / 3 * Δω²ₜ * τₜ + 10 / 3 * Δω²ᵣ * τᵣ)

plot(Gₜ ./ Gₜ[1])


plot(Gᵣ)