#using Plots
using DelimitedFiles
using LinearAlgebra
using FFTW
using Statistics
using Trapz
using StaticArrays

# Define physical constants
const γ = 267.52218744e6; # (rad s^-1 T^-1)
const ħ = 1.054571817e-34;  # (J s)
const μ₀ = 1.25663706212e-6; #N A^-2
const prefactor = (3 / 16) * (μ₀/(4π))^2 * ħ^2 * γ^4

"""
fucntion to calculate spectral density
"""
function J(G::Vector, t::Vector, ω::Real)::Real
    return 2 * trapz(t, (G .* cos.(ω .* t)))
end

"""
function which takes an array of positions and modifies the Hpairs array to
store all the combinations of vectors connecting the hydrogens
(Hpairs is allocated within calculateF function)
"""
function getpairs!(Hpairs::Vector{SVector{3, Float32}},
                   positions::Vector{SVector{3, Float32}},
                   combinations::String)

    counter::Int64 = 1

    if combinations == "all"
        for k in axes(positions, 1)
            for j in axes(positions, 1)
                j > k || continue
                Hpairs[counter] = positions[k] - positions[j]
                counter += 1
            end
        end

    elseif combinations == "intra"
        for i in 1:2:length(positions)
            Hpairs[counter] = positions[i] - positions[i+1]
            counter += 1
        end

    elseif combinations == "inter"
        for i in axes(positions, 1)
            for j in axes(positions, 1)
                (isodd(i) && j > (i + 1)) || (iseven(i) && j > i) || continue
                Hpairs[counter] = positions[i] - positions[j]
                counter += 1
            end
        end

        ## Just a sanity check of how the loop above works, not actual part of the code
        #=for i in 1:10=#
        #=     for j in 1:10=#
        #=        (isodd(i) && j > (i + 1)) || (iseven(i) && j > i) || continue=#
        #=        display("$i with $j")=#
        #=     end=#
        #=end=#
    end

end

"""
applies periodic boundary to the coordinates of the atoms
according to the box size
"""
function periodicboundary!(xyz::Vector{SVector{3, Float32}}, box::MVector{3, Float32})
    displacement_vector::MVector{3, Float32} = @MVector zeros(Float32, 3)
    for i in 1:length(xyz)
        for l in 1:3
            if xyz[i][l] > box[l] / 2
                displacement_vector[l] -= box[l]
            elseif xyz[i][l] < -box[l] / 2
                displacement_vector[l] += box[l]
            end
        end
        xyz[i] += displacement_vector
        displacement_vector .= 0.0
    end
end


"""
converts [x, y, z] to [r, θ, φ]
"""
function cart2sph!(rtp::Vector{SVector{3, Float32}}, xyz::Vector{SVector{3, Float32}})
    for i in eachindex(rtp)
        rtp[i] = @SVector [ norm(xyz[i]), acos(xyz[i][3] / norm(xyz[i])), atan(xyz[i][2], xyz[i][1]) ]
    end
end


"""
F in this case is the quantity (3cos(θ)^2 - 1) / r^3
The output of this function is a matrix which contains this quantity 
for all pairs of hydrogens (columns)
and each time step (rows).
"""
function calculateF(dumpfilepath, contributions::String)

    num_lines::Int = countlines(dumpfilepath)

    natoms::Int = open(dumpfilepath) do io
        for _ in 1:3
            readline(io)
        end
        return parse(Int, readline(io))
    end

    totalsteps::Int = num_lines ÷ (natoms + 9) 

    exclude_atoms = natoms ÷ 2
    natoms -= exclude_atoms
    nhydrogens::Int = 2 * natoms ÷ 3

    npairs::Int64 = 0 

    if contributions == "all"
        npairs += nhydrogens * (nhydrogens - 1) ÷ 2
    elseif contributions == "inter"
        npairs += nhydrogens * (nhydrogens - 1) ÷ 2 - nhydrogens ÷ 2 
    elseif contributions == "intra"
        npairs += nhydrogens ÷ 2
    end

    @show npairs
    @show nhydrogens

    # Initialise arrays
    boxlengths::MVector{3, Float32} = @SVector zeros(Float32, 3)
    positions::Vector{SVector{3, Float32}} = [@SVector zeros(Float32, 3) for _ in 1:nhydrogens]
    Hpairs::Vector{SVector{3, Float32}} = [@SVector zeros(Float32, 3) for _ in 1:npairs]
    F::Matrix{Float32} = zeros(Float32, npairs, totalsteps)
    zvec::SVector{3, Float32} = SA_F32[0.0, 0.0, 1.0]
    vecnorm::Float32 = 1.0

    # Open the dump file
    open(dumpfilepath) do io

        # Loop over time steps
        for s in 1:totalsteps

            # Print progress (optional)
            if isinteractive()
                if s in floor.(Int, collect((totalsteps/10):(totalsteps/10):totalsteps))
                    progresspercent = ceil(s * 100 / totalsteps)
                    display("Calculation progress: $progresspercent %")
                end
            end

            # Skip headers
            for _ in 1:5
                readline(io)
            end

            # Read box bounds
            for i in 1:3
                boxlengths[i] = sum(abs.(parse.(Float32, split(readline(io), ' '))))
            end

            # Skip header
            readline(io)

            # Read positions
            for i in 1:2:nhydrogens
                readline(io) # skip the oxygen
                positions[i] = SVector( parse.(Float32, split(readline(io), ' '))[3:5]...)
                positions[i+1] = SVector( parse.(Float32, split(readline(io), ' '))[3:5]...)
            end

            readuntil(io,"TIME")

            # Do calculations
            getpairs!(Hpairs, positions, contributions)
            periodicboundary!(Hpairs, boxlengths)

            for (i, p) in enumerate(Hpairs)
                vecnorm = norm(p)
                F[i, s] = ( 3 * ( dot(zvec, p)/vecnorm)^2 -1 ) / vecnorm ^3
            end


            # Go to next timestep
        end

        # Close the IO file
    end

    return F

    # Exit the function
end


"""
Radial distribution function
(to see if the model makes sense)


(NEEDS TO BE UPDATED FOR SVECTORS LIKE ABOVE)


"""
function calculate_rdf(dumpfilepath)

    n_lines = countlines(dumpfilepath)

    n_atoms, r_max = open(dumpfilepath) do io
        for _ in 1:3
            readline(io)
        end
        n_atoms = parse(Int, readline(io))
        readline(io)
        r_max = (parse.(Float32, split(readline(io), ' ')))[2] * 0.7 # nk = less than half box length
        return n_atoms, r_max 
    end

    n_oxygen = n_atoms ÷ 3 
    totalsteps = n_lines ÷ (n_atoms + 9)
    n_pairs = n_oxygen * (n_oxygen - 1) ÷ 2

    # Initialise arrays
    boxlengths::Vector{Float32} = zeros(3)
    positions::Vector{Vector{Float32}} = [zeros(3) for _ in 1:n_oxygen]
    o_xyz::Vector{Vector{Float32}} = [[1.0, 1.0, 1.0] for _ in 1:n_pairs]
    o_rtp::Vector{Vector{Float32}} = [[1.0, 1.0, 1.0] for _ in 1:n_pairs]
    h::Vector{Int} = zeros(Int, 1000)

    nk = length(h)
    dr = r_max / nk
    r = [dr * i for i in 0:nk-1]

    # Open the dump file
    open(dumpfilepath) do io
        
        # Loop over time steps
        for s in 1:totalsteps

            # Skip headers
            for _ in 1:5
                readline(io)
            end

            # Read box bounds
            for i in 1:3
                boxlengths[i] .= sum(abs.(parse.(Float32, split(readline(io), ' '))))
            end

            # Skip header
            readline(io)

            # Read positions
            for i in 1:n_oxygen
                positions[i] = parse.(Float32, split(readline(io), ' '))[3:5]
                readline(io) # skip the hydrogen
                readline(io) # skip the hydrogen
            end

            # Do calculations
            getpairs!(o_xyz, positions, "all")
            periodicboundary!(o_xyz, boxlengths)
            cart2sph!(o_rtp, o_xyz)

            for rij in getindex.(o_rtp, 1)

                k = floor(Int, rij/dr) + 1
                if k <= nk
                    h[k] += 2
                end

            end

            # Go to next timestep
        end

        # Close the IO file
    end
    
    ri = r # lower shells
    rj = push!(r[2:end], r[end]+dr) # upper shells
    g = h ./(901 * 501 * (4π /3) * (6.022/180) .* (rj.^3 .- ri.^3)) 

    return r, g  
    # Exit the function
end


"""
extract a time array from the dump file
"""
function time_array(dumpfilepath, timestep)

    num_lines::Int32 = countlines(dumpfilepath)

    natoms = open(dumpfilepath) do io
        for _ in 1:3
            readline(io)
        end
        return parse(Int, readline(io))
    end

    totalsteps = num_lines ÷ (natoms + 9) 

    dump_every_how_many_steps = open(dumpfilepath) do io
        readuntil(io, "TIMESTEP")
        readuntil(io, "TIMESTEP")
        readline(io)
        return  parse(Int,readline(io))
    end
    
    t = collect(Float32, 0:totalsteps-1) .* (dump_every_how_many_steps * timestep) .* 1e-15;

    return t
end


"""
Autocorrelation function
 FFT based recipe:
1.	Pad vector-a by an equal number of zeros. Thus, [1 2 3 0 0 0]
2.	Take the discrete FFT of the new array. Call this F.
3.	Take the conjugate. Call this F*
4.	Compose F times F* (you should do term by term multiplication). Call this Cff.
5.	Take the inverse FFT of Cff, and take only the first 3 terms. Call this ACF.
6.	Now normalize by the vector [3, 2, 1]. That is your answer.
"""
function ACF(v::AbstractVector)

    a = [v ; zeros(length(v))]
    F = fft(a)
    Cff = F .* vec(F')
    ACF = real.(ifft(Cff)[1:length(v)] ./ reverse(collect(1:length(v))))

end

function calculateG(F::Matrix)
    G_ens_av = mean(ACF.(eachrow(F))) # Ensemble average (no prefactors)
    return prefactor * G_ens_av * 1e60 
end

function delta_omega(F::Matrix)
    return 3 * prefactor * mean(mean(F.^2, dims=2)) * 1e60
end
