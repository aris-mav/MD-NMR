using Plots
using LinearAlgebra
using FFTW
using Statistics
using CSV
using LsqFit

# Define physical constants
const γ = 267.52218744e6; # (rad s^-1 T^-1)
const ħ = 1.054571817e-34;  # (J s)
const μ₀ = 1.25663706212e-6; #N A^-2

# Define the dump file before running script
# if !@isdefined(dumpfilepath)

#     display("Please provide a dump file path below:")
#     global dumpfilepath = readline()
# end

dumpfilepath = "/lustre03/other/9847zg/dump_files/dump.lammpstrj"

nblocks = 3
timestep = 3 #femptosec

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


function getpairs!(Hpairs::Vector{<:AbstractVector{<:Real}}, positions::Vector{<:AbstractVector{<:Real}}, combinations)
    counter::Int32 = 1
    if combinations == "all"
        for k in axes(positions, 1)
            for j in axes(positions, 1)
                j > k || continue
                Hpairs[counter] .= @views (positions[k] .- positions[j])
                counter += 1
            end
        end

    elseif combinations == "intra"
        for i in 1:2:length(positions)
            Hpairs[counter] = positions[i] - positions[i+1]
            counter += 1
        end

    elseif combinations == "inter"
        for k in axes(positions, 1)
            for j in axes(positions, 1)
                (isodd(k) && j > (k + 1)) || (iseven(k) && j > k) || continue
                Hpairs[counter] .= @views (positions[k] .- positions[j])
                counter += 1
            end
        end

    end
end


function periodicboundary!(xyz::Vector{<:AbstractVector{<:Real}}, box::Vector{<:Real})
    for i in 1:length(xyz)
        for l in 1:3
            if xyz[i][l] > box[l] / 2
                xyz[i][l] -= box[l]
            elseif xyz[i][l] < -box[l] / 2
                xyz[i][l] += box[l]
            end
        end
    end
end



function cart2sph!(rtp::Vector{<:AbstractVector{<:Real}}, xyz::Vector{<:AbstractVector{<:Real}})
    for i in eachindex(rtp)
        rtp[i][1] = norm(xyz[i])
        rtp[i][2] = acos(xyz[i][3] / norm(xyz[i]))
        rtp[i][3] = atan(xyz[i][2], xyz[i][1])
    end
end



function F012!(F::Vector{<:AbstractVector{<:Complex}}, rtp::Vector{<:AbstractVector{<:Real}})
    for i in 1:length(rtp)
        F[i][1] = complex((3 * cos(rtp[i][2])^2 - 1) / (rtp[i][1]^3))
        F[i][2] = (sin(rtp[i][2]) * cos(rtp[i][2]) * exp(-1im * rtp[i][3])) / rtp[i][1]^3
        F[i][3] = (sin(rtp[i][2])^2 * exp(-2im * rtp[i][3])) / rtp[i][1]^3
    end
end




function mainloop(steps, blocks, combinations)

    if combinations == "all"
        npairs = floor(Int, nhydrogens * (nhydrogens - 1) / 2)
    elseif combinations == "inter"
        npairs = floor(Int, nhydrogens * (nhydrogens - 1) / 2) - nhydrogens / 2
    elseif combinations == "intra"
        npairs = nhydrogens / 2
    end

    # Initialise arrays
    boxlengths::Vector{Float32} = zeros(3)
    positions::Vector{Vector{Float32}} = [zeros(3) for _ in 1:nhydrogens]
    Hpairs::Vector{Vector{Float32}} = [[1.0, 1.0, 1.0] for _ in 1:npairs]
    rtp::Vector{Vector{Float32}} = [[1.0, 1.0, 1.0] for _ in 1:npairs]
    F::Vector{Vector{ComplexF32}} = [ones(ComplexF32, 3) for _ in 1:npairs]
    F0::Vector{Vector{ComplexF32}} = [ones(ComplexF32, 3) for _ in 1:npairs]
    G::Vector{Vector{ComplexF32}} = [ones(ComplexF32, 3) for _ in 1:steps]

    # Open the dump file
    open(dumpfilepath) do io

        # Loop over time steps
        for s in 1:steps

            # Print progress (optional)

            if s in floor.(Int, collect((steps/10):(steps/10):steps))
                progresspercent = round(s * 100 / steps, digits=2)
                display("Calculation progress: $progresspercent %")
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
                positions[i] = parse.(Float32, split(readline(io), ' '))[3:5]
                positions[i+1] = parse.(Float32, split(readline(io), ' '))[3:5]
            end

            # Do calculations
            getpairs!(Hpairs, positions, combinations)
            periodicboundary!(Hpairs, boxlengths)
            cart2sph!(rtp, Hpairs)
            F012!(F, rtp)

            if isinteger((s - 1) / (steps / blocks)) # this number will only be an integer at the beginning of each block
                F0[:] = deepcopy(F)
            end

            G[s][:] .= @views vec(mean(Base.stack(F0) .* Base.stack(conj.(F)), dims=2))

            # Go to next timestep
        end

        # Close the IO file
    end

    return G

    # Exit the function
end


function calculateF(contributions)

    if contributions == "all"
        npairs = floor(Int, nhydrogens * (nhydrogens - 1) / 2)
    elseif contributions == "inter"
        npairs = floor(Int, nhydrogens * (nhydrogens - 1) / 2) - nhydrogens / 2
    elseif contributions == "intra"
        npairs = floor(Int, nhydrogens / 2)
    end

    # Initialise arrays
    boxlengths::Vector{Float32} = zeros(3)
    positions::Vector{Vector{Float32}} = [zeros(3) for _ in 1:nhydrogens]
    Hpairs::Vector{Vector{Float32}} = [[1.0, 1.0, 1.0] for _ in 1:npairs]
    rtp::Vector{Vector{Float32}} = [[1.0, 1.0, 1.0] for _ in 1:npairs]
    F::Matrix{Float32} = zeros(npairs, totalsteps)

    # Open the dump file
    open(dumpfilepath) do io

        # Loop over time steps
        for s in 1:totalsteps

            # Print progress (optional)

            if s in floor.(Int, collect((totalsteps/10):(totalsteps/10):totalsteps))
                progresspercent = round(s * 100 / totalsteps, digits=2)
                display("Calculation progress: $progresspercent %")
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
                positions[i] = parse.(Float32, split(readline(io), ' '))[3:5]
                positions[i+1] = parse.(Float32, split(readline(io), ' '))[3:5]
            end

            # Do calculations
            getpairs!(Hpairs, positions, contributions)
            periodicboundary!(Hpairs, boxlengths)
            cart2sph!(rtp, Hpairs)


            for (i, j) in enumerate(rtp)
                F[i, s] = (3 * cos(j[2])^2 - 1) / j[1]^3

            end

            # Go to next timestep
        end

        # Close the IO file
    end

    return F

    # Exit the function
end


function calc_msd(steps)

    MSD::Vector{Float32} = zeros(steps)
    boxlengths::Vector{Float32} = zeros(3)
    r::Vector{Vector{Float32}} = [zeros(3) for _ in 1:nhydrogens]
    r0 = deepcopy(r)
    r_prev = deepcopy(r)

    open(dumpfilepath) do io

        # Loop over time steps
        for s in 1:steps

            # Print progress (optional)
            if s in floor.(Int, collect((steps/10):(steps/10):steps))
                progresspercent = round(s * 100 / steps, digits=2)
                display("Calculation progress: $progresspercent %")
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
                r[i] = parse.(Float32, split(readline(io), ' '))[3:5]
                r[i+1] = parse.(Float32, split(readline(io), ' '))[3:5]
            end


            if s == 1
                r0 = deepcopy(r)
                r_prev = deepcopy(r)

            else

                for i in 1:nhydrogens
                    for j in 1:3

                        if (r[i][j] - r_prev[i][j]) < -boxlengths[j] / 2
                            r[i][j] += boxlengths[j]
                        elseif (r[i][j] - r_prev[i][j]) > boxlengths[j] / 2
                            r[i][j] -= boxlengths[j]
                        end
                    end


                end
                MSD[s] = mean((norm.(r - r0)) .^ 2)
                r_prev = deepcopy(r)

            end

            # Go to next timestep
        end

        # Close the IO file
    end


    t = collect(range(0, (steps - 1) * 100 * 1, steps))

    eq(x, D) = 6 * D .* x
    fit = curve_fit(eq, t, MSD, [2e-9])

    p = scatter(t, MSD)
    plot!(t, eq(t, coef(fit)))
    display(p)

    D = coef(fit)[1] * 1e-5 #convert to m^2s^-1

    println("The diffusion coefficient is $D m^2s^-1")
    return D
end


function plotG(Gmean)
    a = plot(t[1:size(Gmean)[2]], real.(Gmean[1, :]), label="G0")
    a = plot!(t[1:size(Gmean)[2]], real.(Gmean[2, :]), label="G1")
    a = plot!(t[1:size(Gmean)[2]], real.(Gmean[3, :]), label="G2")
    # plot!(background_color=:transparent)
    a = plot!(xaxis="time (fs)", yaxis="G (Angstrom^-6)")
    display(a)
end


## Method for selecting inta or inter hydrogens exclusively
## Just a demonstration, not actual part of the code
#for i in 1:10
#     for j in 1:10
#        (isodd(i) && j > (i + 1)) || (iseven(i) && j > i) || continue
#        display("$i with $j")
#     end
#end
#
#for i in 1:2:10
#    display("$i with $(i+1)")
#end
