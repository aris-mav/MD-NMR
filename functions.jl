#using Plots
using DelimitedFiles
using LinearAlgebra
using FFTW
using Statistics
using LsqFit
using Trapz

# Define physical constants
const γ = 267.52218744e6; # (rad s^-1 T^-1)
const ħ = 1.054571817e-34;  # (J s)
const μ₀ = 1.25663706212e-6; #N A^-2

# fucntion to calculate spectral density
J(G, t, ω) = 2 * trapz(t, (G .* cos.(ω .* t)))

function getpairs!(Hpairs::Vector{<:AbstractVector{<:Real}}, positions::Vector{<:AbstractVector{<:Real}}, combinations)

    counter::Int32 = 1

    if combinations == "all"
        for k in axes(positions, 1)
            for j in axes(positions, 1)
                j > k || continue
                Hpairs[counter] .= positions[k] .- positions[j]
                counter += 1
            end
        end

    elseif combinations == "intra"
        for i in 1:2:length(positions)
            Hpairs[counter] .= positions[i] - positions[i+1]
            counter += 1
        end

    elseif combinations == "inter"
        for k in axes(positions, 1)
            for j in axes(positions, 1)
                (isodd(k) && j > (k + 1)) || (iseven(k) && j > k) || continue
                Hpairs[counter] .= positions[k] - positions[j]
                counter += 1
            end
        end

        ## Just a demonstration, not actual part of the code
        #for i in 1:10
        #     for j in 1:10
        #        (isodd(i) && j > (i + 1)) || (iseven(i) && j > i) || continue
        #        display("$i with $j")
        #     end
        #end
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

function calculateF(dumpfilepath,contributions)

    num_lines::Int32 = countlines(dumpfilepath)

    natoms = open(dumpfilepath) do io
        for _ in 1:3
            readline(io)
        end
        return parse(Int, readline(io))
    end

    nhydrogens = 2 * natoms ÷ 3
    totalsteps = num_lines ÷ (natoms + 9) 

    if contributions == "all"
        npairs = nhydrogens * (nhydrogens - 1) ÷ 2
    elseif contributions == "inter"
        npairs = nhydrogens * (nhydrogens - 1) ÷ 2 - nhydrogens ÷ 2 
    elseif contributions == "intra"
        npairs = nhydrogens ÷ 2
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
                progresspercent = ceil(s * 100 / totalsteps)
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
    g = h ./(901* 501* (4π /3) * (6.022/180) .* (rj.^3 .- ri.^3)) 

    return r, g  
    # Exit the function
end


function calc_msd(dumpfilepath, steps)

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


## Autocorrelation function
#= FFT based recipe:
1.	Pad vector-a by an equal number of zeros. Thus, [ 1 2 3 0 0 0]
2.	Take the discrete FFT of the new array. Call this F.
3.	Take the conjugate. Call this F*
4.	Compose F \times F* (you should do term by term multiplication). Call this Cff.
5.	Take the inverse FFT of Cff, and take only the first 3 terms. Call this ACF.
6.	Now normalize by the vector [3, 2, 1]. That is your answer.
=#

function ACF(v::AbstractVector)

    a = [v ; zeros(length(v))]
    F = fft(a)
    Cff = F .* vec(F')
    ACF = real.(ifft(Cff)[1:length(v)] ./ reverse(collect(1:length(v))))

end
