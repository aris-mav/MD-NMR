# %% Import stuff
import numpy as np
from scipy.fft import fft, fftfreq
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import time

BeginScriptTime = time.time()

# %% Define Functions


def cart2sph(xyz):
    r = np.linalg.norm(xyz, axis=0)
    el = np.arccos(xyz[2]/r) 
    az = np.arctan2(xyz[1], xyz[0]) 
    return np.array([r, el, az])
    # return np.array([r, az, el]) 


def F012(rtp):  # input must be [ρ,θ,φ]
    F0 = (1 - 3*np.cos(rtp[1])**2) / rtp[0]**3
    F1 = (np.sin(rtp[1]) * np.cos(rtp[1]) * np.exp(-1j*rtp[2])) / rtp[0]**3
    F2 = (np.sin(rtp[1])**2 * np.exp(-2j*rtp[2])) / rtp[0]**3
    return np.array([F0, F1, F2])

def msdfit(t,D):
    return 6*D*t


# %% Setup parameters
# dumpfilepath = "dump.lammpstrj"
dumpfilepath = "dump_low_density.lammpstrj"

# Count lines in dump
with open(dumpfilepath, 'r') as fp:
    num_lines = sum(1 for line in fp)

# Read how many atoms are in the simulation
natoms = int(np.loadtxt(dumpfilepath, dtype='str', skiprows=3, max_rows=1))

# How many configurations are in dump (how many time steps)
nsteps = int(num_lines/(natoms+9))

del num_lines, fp

# Correlation functions will be averaged between blocks
nGblocks = int(2)


# Modify nsteps so that " nsteps%nGblocks=0 " (remove the remainder)
nsteps = (nsteps//nGblocks)*nGblocks

# Preallocate G (correlation function) and t (time array)
G = np.zeros((3, nsteps//nGblocks, nGblocks), dtype=complex)

effective_steps = nsteps//nGblocks
MSD = np.zeros((effective_steps,nGblocks))
t = np.zeros(nsteps)
timestep = 0.3

# %% Loop to calculate relative positions for all pairs

for i in range(nsteps):

    # Import things from dump: trajectories (r), time steps and box bounds
    r_raw = np.loadtxt(dumpfilepath, skiprows=i*(natoms+9)+9, max_rows=natoms)
    boxbounds = np.loadtxt(dumpfilepath, skiprows=i*(natoms+9)+5, max_rows=3)
    boxlength = boxbounds[:, 1] - boxbounds[:, 0]

    t[i] = np.loadtxt(dumpfilepath, skiprows=i*(natoms+9)+1, max_rows=1)*timestep

    # Remove r for oxygens and M points (we only care about hydrogens)
    r_raw2 = np.delete(r_raw, np.where(r_raw[:, 1] != 2), 0)
    # Keep only coordinates, delete the indices
    r = np.delete(r_raw2, [0, 1], axis=1)

    # if i == 0:
    #     r0 = r
    #     r_prev = r
    # else:

    #     for m in range(3):
    #         r[:, m][np.argwhere(
    #             (r[:, m]-r_prev[:, m]) < -boxlength[m]/2)] += boxlength[m]

    #         r[:, m][np.argwhere(
    #             (r[:, m]-r_prev[:, m]) > boxlength[m]/2)] -= boxlength[m]

    #     MSD[i] = np.mean(np.sum((r-r0)**2, axis=1))
    #     r_prev = r

    # Make empty array for all the possible pairs of hydrogens
    nhydrogens = np.shape(r)[0]
    npairs = int(nhydrogens*(nhydrogens-1)/2)
    rpairs = np.zeros((npairs, 3))

    rpairsintra = np.zeros((int(nhydrogens/2), 3))
    ni = 0
    n = 0

    for k in range(nhydrogens):
        rpairs[n: n + nhydrogens-1-k, :] = r[k+1: , :] - r[k, :]

        # Collect intramolecular pairs (optional)
        if (k % 2) == 0:
            rpairsintra[ni] = rpairs[n]
            ni += 1

        n += nhydrogens-1-k

    for m in range(3):
        rpairs[:, m][np.argwhere(
            rpairs[:, m] < -boxlength[m]/2)] += boxlength[m]
        rpairs[:, m][np.argwhere(
            rpairs[:, m] > boxlength[m]/2)] -= boxlength[m]

    # The "blockindex" below will be an integer only at the 1st time step of each block
    # It is used to index the 3rd dimension of the G array below
    blockindex = i/(nsteps/nGblocks)

    # Compute F values for all pairs of spins at time step i
    F = F012(cart2sph(rpairs.T))

    # Compute the ensebmble average for G at time step i


    effective_index = i-int(blockindex)*effective_steps
    if blockindex % 1 == 0:   # if it is the 1st time step of the block
        F0 = F
        G[:,effective_index , int(blockindex)
          ] = np.mean(F0 * np.conjugate(F), axis=1)

        print("Calculating block " + str(int(blockindex)+1) +
              " out of " + str(nGblocks))
        
        r0 = r
        r_prev = r

    else:                     # if it's not the 1st time step of the block
        G[:, effective_index, int(blockindex)
          ] = np.mean(F0 * np.conjugate(F), axis=1)
        
        #MSD remove periodic boundary
        for m in range(3):
            r[:, m][np.argwhere(
                (r[:, m]-r_prev[:, m]) < -boxlength[m]/2)] += boxlength[m]

            r[:, m][np.argwhere(
                (r[:, m]-r_prev[:, m]) > boxlength[m]/2)] -= boxlength[m]

        # MSD compute 
        MSD[effective_index,int(blockindex)] = np.mean(np.sum((r-r0)**2, axis=1))
        r_prev = r

# %% Compute the spectral density function

# Do the fourier transform of the correlation function, and save half of it as J
J = 2*fft(np.mean(G, axis=2))[:, range(int(np.shape(G)[1])//2)]

# Covnert to real numbers
J = np.real(J)

# Compute frequency bins (cycles per unit of the sample spacing), and cut it in half
# (*1e15 converts to Hz)
omega = fftfreq(nsteps, t[1]-t[0])[range(np.shape(J)[1])] * 1e15


# Set constants
gamma = 267.52218744e6  # (rad s^-1 T^-1)
hbar = 1.054571817e-34  # (J s)


# Calculate T1 and T2
T1 = 1/((9/8) * gamma**4 * hbar**2 * (J[1, 0]+J[2, 0]))
T2 = 1/((9/24)* gamma**4 * hbar**2 * (J[0, 0]+J[1, 0]+J[2, 0]))

# T1s = 1./((10**-7)**2* (6/20 * hbar**2 *gamma**4) * (J[1, 0]  + 4*J[2,0] ) )
# T2s = 1./((10**-7)**2* (3/20 * hbar**2 *gamma**4) * (3*J[0,0] + 5*J[1,0] + 2*J[2,0]) )


# %% Plots
fig1, ax1 = plt.subplots(1, 2, figsize=(12, 5))

for i in range(3):
    ax1[0].plot(t[0:np.shape(G)[1]]/1e3, np.abs(np.mean(G, axis=2)[i, :]),
                # (/1e3 converts fempto to picosec)
                label="$G_{" + str(i)*2 + "}$")

ax1[0].legend(loc="upper right")
ax1[0].set_title('Correlation functions')
ax1[0].set(xlabel="t-τ (ps)",
           ylabel=r"$\langle F_q(τ)$ $F_q^*(τ+t) \rangle$  (unit length$^{-6})$")



for i in range(3):
    ax1[1].plot(omega, np.abs(J[i, :]), label="$J_{" + str(i)*2 + "}$")

ax1[1].legend(loc="upper right")
ax1[1].set_title("Spectral densities")
ax1[1].set(xlabel="Hz", ylabel=r"$ |\mathcal{F}(G)|$")


plt.savefig("Correlations.png")


fig2, ax2 = plt.subplots()
ax2.scatter(t[range(effective_steps)], np.mean(MSD,axis=1),marker=',',s=0.6)
ax2.set(xlabel="time (femptoseconds)", ylabel="msd (angstrom^2)")

D,_ = curve_fit(msdfit,t[range(effective_steps)],np.mean(MSD,axis=1)) #converting to m2s-1
ax2.plot(t[range(effective_steps)],msdfit(t[range(effective_steps)],D),color='red')

D = D * 1e-5 #convert from Ang2fs-1 to m2s-1

# Print timer
EndScriptTime = time.time()
totaltime = EndScriptTime-BeginScriptTime
print("Runtime: " + str(round(EndScriptTime-BeginScriptTime))+" seconds.")

print("T1: " + str(T1)+" seconds.")
print("T2: " + str(T2)+" seconds.")
print("D: " + str(D)+" m2s-1.")

