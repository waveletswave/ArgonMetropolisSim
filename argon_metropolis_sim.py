import numpy as np

# Constants
N = 256             # Number of atoms
L = 6.8             # Box side length (assumed to be in reduced units)
epsilon = 1.0       # Depth of potential well (in reduced units)
sigma = 1.0         # Distance at which potential is zero (in reduced units)
T = 100.0           # Temperature in Kelvin
kB = 1.380649e-23   # Boltzmann constant in J/K
epsilon_k = 120.0   # Depth of potential well in K (specific to Argon)
sigma_A = 3.4       # Distance at which potential is zero in Angstroms (for Argon)
cutoff = 2.5        # Cutoff distance for potential calculation (in reduced units)
beta = 1 / (kB * T) # Inverse temperature in reduced units

# Temperature in reduced units specific to Argon
T_star = T / epsilon_k
beta = 1 / T_star

# Function to generate initial positions for atoms
def generate_initial_positions(N, L):
    """
    Generate initial positions for N atoms within a box of size L.
    """
    positions = np.random.rand(N, 3) * L
    return positions

# Generate initial positions for atoms
positions = generate_initial_positions(N, L)

# Ensure positions array is correctly shaped
assert positions.shape == (N, 3), "Position data is not in the correct shape."

# Function to calculate the Lennard-Jones potential
def calculate_LJ_potential(r):
    """
    Calculate the Lennard-Jones potential for a given distance.
    """
    r6 = (sigma / r)**6
    r12 = r6**2
    return 4 * epsilon * (r12 - r6)

# Function to calculate the distance between two atoms
def calculate_distance(atom1, atom2):
    """
    Calculate the distance between two atoms, considering periodic boundary conditions.
    """
    delta = np.abs(atom1 - atom2)
    delta = np.where(delta > 0.5 * L, L - delta, delta)
    return np.sqrt((delta**2).sum())

# Function to calculate the total potential energy of the system
def calculate_total_potential():
    """
    Calculate the total potential energy of the system.
    """
    U_total = 0.0
    for i in range(N):
        for j in range(i + 1, N):
            r = calculate_distance(positions[i], positions[j])
            if r < cutoff:
                U_total += calculate_LJ_potential(r)
    return U_total

# Function to calculate the pressure using the virial theorem in reduced units
def calculate_pressure():
    """
    Calculate the pressure of the system using the virial theorem in reduced units.
    """
    volume = L**3
    ideal_gas_pressure = N * T_star / volume
    virial = 0.0

    for i in range(N):
        for j in range(i + 1, N):
            r_ij_vec = positions[i] - positions[j]
            r_ij_vec = r_ij_vec - np.rint(r_ij_vec / L) * L
            r_ij = np.linalg.norm(r_ij_vec)
            if r_ij < cutoff:
                force_magnitude = 24 * ((2 / r_ij**13) - (1 / r_ij**7))
                virial += force_magnitude * r_ij

    total_pressure = ideal_gas_pressure - virial / (3 * volume)
    return total_pressure

# Function to perform one step of the Metropolis algorithm
def metropolis_step():
    """
    Perform one step of the Metropolis algorithm.
    """
    atom_index = np.random.randint(N)
    displacement = (np.random.rand(3) - 0.5) * 0.1
    old_position = positions[atom_index].copy()
    new_position = (old_position + displacement) % L

    delta_U = 0.0
    for i in range(N):
        if i != atom_index:
            r_old = calculate_distance(old_position, positions[i])
            r_new = calculate_distance(new_position, positions[i])
            if r_old < cutoff:
                delta_U -= calculate_LJ_potential(r_old)
            if r_new < cutoff:
                delta_U += calculate_LJ_potential(r_new)

    if delta_U < 0 or np.random.rand() < np.exp(-beta * delta_U):
        positions[atom_index] = new_position

# Run the Metropolis algorithm for a set number of steps
steps = 10000
for step in range(steps):
    metropolis_step()

# Calculate and display the potential energy and pressure
U_potential = calculate_total_potential()
total_pressure = calculate_pressure()
print("Potential Energy:", U_potential)
print("Pressure:", total_pressure)
