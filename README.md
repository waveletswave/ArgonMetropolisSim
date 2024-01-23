# ArgonMetropolisSim

ArgonMetropolisSim is a Python-based tool I developed to practice Monte Carlo simulations as part of my coursework in molecular dynamics simulation. It models the behavior of Argon atoms within a Lennard-Jones potential field using the Metropolis algorithm. This humble project was a significant part of my learning journey, and I hope it can be of use or interest to others studying molecular dynamics simulation or Monte Carlo simulation.

## Features

- Simulates Argon atoms in a cubic box using the Lennard-Jones potential.
- Implements the Metropolis algorithm for Monte Carlo simulations.
- Calculates total potential energy and pressure of the Argon system.
- Allows customization of simulation parameters for various setups.

## Requirements

- Python 3.x
- NumPy library

## Installation

No special installation is needed beyond Python and NumPy. You can clone this repository to your local machine using:

```bash
git clone https://github.com/waveletswave/ArgonMetropolisSim.git
```

## Usage

To run the simulation, simply execute the script with Python:

```bash
python argon_metropolis_sim.py
```

This script performs a Monte Carlo simulation using the Metropolis algorithm and outputs the final potential energy and pressure of the Argon system.

## Configuration

Feel free to adjust the simulation parameters in the script:

- `N`: Number of Argon atoms.
- `L`: Length of the side of the cubic box.
- `T`: Temperature in Kelvin.
- Other physical constants as needed.

## Contributing

Any contributions or suggestions for improvement to ArgonMetropolisSim are warmly welcomed. I believe in collaborative learning and improvement, so please feel free to share your thoughts or best practices.

## License

This project is open-source and available under the [MIT License](LICENSE).

## Acknowledgments

I am grateful for the opportunity to engage in this class learning journey. My sincere appreciation goes to my instructors for their invaluable guidance and to my classmates for the collaborative learning environment we created together. This exercise has been a crucial step in enhancing my understanding of Molecular Dynamics and Monte Carlo Simulations.