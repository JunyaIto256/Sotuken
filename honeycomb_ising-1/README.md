# Honeycomb Ising Model Simulation

This project implements a 2D Ising model on a honeycomb lattice using the Metropolis Monte Carlo method. The simulation allows for the study of magnetic properties and phase transitions in statistical mechanics.

## Project Structure

- **src/HoneycombIsing.jl**: Defines the main module for the honeycomb lattice Ising model, including functions for spin initialization, energy calculations, and spin updates using the Metropolis algorithm.
  
- **src/lattice.jl**: Contains functions for defining the structure of the honeycomb lattice, including lattice generation and obtaining neighboring spins.
  
- **src/metropolis.jl**: Implements functions related to the Metropolis algorithm, including spin updates and energy difference calculations.
  
- **src/observables.jl**: Provides functions for calculating physical observables such as average magnetization and spin correlations.
  
- **scripts/run_simulation.jl**: A script to execute the simulation, setting necessary parameters and starting the simulation process.
  
- **test/runtests.jl**: A script for running unit tests to verify the correctness of functions in each module.
  
- **Project.toml**: Defines the project's dependencies and settings, listing the packages and their versions used in the project.
  
- **.gitignore**: Specifies files and directories to be ignored by Git.

## Usage

To run the simulation, execute the `run_simulation.jl` script located in the `scripts` directory. You can modify the parameters within the script to customize the simulation settings.

## Requirements

Ensure you have Julia installed along with the necessary packages specified in the `Project.toml` file. 

## License

This project is licensed under the MIT License.