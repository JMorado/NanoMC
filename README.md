     _   _                   __  __  _____ 
    | \ | |                 |  \/  |/ ____|
    |  \| | __ _ _ __   ___ | \  / | |     
    | . ` |/ _` | '_ \ / _ \| |\/| | |     
    | |\  | (_| | | | | (_) | |  | | |____ 
    |_| \_|\__,_|_| |_|\___/|_|  |_|\_____|
                                            
    
    
NanoMC is a software that simulates Lennard-Jones fluids inside nanotubes.

                                    
### Features:

Monte Carlo:

- Canonical Ensemble (NVT)
- Grand Canonical Ensemble (uVT)
    
### Installation

To compile the code simply type the following line at the main NanoMC directory
    
    make 
    
To clean up the code use 
 
    make clean
    
 or
    
    make veryclean

### Code structure:

- __nanomc_nvt.f90__: main program for NVT MC simulation
- __nanomc_uvt.f90__: main program for uVT MC simulation
- __monte_carlo_module.f90__: subroutines that perform core parts of uVT and NVT MC simulations
- __simulation_module.f90__: simulation object module
- __pbc_module.f90__: periodic boundary conditions module
- __displacement.f90__: MC displacements module
- __time_module.f90__: time-dependent MC module
- __cell.f90__: simulation cell object module
- __io_module.f90__: subroutines that perform input/output operations
- __energy_module.f90__: subroutines and functions to calculate energies
- __std_output_module.f90__: subroutines to print standard output 
- __constants_module.f90__: definition of general-purpose constants

### References

To be included soon.



