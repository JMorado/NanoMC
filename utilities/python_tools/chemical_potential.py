def chemical_potential(temperature, pressure):
    """
    Temperature is in Kelvin and pressure is in Pascal.
    """
    import numpy as np

    kb = 1.38064852e-23                  # Boltzmann constant (J K^{-1})
    planck = 6.62607004e-34              # Planck constant
    mass = 2.0 * 1.00784 * 1.6737236e-27 # Mass of H2 in Kg
    beta = 1.0 / (kb * T)                # Thermodynamic beta

    # Thermal wavelength
    thermal_wav = planck * planck * beta / (2 * math.pi * mass)
    thermal_wav = thermal_wav ** (1/2.)

    # Calculate chemical potential 
    chemical_pot = np.log(pressure * thermal_wav ** 3 * beta) * T

    print("The chemical potential is: " + str(chemical_pot))
    
    return chemical_pot
