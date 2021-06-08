def fugacity_coeff(temperature, pressure):
    """
    Fugacity coefficients for hydrogen gas between 0C and 1000C, for pressures up to 3000 atm.
    Temperature is in Kelvin and pressure is in atm.

    Note than in the original paper, c3 coefficient is wrong.
    Originally, it is written as:

    c3 = 300.0*np.exp(-0.11901*temperature - 5.941)

    But for the results of this equation to match the table in the paper one has to use:
   
    c3 = 300.0*np.exp(-0.011901*temperature - 5.941)

    Reference:
    Shaw, H. R.; Wones, D. R. 
    Fugacity Coefficients for Hydrogen Gas between O Degrees and 1000 Degrees C, for Pressures to 3000 Atm. 
    American Journal of Science 1964, 262 (7), 918–929. 
    https://doi.org/10.2475/ajs.262.7.918.
    """
    import numpy as np

    if pressure > 3000:
        print("Pressure {} is greater than 3000 atm.".format(pressure))
        exit()

    # Define coefficients
    c1 = np.exp(-3.8402*temperature**(1.0/8.0) + 0.5410)
    c2 = np.exp(-0.1263*temperature**(1.0/2.0) - 15.980)
    c3 = 300.0*np.exp(-0.011901*temperature - 5.941)

    # Calculate fugacity coefficient
    fugacity = c1*pressure - c2*pressure*pressure + c3*(np.exp(-pressure/300.0)-1.0)

    return np.exp(fugacity)



